"""B03 Phase 5 (step B): E2 cross-platform coverage — Proseg replication on the crop.

Faithful port of b03_phase2_analysis.py (fixture-validated tissue-1 machinery) with:
  - crop paths (xenium_lung/crop/data, xenium_lung/crop/proseg_out)
  - 18 G1-testable pairs / 54 rows; A3-frozen TYPE_MAP/RECV_MAP
  - certified set from scout_final_lung.csv (cert_pos_final / cert_neg_final)
All schema assertions from tissue 1 retained (amendment record naming, noise-slot handling,
gene-alignment proof, cell-id/mtx mapping, non-evaluable cap).

Endpoint (prereg E2): fraction of certified rows with T(A_proseg) inside the
certified interval. Gate >= 90% PASS, 75-90% PARTIAL, < 75% FAIL.
Secondary: sign agreement on every certified row.
Output: xenium_lung/crop/data/phase5_e2_coverage.csv + phase5_e2_summary.json
"""
import gzip, json, os, sys, time
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_scout import EPS, t_log2

D = "xenium_lung/crop/data"
PRO = "xenium_lung/crop/proseg_out"
PAIRS = [("CD274","PDCD1"),("CXCL12","CXCR4"),("PECAM1","KDR"),
 ("ERBB2","EGFR"),("ERBB2","PDCD1"),("ESR1","PGR"),("PGR","ESR1"),
 ("CD274","CTLA4"),("CD274","CD8A"),("CXCL12","CCR7"),("PECAM1","PDCD1"),
 ("PECAM1","CTLA4"),("CD163","CCR7"),("CD163","CXCR4"),("CDH1","EGFR"),
 ("CDH1","ERBB2"),("MS4A1","CD274"),("CD68","CD274")]
TYPE_MAP = {
 "CD274": "Myeloid", "CXCL12": "Fibroblast", "PECAM1": "Endothelial",
 "CD163": "Myeloid", "ERBB2": "Epithelial", "ESR1": "Epithelial",
 "PGR": "Epithelial", "CDH1": "Epithelial", "MS4A1": "B_cells", "CD68": "Myeloid",
}
RECV_MAP = {
 "PDCD1": "T_cells", "CXCR4": "Myeloid", "CCR7": "Myeloid",
 "KDR": "Endothelial", "EGFR": "Epithelial", "PGR": "Epithelial",
 "ESR1": "Epithelial", "CTLA4": "T_cells", "CD8A": "T_cells",
 "ERBB2": "Epithelial", "CD274": "Epithelial",
}
t0 = time.time()
tested = sorted({g for p in PAIRS for g in p})
gidx = {g: i for i, g in enumerate(tested)}

# ---------- frozen vendor-side artifacts ----------
vm = np.load(f"{D}/vendor_markers.npz", allow_pickle=True)
centroids = vm["centroids"].astype(np.float64)
type_names = [str(t) for t in vm["type_names"]]
marker_genes = [str(g) for g in vm["marker_genes"]]
scale = float(vm["scale"])
vendor_xy = vm["vendor_centroids_xy"]
regions_vendor = vm["regions_row"].astype(int)

# ---------- certified intervals + certified set ----------
b = pd.read_csv(f"{D}/scout_bounds.csv")
assert len(b) == 54, len(b)
fin = pd.read_csv(f"{D}/scout_final_lung.csv")
assert len(fin) == 54

# ---------- proseg outputs ----------
with gzip.open(f"{PRO}/counts", "rb") as f:
    Coo = mmread(f).tocoo()
cells_meta_p = pd.read_parquet(f"{PRO}/cells")
genes_p = pd.read_parquet(f"{PRO}/genes")
print(f"[{time.time()-t0:.0f}s] mtx {Coo.shape} nnz={Coo.nnz:,}; cells "
      f"{cells_meta_p.shape}; genes {genes_p.shape}", flush=True)
assert Coo.shape[0] == len(cells_meta_p)
assert Coo.shape[1] in (len(genes_p), len(genes_p) - 1)
gene_names = genes_p["gene"].astype(str).to_numpy()
assert len(set(gene_names)) == len(gene_names)

pc = cells_meta_p
is_bg = pc["original_cell_id"].astype(str) == "0"
n_bg_cells = int(is_bg.sum())
if n_bg_cells:
    keep = (~is_bg).to_numpy().nonzero()[0]
    Coo = Coo.tocsr()[keep].tocoo()
    pc = pc[~is_bg].reset_index(drop=True)
n_pc = len(pc)
print(f"[{time.time()-t0:.0f}s] background rows dropped: {n_bg_cells}; universe {n_pc}",
      flush=True)
ids_kept = pc["cell"].to_numpy(np.int64)
assert (np.diff(ids_kept) > 0).all()
assert (n_bg_cells > 0) or (ids_kept == np.arange(n_pc)).all()

mtx_n_genes = Coo.shape[1]
dropped_noise_gene = None
if mtx_n_genes == len(gene_names):
    gene_names_eff = gene_names
elif mtx_n_genes == len(gene_names) - 1:
    dropped_noise_gene = str(gene_names[-1])
    assert dropped_noise_gene not in set(marker_genes) | set(tested)
    gene_names_eff = gene_names[:-1]
else:
    raise AssertionError("mtx gene cols match neither schema")

need = marker_genes + [g for g in tested if g not in set(marker_genes)]
gene_pos = {g: i for i, g in enumerate(gene_names_eff)}
missing = [g for g in need if g not in gene_pos]
assert not missing, f"genes missing from proseg output: {missing[:5]}"
gene_totals = genes_p.set_index(genes_p["gene"].astype(str))["total_count"]
cols = np.array([gene_pos[g] for g in need])
Sub = Coo.tocsr()[:, cols]
col_sums = np.asarray(Sub.sum(axis=0)).ravel()
viol = []
for k, g in enumerate(need):
    tot = float(gene_totals[g])
    if col_sums[k] > tot * 1.02 + 1e-6:
        viol.append((g, float(col_sums[k]), tot))
assert not viol, f"GENE MISALIGNMENT: {viol[:5]}"
print(f"[{time.time()-t0:.0f}s] alignment proof: {len(need)} genes, 0 violations",
      flush=True)

Sub = Sub.tocoo()
Mall = np.zeros((n_pc, len(need)), dtype=np.int32)
np.add.at(Mall, (Sub.row.astype(np.int64), Sub.col.astype(np.int64)),
          np.rint(Sub.data).astype(np.int64))
need_pos = {g: i for i, g in enumerate(need)}
M = Mall[:, [need_pos[g] for g in marker_genes]].astype(np.int32)
Gm = Mall[:, [need_pos[g] for g in tested]].astype(np.int64)
assert (Gm.sum(axis=0) > 0).all()

# ---------- label transfer (identical transform) ----------
Ctot = M.sum(axis=1).astype(np.int64)
has_marker = Ctot > 0
Xp = np.log1p((M / np.maximum(Ctot, 1)[:, None]) * scale).astype(np.float64)
Xn = Xp / np.maximum(np.linalg.norm(Xp, axis=1, keepdims=True), 1e-12)
Cn = centroids / np.maximum(np.linalg.norm(centroids, axis=1, keepdims=True), 1e-12)
labels = np.full(n_pc, "Unknown", dtype=object)
sim = Xn[has_marker] @ Cn.T
labels[has_marker] = np.array(type_names, dtype=object)[sim.argmax(axis=1)]
print(f"[{time.time()-t0:.0f}s] labels: {pd.Series(labels).value_counts().to_dict()}",
      flush=True)

# ---------- region mapping (nearest vendor centroid, tissue-1 procedure) ----------
pc_xy = pc[["centroid_x", "centroid_y"]].to_numpy(np.float64)
tree = cKDTree(vendor_xy)
dist, nn = tree.query(pc_xy, k=1)
regions = regions_vendor[nn]
print(f"[{time.time()-t0:.0f}s] regions: {np.bincount(regions).tolist()}; "
      f"nn-dist median {np.median(dist):.2f} um", flush=True)
assert np.median(dist) < 15.0

# ---------- T(A_proseg) for the 54 rows ----------
types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
n_types = len(types_needed)
SHIFT = 1
type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
n_slots = n_types + SHIFT
cell_slot = regions * n_slots + np.array(
    [type_shift.get(l, 0) for l in labels], dtype=np.int64)

rows = []
for (Lg, Rg) in PAIRS:
    S, R = TYPE_MAP[Lg], RECV_MAP[Rg]
    for region_id in range(3):
        sS = cell_slot == region_id * n_slots + type_shift[S]
        sR = cell_slot == region_id * n_slots + type_shift[R]
        nL, nR = int(sS.sum()), int(sR.sum())
        if nL == 0 or nR == 0:
            rows.append(dict(pair=f"{Lg}->{Rg}", region=region_id, nL=nL, nR=nR,
                             N_L=0, N_R=0, T=None))
            continue
        N_L = int(Gm[sS, gidx[Lg]].sum())
        N_R = int(Gm[sR, gidx[Rg]].sum())
        rows.append(dict(pair=f"{Lg}->{Rg}", region=region_id, nL=nL, nR=nR,
                         N_L=N_L, N_R=N_R, T=t_log2(N_L, N_R, nL, nR)))
P = pd.DataFrame(rows)
assert len(P) == 54
n_excl = int(P["T"].isna().sum())
n_eval = 54 - n_excl
print(f"[{time.time()-t0:.0f}s] evaluable rows: {n_eval}/54 (excluded {n_excl})",
      flush=True)
assert n_excl <= 5, f"{n_excl}/54 non-evaluable (>10%) - ABORT"

# ---------- endpoints ----------
PE = P.dropna(subset=["T"])
out = b[["pair","region","control","T0","T_lo","T_hi"]].merge(
    PE[["pair","region","T","nL","nR"]], on=["pair","region"], how="left",
    validate="1:1")
assert len(out) == 54
assert ((out.T0 + 1e-9 >= out.T_lo) & (out.T0 - 1e-9 <= out.T_hi)).all()
out["inside"] = (out["T"] >= out.T_lo - 1e-9) & (out["T"] <= out.T_hi + 1e-9)
fin2 = fin[["pair","region","cert_pos_final","cert_neg_final"]]
out = out.merge(fin2, on=["pair","region"], validate="1:1")
cert_mask = out.cert_pos_final | out.cert_neg_final
cert_eval = cert_mask & out["T"].notna()
n_in = int(out.loc[cert_eval, "inside"].sum())
frac = n_in / int(cert_eval.sum())
sign_ok = int((((out["T"] > 0) == out.cert_pos_final)[cert_eval]).sum())

summary = dict(
    prereg="B03_PHASE5_PREREG.md v1.0 + ADDENDUM v1.0 (A2 crop)",
    n_proseg_cells=int(n_pc), n_vendor_cells=int(len(b)) ,
    cells_with_marker_signal=int(has_marker.sum()),
    noise_gene_slot=dropped_noise_gene,
    primary=dict(n_inside=n_in, n_certified_evaluable=int(cert_eval.sum()),
                 fraction=frac, gate=("PASS" if frac >= 0.90 else
                                      "PARTIAL" if frac >= 0.75 else "FAIL")),
    certified_rows=int(cert_mask.sum()),
    sign_agreement_on_certified=f"{sign_ok}/{int(cert_eval.sum())}",
)
out.to_csv(f"{D}/phase5_e2_coverage.csv", index=False)
with open(f"{D}/phase5_e2_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2), flush=True)
print("PHASE5 E2 COMPLETE", flush=True)
