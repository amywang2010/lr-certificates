"""B03 Phase 7 (step B): E2 cross-segmentation coverage — Proseg replication on s6.

Faithful port of b03_phase5_e2.py (fixture-validated machinery) with:
  - crop paths (workspace xenium_breast_s6/crop/data, xenium_breast_s6/crop/proseg_out)
  - 10 vendor-evaluable pairs / 30 rows; A3-frozen TYPE_MAP/RECV_MAP (imported)
  - marker centroids from vendor_markers.npz (A1-h phase2-space targets)
  - Endpoints: the certified-set primary endpoint is UNDEFINED at 0
    certified rows (reported as such); registered replacements are (a) interval
    coverage on all 30 vendor-evaluable rows with the standing T0-in-interval
    assertion, (b) sign agreement on all 30 rows.
All schema assertions from tissues 1/2 retained (type-naming, noise-slot
handling, gene-alignment proof, cell-id/mtx mapping).

Output: the s6 crop data dir's phase7_e2_coverage.csv + phase7_e2_summary.json
"""
import gzip, json, os, sys, time
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_scout import EPS, t_log2
from b03_lung_config import PAIRS, TYPE_MAP, RECV_MAP

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
D = os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "data")
PRO = os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "proseg_out")
t0 = time.time()
tested = sorted({g for p in PAIRS for g in p})
gidx = {g: i for i, g in enumerate(tested)}

# ---------- frozen vendor-side artifacts ----------
vm = np.load(f"{D}/vendor_markers.npz", allow_pickle=True)
centroids = vm["centroids"].astype(np.float64)
type_names = [str(t) for t in vm["type_names"]]
marker_genes = [str(g) for g in vm["marker_genes"]]
scale = float(vm["scale"])
assert scale > 0, "A1-h violated: scale == 0"
vendor_xy = vm["vendor_centroids_xy"]
regions_vendor = vm["regions_row"].astype(int)

# ---------- certified intervals ----------
b = pd.read_csv(f"{D}/scout_bounds.csv")
assert len(b) == 30, len(b)

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
# A1-j: proseg dropped 4 ultra-rare custom probes (18/4.91M tx, 0.0004%) and the
# A1-f fill added 3 codeword pseudo-genes. The marker space is intersected with
# the proseg output (only A1-j probes may drop, counted); tested-gene completeness
# and the alignment proof below remain binding.
missing = [g for g in need if g not in gene_pos]
assert set(missing) <= {"Human_CTNNB1_S33C_ALT:G", "Human_FGFR3_R248H_ALT:A",
                        "Human_PTEN_G129E_WT:A", "SRY"}, \
    f"unexpected genes missing from proseg output: {missing[:5]}"
n_dropped_probes = len(missing)
need = [g for g in need if g in gene_pos]
missing_tested = [g for g in tested if g not in gene_pos]
assert not missing_tested, f"tested genes missing from proseg output: {missing_tested}"
print(f"[{time.time()-t0:.0f}s] A1-j: {n_dropped_probes} ultra-rare probes dropped "
      f"from marker space; all 17 tested genes present", flush=True)
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
# M must carry ALL marker-space columns in centroid order (A1-h alignment);
# the 4 A1-j probes are absent from proseg -> all-zero columns (they summed to
# 18 transcripts section-wide), so the transfer space is unchanged in effect.
M = np.zeros((n_pc, len(marker_genes)), dtype=np.int32)
for j, g in enumerate(marker_genes):
    if g in need_pos:
        M[:, j] = Mall[:, need_pos[g]]
Gm = Mall[:, [need_pos[g] for g in tested]].astype(np.int64)
assert (Gm.sum(axis=0) > 0).all()

# ---------- label transfer (identical transform; A1-h target space) ----------
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

# ---------- region mapping (nearest vendor centroid, sealed procedure) ----------
pc_xy = pc[["centroid_x", "centroid_y"]].to_numpy(np.float64)
tree = cKDTree(vendor_xy)
dist, nn = tree.query(pc_xy, k=1)
regions = regions_vendor[nn]
print(f"[{time.time()-t0:.0f}s] regions: {np.bincount(regions).tolist()}; "
      f"nn-dist median {np.median(dist):.2f} um", flush=True)
assert np.median(dist) < 15.0

# ---------- T(A_proseg) for the 30 vendor-evaluable rows ----------
types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
n_types = len(types_needed)
SHIFT = 1
type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
n_slots = n_types + SHIFT
cell_slot = regions * n_slots + np.array(
    [type_shift.get(l, 0) for l in labels], dtype=np.int64)

rows = []
n_type_absent = 0
# The registered E2 universe is the 30 vendor-evaluable rows in scout_bounds.csv
# (the 8 Myeloid-dependent pairs are vendor-side absent by the frozen class floor and are not
# part of any registered endpoint on this leg).
for _, r in b.iterrows():
    Lg, Rg = r.pair.split("->")
    region_id = int(r.region)
    S, R = TYPE_MAP[Lg], RECV_MAP[Rg]
    sS = cell_slot == region_id * n_slots + type_shift[S]
    sR = cell_slot == region_id * n_slots + type_shift[R]
    nL, nR = int(sS.sum()), int(sR.sum())
    if nL == 0 or nR == 0:
        # proseg-side type absence on a vendor-evaluable row (counted in the all-rows endpoint)
        n_type_absent += 1
        rows.append(dict(pair=r.pair, region=region_id, nL=nL, nR=nR,
                         N_L=0, N_R=0, T=None))
        continue
    N_L = int(Gm[sS, gidx[Lg]].sum())
    N_R = int(Gm[sR, gidx[Rg]].sum())
    rows.append(dict(pair=r.pair, region=region_id, nL=nL, nR=nR,
                     N_L=N_L, N_R=N_R, T=t_log2(N_L, N_R, nL, nR)))
P = pd.DataFrame(rows)
assert len(P) == 30 and set(P.pair) == set(b.pair)
n_excl = int(P["T"].isna().sum())
n_eval = 30 - n_excl
print(f"[{time.time()-t0:.0f}s] evaluable rows: {n_eval}/30 (excluded {n_excl}; "
      f"proseg-side type-absent {n_type_absent})", flush=True)
# Registered cap: <=10% over the 30 vendor-evaluable rows
assert n_excl <= 3, f"{n_excl}/30 non-evaluable (>10%) - ABORT"

# ---------- endpoints ----------
PE = P.dropna(subset=["T"])
out = b[["pair","region","control","T0","T_lo","T_hi"]].merge(
    PE[["pair","region","T","nL","nR"]], on=["pair","region"], how="left",
    validate="1:1")
assert len(out) == 30
assert ((out.T0 + 1e-9 >= out.T_lo) & (out.T0 - 1e-9 <= out.T_hi)).all()
out["inside"] = (out["T"] >= out.T_lo - 1e-9) & (out["T"] <= out.T_hi + 1e-9)
n_in_all = int(out.loc[out["T"].notna(), "inside"].sum())
sign_agree_all = int((((out["T"] > 0) == (out.T0 > 0))[out["T"].notna()]).sum())

summary = dict(
    design="E2 Proseg-interval coverage for the second breast section; frozen before compute",
    n_proseg_cells=int(n_pc),
    cells_with_marker_signal=int(has_marker.sum()),
    noise_gene_slot=dropped_noise_gene,
    primary_certified_endpoint="UNDEFINED: 0 certified rows on s6",
    interval_coverage_all_rows=dict(
        n_inside=n_in_all, n_evaluable=n_eval, fraction=(n_in_all / n_eval)),
    sign_agreement_all_rows=f"{sign_agree_all}/{n_eval}",
    proseg_side_type_absent_rows=n_type_absent,
)
out.to_csv(f"{D}/phase7_e2_coverage.csv", index=False)
with open(f"{D}/phase7_e2_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2), flush=True)
print("PHASE7 E2 COMPLETE", flush=True)
