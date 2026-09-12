"""B03 Phase 2 (step B): Proseg replication analysis per B03_PHASE2_PREREG.md.

Endpoint (prereg, primary): fraction of the 75 (pair, region) rows with
T(A_proseg, K0) inside the scout's certified interval [T_lo, T_hi]. Gate >= 90%.
Secondary: sign agreement on rows certified by the matched null (q <= 0.10).

Input schema (verified from proseg 3.2.0 source, output.rs/main.rs — not assumed):
  - counts:   gzipped MatrixMarket coordinate (cells x genes), 1-indexed, written by
              write_sparse_mtx from params.foreground_counts ("point estimate of
              transcript counts per cell"). Values are x.xx floats; counts are
              converted with np.rint (definitional: they ARE counts).
  - genes:    parquet with column `gene`, row i = gene index i (0-based) of the mtx
              columns (written in gene-index order by write_gene_metadata).
  - cells:    parquet with columns cell (0..ncells-1), original_cell_id,
              centroid_x, centroid_y, ...; row order = cell id order (writer
              enumerates 0..ncells-1), so mtx row i (1-based) = cells row i (0-based).

Machinery reuse (prereg section 4): EPS and t_log2 IMPORTED from b03_scout.py.
Labels: frozen vendor centroids (vendor_markers.npz, identical transform, step A).
Regions: nearest vendor centroid (prereg section 1). Background cells (id <= 0),
if any exist, are excluded per prereg Amendment A1 and the count is reported.
"""

import gzip
import json
import time
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial import cKDTree

sys.path.insert(0, "B03_project/src")
from b03_scout import ALL_PAIRS, TYPE_MAP, RECV_MAP, EPS, t_log2

DATA = "B03_project/data"
RES = "B03_project/results"
t0 = time.time()

tested = sorted({g for p in ALL_PAIRS for g in p})
gidx = {g: i for i, g in enumerate(tested)}

# ---------- frozen vendor-side artifacts ----------
vm = np.load(f"{DATA}/vendor_markers.npz", allow_pickle=True)
centroids = vm["centroids"].astype(np.float64)
type_names = [str(t) for t in vm["type_names"]]
marker_genes = [str(g) for g in vm["marker_genes"]]
scale = float(vm["scale"])
vendor_xy = vm["vendor_centroids_xy"]
regions_vendor = vm["regions_row"].astype(int)
assert centroids.shape[0] == len(type_names) and centroids.shape[1] == len(marker_genes)

# ---------- scout certified intervals ----------
b = pd.read_csv(f"{RES}/scout_bounds.csv")
if b["control"].dtype != bool:
    b["control"] = b["control"].astype(str).map({"True": True, "False": False})
assert len(b) == 75
fin = pd.read_csv(f"{RES}/scout_final_matched.csv")

# ---------- proseg outputs (schema per source, asserted at load) ----------
# amendment record: proseg 3.2.0 writes metadata files named EXACTLY as the --output-*
# arguments (cells/genes, no .parquet suffix) — verified on disk (PAR1 magic).
with gzip.open(f"{DATA}/proseg_out/counts", "rb") as f:
    Coo = mmread(f).tocoo()
cells_meta_p = pd.read_parquet(f"{DATA}/proseg_out/cells")
genes_p = pd.read_parquet(f"{DATA}/proseg_out/genes")
print(f"[{time.time()-t0:.0f}s] mtx {Coo.shape} nnz={Coo.nnz:,}; cells "
      f"{cells_meta_p.shape}; genes {genes_p.shape}", flush=True)
assert Coo.shape[0] == len(cells_meta_p), "mtx rows != cells rows"
# A2b: gene columns may be n_named (all named) or n_named-1 (internal noise slot,
# handled below); anything else is an unrecognized schema change -> abort.
assert Coo.shape[1] in (len(genes_p), len(genes_p) - 1), \
    f"mtx cols {Coo.shape[1]} vs genes metadata {len(genes_p)} - schema changed"
gene_names = genes_p["gene"].astype(str).to_numpy()
assert len(set(gene_names)) == len(gene_names), "duplicate gene names"

# ---------- prereg Amendment A2: cell universe + gene-alignment verification ----------
# (a) counts mtx is documented (sampler.rs:249) as foreground (non-noise) counts, so
# cells.parquet should contain only real cells; drop rows with the vendor unassigned
# marker original_cell_id == "0" iff present (expected 0) and report the count.
pc = cells_meta_p
is_bg = pc["original_cell_id"].astype(str) == "0"
n_bg_cells = int(is_bg.sum())
if n_bg_cells:
    keep = (~is_bg).to_numpy().nonzero()[0]
    Coo = Coo.tocsr()[keep].tocoo()
    pc = pc[~is_bg].reset_index(drop=True)
n_pc = len(pc)
print(f"[{time.time()-t0:.0f}s] background-cell rows dropped (A2a): {n_bg_cells} "
      f"(expected 0); universe {n_pc} cells", flush=True)
ids_kept = pc["cell"].to_numpy(np.int64)
assert (np.diff(ids_kept) > 0).all(), "kept cell ids not strictly increasing"
assert (n_bg_cells > 0) or (ids_kept == np.arange(n_pc)).all(), \
    "no background rows dropped yet ids != 0..n-1 - unexpected schema"
print(f"[{time.time()-t0:.0f}s] cell-id/mtx-row mapping verified "
      f"(ids {ids_kept[0]}..{ids_kept[-1]})", flush=True)

# (b) mtx gene-column count: n_named, or n_named-1 if an internal noise slot exists
# (sampler.rs:445). The n-1 case drops the LAST named gene (noise slot is the last
# internal index); ABORT if that gene is needed for the analysis.
mtx_n_genes = Coo.shape[1]
dropped_noise_gene = None
if mtx_n_genes == len(gene_names):
    gene_names_eff = gene_names
elif mtx_n_genes == len(gene_names) - 1:
    dropped_noise_gene = str(gene_names[-1])
    assert dropped_noise_gene not in set(marker_genes) | set(tested), \
        f"noise-slot drop would remove analysis gene {dropped_noise_gene} - ABORT"
    gene_names_eff = gene_names[:-1]
    print(f"[{time.time()-t0:.0f}s] mtx has ngenes-1 columns; noise slot = last named "
          f"gene {dropped_noise_gene!r} (not analysis-needed) - aligned", flush=True)
else:
    raise AssertionError(f"mtx gene cols {mtx_n_genes} matches neither n_named "
                         f"({len(gene_names)}) nor n_named-1 - schema changed, ABORT")

# (c) per-gene alignment proof: genes.parquet total_count is computed from ALL input
# transcripts (write_gene_metadata), while the mtx holds foreground (non-noise)
# counts. Therefore mtx col total <= gene total is a HARD invariant; exceeding it
# (beyond float/point-estimate slack) proves column misalignment. Check all needed
# genes and abort on any violation.
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
assert not viol, f"GENE MISALIGNMENT DETECTED (mtx total > all-tx total): {viol[:5]}"
ratios = col_sums / np.maximum([float(gene_totals[g]) for g in need], 1.0)
print(f"[{time.time()-t0:.0f}s] alignment proof: {len(need)} needed genes, 0 invariant "
      f"violations; foreground/all-tx ratio range [{ratios.min():.3f}, "
      f"{ratios.max():.3f}]", flush=True)

n_need = len(need)
Sub = Sub.tocoo()

# integer count matrix (rounding point estimates to counts)
Mall = np.zeros((n_pc, n_need), dtype=np.int32)
np.add.at(Mall, (Sub.row.astype(np.int64), Sub.col.astype(np.int64)),
          np.rint(Sub.data).astype(np.int64))
need_pos = {g: i for i, g in enumerate(need)}
M = Mall[:, [need_pos[g] for g in marker_genes]].astype(np.int32)   # markers, frozen order
Gm = Mall[:, [need_pos[g] for g in tested]].astype(np.int64)        # tested genes
assert (Gm.sum(axis=0) > 0).all(), "some tested gene has zero total counts"
print(f"[{time.time()-t0:.0f}s] submatrix built: {Mall.shape}, total "
      f"{int(Mall.sum()):,}", flush=True)

# ---------- label transfer (prereg section 2) ----------
Ctot = M.sum(axis=1).astype(np.int64)
has_marker = Ctot > 0
Xp = np.log1p((M / np.maximum(Ctot, 1)[:, None]) * scale).astype(np.float64)
Xn = Xp / np.maximum(np.linalg.norm(Xp, axis=1, keepdims=True), 1e-12)
Cn = centroids / np.maximum(np.linalg.norm(centroids, axis=1, keepdims=True), 1e-12)
labels = np.full(n_pc, "Unknown", dtype=object)
sim = Xn[has_marker] @ Cn.T
labels[has_marker] = np.array(type_names, dtype=object)[sim.argmax(axis=1)]
print(f"[{time.time()-t0:.0f}s] labels: "
      f"{pd.Series(labels).value_counts().head(8).to_dict()}", flush=True)

# ---------- region mapping (prereg section 1) ----------
pc_xy = pc[["centroid_x", "centroid_y"]].to_numpy(np.float64)
tree = cKDTree(vendor_xy)
dist, nn = tree.query(pc_xy, k=1)
regions = regions_vendor[nn]
print(f"[{time.time()-t0:.0f}s] regions: {np.bincount(regions).tolist()}; "
      f"nn-dist median {np.median(dist):.2f} um (p95 {np.quantile(dist, 0.95):.2f})",
      flush=True)
assert np.median(dist) < 15.0, "proseg cells unreasonably far from vendor cells"

# ---------- T(A_proseg) for the 75 rows (identical machinery) ----------
types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
n_types = len(types_needed)
SHIFT = 1
type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
n_slots = n_types + SHIFT
UNKNOWN_SLOT = 0

cell_slot = regions * n_slots + np.array(
    [type_shift.get(l, UNKNOWN_SLOT) for l in labels], dtype=np.int64)

rows = []
for (Lg, Rg) in ALL_PAIRS:
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
assert len(P) == 75
# prereg Amendment A3: non-evaluable rows (empty sender or receiver population under
# A_proseg) are excluded from the denominator and reported; >10% non-evaluable aborts.
n_excl = int(P["T"].isna().sum())
n_eval = 75 - n_excl
print(f"[{time.time()-t0:.0f}s] evaluable rows: {n_eval}/75 (excluded {n_excl}, "
      f"empty sender/receiver population)", flush=True)
assert n_excl <= 7, f"{n_excl}/75 rows non-evaluable (>10%) - mapping failure, ABORT"
assert P["T"].notna().all() or n_eval > 0
print(f"[{time.time()-t0:.0f}s] T(A_proseg) computed for 75 rows", flush=True)

# ---------- endpoints ----------
PE = P.dropna(subset=["T"])
out = b[["pair", "region", "control", "T0", "T_lo", "T_hi"]].merge(
    PE[["pair", "region", "T", "nL", "nR"]], on=["pair", "region"],
    how="left", validate="1:1")
assert len(out) == 75
assert ((out.T0 + 1e-9 >= out.T_lo) & (out.T0 - 1e-9 <= out.T_hi)).all(), \
    "vendor T0 outside its own bounds - bounds file mismatch"
out["inside"] = (out["T"] >= out.T_lo - 1e-9) & (out["T"] <= out.T_hi + 1e-9)
inside_eval = out.loc[out["T"].notna(), "inside"]
n_in = int(inside_eval.sum())
frac = n_in / n_eval
fin2 = fin[["pair", "region", "certified_pos_matched", "certified_neg_matched"]]
out = out.merge(fin2, on=["pair", "region"], validate="1:1")
cert_mask = out.certified_pos_matched | out.certified_neg_matched
cert_eval = cert_mask & out["T"].notna()
sign_ok = int((((out["T"] > 0) == out.certified_pos_matched)[cert_eval]).sum())

summary = dict(
    prereg="B03_PHASE2_PREREG.md v1.0 (+ Amendments A1, A2)",
    proseg_version="3.2.0",
    command=open("B03_project/logs/proseg_cmd.txt").read().strip(),
    counts_format="gzip MatrixMarket (write_sparse_mtx), rint->int64",
    n_background_cells_excluded=n_bg_cells,
    noise_gene_slot=dropped_noise_gene,
    n_proseg_cells=int(n_pc),
    cells_with_marker_signal=int(has_marker.sum()),
    primary=dict(n_inside=n_in, n_evaluable=n_eval, n_excluded=n_excl,
                 fraction=frac, gate_pass=bool(frac >= 0.90)),
    certified_rows=int(cert_mask.sum()),
    sign_agreement_on_certified=f"{sign_ok}/{int(cert_mask.sum())}",
)
out.to_csv(f"{RES}/phase2_proseg_coverage.csv", index=False)
with open(f"{RES}/phase2_proseg_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2), flush=True)
print("PHASE2 ANALYSIS COMPLETE", flush=True)
