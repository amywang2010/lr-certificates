"""B03 Phase 2 prep (step A): vendor-side, Proseg-INDEPENDENT primitives.

Pre-registered in B03_PHASE2_PREREG.md (frozen BEFORE any Proseg output exists).
Computes only artifacts that do not depend on Proseg:

  1. Vendor per-type marker centroids (label-transfer targets, prereg section 2).
     Marker space: panel genes EXCLUDING all 20 tested LR genes (same discipline
     as the scout labels). Transform, identical for vendor and Proseg cells:
        p = count / cell_total_panel_counts      (relative profile)
        x = log1p(p * scale),  scale = median vendor cell panel total
     Type centroid = mean of x over vendor cells of that type.

  2. Region-mapping inputs (prereg section 1): vendor cell centroids (microns)
     plus each cell's frozen region, in cells_meta ROW space, so the analysis
     step can build a cKDTree and assign each Proseg cell the region of its
     nearest vendor cell (geometry only, no expression).

Proseg-independence argument (prereg discipline): the ONLY inputs are the
vendor cell-by-gene h5, the panel TSV, cells_meta.parquet, and the scout's
frozen labels/regions/donor_map. The transcript file (and therefore anything
Proseg produces) is never read here.
"""

import numpy as np
import pandas as pd
import pickle
import h5py
import sys
import time

sys.path.insert(0, "B03_project/src")
from b03_scout import ALL_PAIRS  # frozen pair list (single source of truth)

DATA = "B03_project/data"
OUT = "B03_project/data/vendor_markers.npz"

t0 = time.time()

tested = sorted({g for p in ALL_PAIRS for g in p})
print(f"[{time.time()-t0:.0f}s] {len(tested)} tested genes (excluded from marker panel)",
      flush=True)

# ---------- frozen scout artifacts ----------
labels_poly = np.load(f"{DATA}/labels.npy", allow_pickle=True)      # polygon-index space
regions_poly = np.load(f"{DATA}/regions.npy")                        # polygon-index space
with open(f"{DATA}/donor_map.pkl", "rb") as f:
    dm = pickle.load(f)
cells_meta = pd.read_parquet(f"{DATA}/cells_meta.parquet")

# polygon idx -> cells_meta row (identical alignment code to the scout)
cid_map = {str(k): v for k, v in dm["cid_map"].items()}
cm_ids = cells_meta.cell_id.astype(str).to_numpy()
poly_row = np.full(max(cid_map.values()) + 1, -1, dtype=np.int64)
for row, cid in enumerate(cm_ids):
    p = cid_map.get(cid, -1)
    if p >= 0:
        poly_row[p] = row
assert (poly_row >= 0).all(), "unmapped polygons present"
n_cells = len(cells_meta)
print(f"[{time.time()-t0:.0f}s] alignment: {n_cells} cells", flush=True)

labels_row = labels_poly[poly_row]      # cells_meta ROW space
regions_row = regions_poly[poly_row]    # cells_meta ROW space
assert set(np.unique(regions_row)) == {0, 1, 2}

# ---------- vendor cell-by-gene: detect layout, never assume ----------
h5_path = f"{DATA}/Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5"
with h5py.File(h5_path, "r") as h5:
    feat_names_h5 = h5["matrix/features/name"][:]
    features = [b.decode() if isinstance(b, bytes) else str(b) for b in feat_names_h5]
    shape = h5["matrix/shape"][:]
    barcodes = [b.decode() if isinstance(b, bytes) else str(b)
                for b in h5["matrix/barcodes"][:]]
    keys = set(h5["matrix"].keys())
    if "indptr" in keys:
        indptr = h5["matrix/indptr"][:]
        indices = h5["matrix/indices"][:]
        data = h5["matrix/data"][:]
    else:
        # feature-blocked layout (e.g. barcodes_0); concatenate all blocks
        blocks = sorted(k for k in keys if k.startswith("barcodes_"))
        assert blocks, f"no CSC arrays in matrix/: {sorted(keys)}"
        ip, ii, dd = [], [], []
        off = 0
        for b in blocks:
            g = h5[f"matrix/{b}"]
            bip = g["indptr"][:]
            bii = g["indices"][:]
            bdd = g["data"][:]
            if ip:
                bip = bip[1:] + off   # drop duplicated boundary, offset indices
            else:
                bip = bip + off
            ip.append(bip.astype(bip.dtype))
            ii.append(bii + off)
            dd.append(bdd)
            off = off + len(features)   # blocked over the OTHER axis; corrected below
        # NOTE: blocking axis is resolved by the indptr-length check below.
        indptr = np.concatenate(ip)
        indices = np.concatenate(ii)
        data = np.concatenate(dd)

n_feat, n_bc = int(shape[0]), int(shape[1])
assert n_bc == n_cells, f"h5 shape {shape} vs cells_meta {n_cells}"
match = np.mean([b == c for b, c in zip(barcodes[:1000], cm_ids[:1000])])
assert match > 0.99, f"barcode order mismatch ({match})"
print(f"[{time.time()-t0:.0f}s] h5: {n_feat} features x {n_bc} cells; "
      f"indptr len {len(indptr)}", flush=True)

if len(indptr) == n_bc + 1:
    # CSC with CELLS as columns (indptr per cell): indices[] = feature ids.
    # (Shape is [features x cells]; the loader's concordance read the same layout.)
    pass
else:
    raise AssertionError(f"indptr len {len(indptr)} != n_cells+1 ({n_bc+1}); "
                         f"layout not recognized")

# ---------- marker panel (prereg section 2) ----------
panel = pd.read_csv(f"{DATA}/Xenium_FFPE_Human_Breast_Cancer_Rep1_panel.tsv", sep="\t")
panel_names = panel.Name.dropna().astype(str).tolist()
feat_set = set(features)
marker_genes = [g for g in panel_names if (g not in tested) and (g in feat_set)]
n_missing = int(sum(1 for g in panel_names if g not in feat_set))
print(f"[{time.time()-t0:.0f}s] marker panel: {len(marker_genes)} genes "
      f"({n_missing} panel genes absent from h5; tested genes excluded)", flush=True)
assert len(marker_genes) >= 100, "marker panel unexpectedly small"

# ---------- accumulate marker submatrix (vectorized) ----------
n_markers = len(marker_genes)
M = np.zeros((n_cells, n_markers), dtype=np.int32)
feat_to_col = np.full(n_feat, -1, dtype=np.int64)
for j, g in enumerate(marker_genes):
    feat_to_col[features.index(g)] = j
data = data.astype(np.int64, copy=False)
col_ids = np.repeat(np.arange(n_cells, dtype=np.int64), np.diff(indptr))
j = feat_to_col[indices]
sel = j >= 0
assert len(col_ids) == len(indices) == len(data)
M[col_ids[sel], j[sel]] = data[sel]  # (cell, gene) unique within a CSC column
assert (M.sum(axis=0) > 0).all(), "some marker genes got zero total counts (read bug)"
print(f"[{time.time()-t0:.0f}s] marker submatrix: {M.shape}, total {int(M.sum()):,}",
      flush=True)

# ---------- transform (prereg section 2, verbatim) ----------
C = M.sum(axis=1).astype(np.int64)
assert (C > 0).mean() > 0.99, "too many cells with zero panel counts"
scale = float(np.median(C))
X = np.log1p((M / np.maximum(C, 1)[:, None]) * scale).astype(np.float32)
print(f"[{time.time()-t0:.0f}s] transform done: scale={scale:.1f}", flush=True)

# ---------- per-type centroids ----------
# Pre-specified (amendment record, frozen before any Proseg output): transfer TARGETS are
# types with >=50 vendor cells. Rationale: (a) centroids from <50 (down to 1) cells
# are statistically meaningless argmax targets; (b) no type with <50 cells appears
# in TYPE_MAP/RECV_MAP for any of the 25 pairs, so the T statistics are unaffected.
types_all = sorted(set(labels_row.tolist()) - {"Unknown"})
counts_all = {t: int((labels_row == t).sum()) for t in types_all}
types_needed = [t for t in types_all if counts_all[t] >= 50]
dropped = {t: counts_all[t] for t in types_all if counts_all[t] < 50}
print(f"[{time.time()-t0:.0f}s] dropped transfer targets (<50 cells): {dropped}",
      flush=True)
assert not ({"Breast cancer", "Fibroblasts", "Macrophages", "Endothelial cells",
             "Smooth muscle cells", "Breast glandular cells", "T cells"}
            - set(types_needed)), "an analysis-relevant type was dropped"
centroids = np.stack([X[labels_row == t].mean(axis=0) for t in types_needed])
counts_by_type = {t: counts_all[t] for t in types_needed}
print(f"[{time.time()-t0:.0f}s] centroids: {counts_by_type}", flush=True)

# ---------- save ----------
np.savez_compressed(
    OUT,
    marker_genes=np.array(marker_genes, dtype=object),
    type_names=np.array(types_needed, dtype=object),
    centroids=centroids.astype(np.float32),
    scale=np.float64(scale),
    vendor_centroids_xy=cells_meta[["x_centroid", "y_centroid"]].to_numpy(np.float64),
    regions_row=regions_row.astype(np.int8),
)
print(f"[{time.time()-t0:.0f}s] saved -> {OUT}", flush=True)
print("VENDOR PREP COMPLETE", flush=True)
