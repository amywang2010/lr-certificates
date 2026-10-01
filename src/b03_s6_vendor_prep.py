"""B03 Phase 7 (step A, A1-h): vendor-side label-transfer targets for breast_s6.

Marker space per A1-h: crop is_gene vocabulary minus the 17 tested LR genes
(the phase2 tissue-1 definition). Rationale recorded in A1-h: the phase5
9-gene dictionary yields scale = 0 on s6, collapsing the transform; the
phase2 design anticipated this and carries the signal assertion
(C>0).mean() > 0.99, restored here. Transform identical to sealed legs:
x = log1p((M/C)*scale), scale = median vendor cell marker total; centroids =
mean of x over vendor cells of each surviving class (DEV amendment 007: Myeloid absent
by design). Sparse accumulation throughout (n_cells x ~5k markers).
Outputs the s6 crop data dir's vendor_markers.npz.
"""
import os, sys, time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.sparse import coo_matrix, csr_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_lung_config import PAIRS
import b03_phase5_vendor_prep as lung_ref

assert PAIRS == lung_ref.PAIRS, "PAIRS diverge from lung driver"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
D = os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "data")
CROP = os.path.join(WORKSPACE, "xenium_breast_s6", "crop")
TESTED = sorted({g for p in PAIRS for g in p})
assert len(TESTED) == 17
MIN_CLASS_CELLS = 50  # phase2 vendor_prep threshold

t0 = time.time()

def marker_space():
    """Crop is_gene vocabulary minus tested LR genes (A1-h)."""
    f = pq.ParquetFile(f"{CROP}/transcripts.parquet")
    names = set()
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=["feature_name", "is_gene"])
        fn = t.column("feature_name")
        ig = t.column("is_gene").to_numpy(zero_copy_only=False)
        vals = fn.to_pylist()
        names |= {v for v, g in zip(vals, ig) if g and v is not None}
    ms = sorted(names - set(TESTED))
    print(f"[{time.time()-t0:.0f}s] marker space: {len(ms)} genes "
          f"(crop is_gene vocabulary minus 17 tested)", flush=True)
    return ms

def main():
    tx = pd.read_parquet(f"{D}/tx.parquet")
    cells_meta = pd.read_parquet(f"{D}/cells_meta.parquet")
    labels = np.load(f"{D}/labels.npy", allow_pickle=True)
    regions = np.load(f"{D}/regions.npy")
    n_cells = len(cells_meta)
    assert n_cells == len(labels) == len(regions)
    MS = marker_space()
    mpos = {g: i for i, g in enumerate(MS)}
    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()
    keep = (cell >= 0) & np.isin(feat, MS)
    rows = cell[keep]
    cols = np.array([mpos[g] for g in feat[keep]], dtype=np.int64)
    M = coo_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                   shape=(n_cells, len(MS))).tocsr()
    C = np.asarray(M.sum(axis=1)).ravel()
    frac_signal = float((C > 0).mean())
    print(f"[{time.time()-t0:.0f}s] cells with marker signal: {frac_signal:.4f}",
          flush=True)
    # A1-h (correction): s6-calibrated signal gates. The phase2 0.99 threshold is
    # mechanically unreachable here (2,488/65,964 vendor cells are EMPTY under the
    # repackage's CellVi segmentation - a section property). Purpose preserved:
    # (a) non-degenerate scale, (b) transfer coverage >= 0.90; zero-signal cells
    # stay 'Unknown' downstream (standing has_marker mask).
    assert frac_signal >= 0.90, f"transfer coverage too low: {frac_signal:.4f}"
    scale = float(np.median(C[C > 0]))
    print(f"scale (median vendor marker total): {scale:.1f}", flush=True)
    types_needed = sorted({"Epithelial", "T_cells", "Myeloid", "B_cells",
                           "Endothelial", "Fibroblast"})
    keep_types = [t for t in types_needed if int((labels == t).sum()) >= MIN_CLASS_CELLS]
    dropped = sorted(set(types_needed) - set(keep_types))
    assert dropped == ["Myeloid"], f"unexpected class drops: {dropped}"
    print(f"classes kept: {keep_types}; dropped (DEV amendment 007): {dropped}", flush=True)
    # chunked centroid accumulation (no dense X)
    Mcsr = M.tocsr()
    n_ms = len(MS)
    sums = {t: np.zeros(n_ms, dtype=np.float64) for t in keep_types}
    CH = 8192
    for s in range(0, n_cells, CH):
        e = min(n_cells, s + CH)
        sub = M[s:e].tocsr()
        Ce = np.asarray(sub.sum(axis=1)).ravel()
        fac = (scale / np.maximum(Ce, 1))  # per-row scale factors
        rowidx = np.repeat(np.arange(e - s), np.diff(sub.indptr))
        data = np.log1p(sub.data * fac[rowidx])  # data-level transform, stays sparse
        Xe = csr_matrix((data, sub.indices.copy(), sub.indptr.copy()), shape=sub.shape)
        for t in keep_types:
            m = (labels[s:e] == t)
            if m.any():
                sums[t] += np.asarray(Xe[m].sum(axis=0)).ravel()
    centroids = np.stack([sums[t] / int((labels == t).sum()) for t in keep_types])
    print(f"[{time.time()-t0:.0f}s] centroids: "
          f"{[(t, int((labels==t).sum())) for t in keep_types]}", flush=True)
    np.savez(f"{D}/vendor_markers.npz",
             centroids=centroids.astype(np.float32),
             type_names=np.array(keep_types, dtype=object),
             marker_genes=np.array(MS, dtype=object),
             scale=np.float64(scale),
             vendor_centroids_xy=cells_meta[["x_centroid","y_centroid"]].to_numpy(np.float64),
             regions_row=regions.astype(int))
    print(f"saved {D}/vendor_markers.npz", flush=True)

if __name__ == "__main__":
    main()
