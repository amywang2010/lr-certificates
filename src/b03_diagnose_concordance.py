"""B03 diagnostic: root-cause the concordance gate result.

Hypotheses:
  H1 (benign): per-cell deltas are small (+-1 boundary semantics, dropped sliver polygons);
               per-gene per-cell count CORRELATION ~ 1.0.
  H2 (fatal):  barcode ORDER mismatch between h5 and cells.parquet -> correlations ~ 0.
  H3 (benign subset): vendor assigns ~2.6% more molecules than PIP (nucleus-expansion /
               sliver drops) -> per-gene totals slightly higher in h5, uniformly.

Decisive outputs:
  - h5 barcode order vs cells.parquet cell_id order (exact string comparison)
  - per-gene per-cell Pearson r (h5 count vs PIP count) for 40 genes
  - delta histogram (h5 - pip) pooled
  - per-gene total ratio distribution
"""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import h5py
import json, time, pickle

DATA = "B03_project/data"

def main():
    t0 = time.time()
    with h5py.File(f"{DATA}/Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5",
                   "r") as f:
        g = f["matrix"]
        shape = g["shape"][:]
        barcodes = [b.decode() if isinstance(b, bytes) else b for b in g["barcodes"][:]]
        features = [b.decode() if isinstance(b, bytes) else b
                    for b in g["features"]["name"][:]]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
    cells_meta = pq.read_table(
        f"{DATA}/Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.parquet").to_pandas()
    vid = cells_meta.cell_id.to_numpy().astype(str)
    bc = np.array(barcodes)
    print(f"barcodes == cells.parquet order: {np.array_equal(bc, vid)}", flush=True)
    print(f"barcodes as sets equal: {set(bc) == set(vid)}", flush=True)
    # build mapping h5 row -> cells.parquet row (by string id) for correct comparison
    pos_in_vid = pd.Series(np.arange(len(vid)), index=vid)
    remap = pos_in_vid.reindex(bc).to_numpy()  # h5 row i corresponds to vid row remap[i]
    print(f"rows needing reorder: {int((remap != np.arange(len(remap))).sum()):,}",
          flush=True)

    tx = pd.read_parquet(f"{DATA}/tx.parquet")
    owner = tx.cell_idx.to_numpy()
    feat = tx.feature_name.to_numpy()
    # invert cid_map: our polygon idx -> cell_id string
    with open(f"{DATA}/donor_map.pkl", "rb") as f:
        dm = pickle.load(f)
    cid_map = dm["cid_map"]  # cell_id -> poly idx
    inv = np.full(len(vid), -1, dtype=np.int64)
    # cell_id in boundaries parquet may differ in format from cells.parquet; map via str
    for cid_str, pidx in cid_map.items():
        s = str(cid_str)
        if s in pos_in_vid.index:
            inv[pos_in_vid[s]] = pidx
    print(f"cell_id->poly mapping coverage: {(inv>=0).mean():.4f} "
          f"({time.time()-t0:.0f}s)", flush=True)

    rng = np.random.default_rng(20260905)
    test_genes = rng.choice(len(features), size=min(40, len(features)), replace=False)
    rs, ratios, deltas = [], [], []
    for gi in test_genes:
        gname = features[gi]
        sel = indices == gi
        cell_rows = np.repeat(np.arange(len(barcodes)), np.diff(indptr))
        cell_rows_g = cell_rows[sel]
        cnts_g = data[sel]
        # h5 per-cell counts mapped to cells.parquet order
        h5v = np.zeros(len(vid))
        rr = remap[cell_rows_g]
        ok = rr >= 0
        np.add.at(h5v, rr[ok], cnts_g[ok])
        # PIP per-cell counts: our poly idx -> cells.parquet row
        m = feat == gname
        pipp = owner[m]
        pipp = pipp[pipp >= 0]
        pipv = np.zeros(len(vid))
        # our poly idx -> cell row: invert inv (sizes: n_polys x n_cells)
        poly2row = -np.ones(int(inv.max()) + 1 if inv.size else 1, dtype=np.int64)
        poly2row[inv[inv >= 0]] = np.flatnonzero(inv >= 0)
        valid = pipp < len(poly2row)
        rr2 = poly2row[pipp[valid]]
        rr2 = rr2[rr2 >= 0]
        np.add.at(pipv, rr2, 1)
        if h5v.sum() == 0 and pipv.sum() == 0:
            continue
        if h5v.std() > 0 and pipv.std() > 0:
            rs.append(float(np.corrcoef(h5v, pipv)[0, 1]))
        if h5v.sum() > 0:
            ratios.append(float(pipv.sum() / h5v.sum()))
        deltas.extend((h5v - pipv)[(h5v > 0) | (pipv > 0)][:2000].tolist())
    rs = np.array(rs); ratios = np.array(ratios); deltas = np.array(deltas)
    print(f"per-gene per-cell Pearson r: median={np.median(rs):.5f} "
          f"min={rs.min():.5f} max={rs.max():.5f} (n={len(rs)})", flush=True)
    print(f"per-gene total ratio pip/h5: median={np.median(ratios):.4f} "
          f"min={ratios.min():.4f} max={ratios.max():.4f}", flush=True)
    hist, edges = np.histogram(deltas, bins=[-10, -3, -2, -1, 0, 1, 2, 3, 10])
    print("delta (h5-pip) histogram:", dict(zip([f"{int(edges[i])}..{int(edges[i+1])-1}"
                                                 for i in range(len(hist))],
                                                hist.tolist())), flush=True)
    out = {"order_identical": bool(np.array_equal(bc, vid)),
           "order_sets_equal": bool(set(bc) == set(vid)),
           "pearson_median": float(np.median(rs)) if len(rs) else None,
           "pearson_min": float(rs.min()) if len(rs) else None,
           "ratio_median": float(np.median(ratios)) if len(ratios) else None,
           "delta_hist": hist.tolist()}
    with open(f"{DATA}/concordance_diagnosis.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"DONE ({time.time()-t0:.0f}s)", flush=True)

if __name__ == "__main__":
    main()
