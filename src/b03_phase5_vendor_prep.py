"""B03 Phase 5 (step A): vendor-side label-transfer targets for the lung crop.

Mirrors b03_phase2_vendor_prep.py exactly (transform x = log1p((M/Ctot)*scale),
scale = median vendor cell marker total; per-class centroids from vendor cells),
with the addendum-A3 marker dictionary (tested-LR genes excluded mechanically).
Outputs xenium_lung/crop/data/vendor_markers.npz: centroids, type_names,
marker_genes, scale, vendor_centroids_xy, regions_row.
"""
import os, sys, pickle
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

D = "xenium_lung/crop/data"
PAIRS = [("CD274","PDCD1"),("CXCL12","CXCR4"),("PECAM1","KDR"),
 ("ERBB2","EGFR"),("ERBB2","PDCD1"),("ESR1","PGR"),("PGR","ESR1"),
 ("CD274","CTLA4"),("CD274","CD8A"),("CXCL12","CCR7"),("PECAM1","PDCD1"),
 ("PECAM1","CTLA4"),("CD163","CCR7"),("CD163","CXCR4"),("CDH1","EGFR"),
 ("CDH1","ERBB2"),("MS4A1","CD274"),("CD68","CD274")]
TESTED = sorted({g for p in PAIRS for g in p})
MARKERS_RAW = {
 "Epithelial": ["EPCAM", "MUC1"],
 "T_cells": ["CD3E", "CD4", "CD8A"],
 "Myeloid": ["CD68", "FCGR3A"],
 "B_cells": ["MS4A1", "CD79A"],
 "Endothelial": ["PECAM1", "CLDN5", "CDH5", "EGFL7"],
 "Fibroblast": ["PDGFRA", "FAP"],
}
MARKERS = sorted({g for gs in MARKERS_RAW.values() for g in gs if g not in TESTED})
MIN_CLASS_CELLS = 50  # phase2 vendor_prep threshold

def main():
    tx = pd.read_parquet(f"{D}/tx.parquet")
    cells_meta = pd.read_parquet(f"{D}/cells_meta.parquet")
    labels = np.load(f"{D}/labels.npy", allow_pickle=True)
    regions = np.load(f"{D}/regions.npy")
    n_cells = len(cells_meta)
    assert n_cells == len(labels) == len(regions)
    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()
    # per-cell marker counts (vendors' assigned molecules)
    M = np.zeros((n_cells, len(MARKERS)), dtype=np.int64)
    mpos = {g: i for i, g in enumerate(MARKERS)}
    sel = (cell >= 0) & np.isin(feat, MARKERS)
    np.add.at(M, (cell[sel], np.array([mpos[g] for g in feat[sel]])), 1)
    C = M.sum(axis=1)
    scale = float(np.median(C))
    X = np.log1p((M / np.maximum(C, 1)[:, None]) * scale).astype(np.float32)
    print(f"scale (median vendor marker total): {scale:.1f}", flush=True)
    types_needed = sorted(MARKERS_RAW.keys())
    keep_types = [t for t in types_needed if int((labels == t).sum()) >= MIN_CLASS_CELLS]
    assert len(keep_types) == len(types_needed), \
        f"class under {MIN_CLASS_CELLS} cells: {set(types_needed)-set(keep_types)}"
    centroids = np.stack([X[labels == t].mean(axis=0) for t in keep_types])
    print(f"centroids: {[(t, int((labels==t).sum())) for t in keep_types]}", flush=True)
    np.savez(f"{D}/vendor_markers.npz",
             centroids=centroids.astype(np.float32),
             type_names=np.array(keep_types, dtype=object),
             marker_genes=np.array(MARKERS, dtype=object),
             scale=np.float64(scale),
             vendor_centroids_xy=cells_meta[["x_centroid","y_centroid"]].to_numpy(np.float64),
             regions_row=regions.astype(int))
    print(f"saved {D}/vendor_markers.npz", flush=True)

if __name__ == "__main__":
    main()
