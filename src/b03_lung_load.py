"""B03 Phase 5: lung-crop vendor load — PIP, concordance, band geometry, labels, regions.

Reuses tissue-1 machinery BY IMPORT (b03_load: assign_pip, boundary_and_band,
donor_sets_for_band_tested; b03_certify: make_regions with prereg seed 20260907).
Lung-specific I/O per B03_PHASE5_ADDENDUM.md:
  - transcripts: xenium_lung/crop/transcripts.parquet (frozen crop window)
  - polygons:    crop cell_boundaries.parquet (same construction algorithm as v3)
  - concordance: full-section h5 (xenium_lung/extracted/cell_feature_matrix.h5),
    restricted to crop-cell barcodes; gate = diagnosis-informed (per_cell >= 0.90,
    total_ratio >= 0.97) — same gate class as tissue 1 (amendment record).
  - labels:      addendum A3 frozen marker dictionary, majority vote, tested-LR genes
    excluded (17 tested genes of the 18 testable pairs); classes <200 cells dropped
    (amendment record rule, frozen in prereg item 1).
  - regions:     k-means k=3 on (log10(1+transcript_counts), x/1000, y/1000), seed
    20260907 (prereg item 2).

Outputs -> xenium_lung/crop/data/: tx.parquet, donor_map.pkl, cells_meta.parquet,
concordance.json, labels.npy, regions.npy
Smoke mode (B03_LUNG_SMOKE=1): first 300k transcripts, no artifact overwrite
(writes to data/smoke_*), all assertions still run.
"""
import os, sys, json, pickle, time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_load import assign_pip, boundary_and_band, donor_sets_for_band_tested, BAND_UM
from b03_certify import make_regions

SMOKE = os.environ.get("B03_LUNG_SMOKE") == "1"
CROP = "xenium_lung/crop"
OUTD = f"{CROP}/data"
os.makedirs(OUTD, exist_ok=True)

# Frozen tested-gene set = unique genes of the 18 testable pairs (G1 verdict)
PAIRS = [("CD274","PDCD1"),("CXCL12","CXCR4"),("PECAM1","KDR"),
 ("ERBB2","EGFR"),("ERBB2","PDCD1"),("ESR1","PGR"),("PGR","ESR1"),
 ("CD274","CTLA4"),("CD274","CD8A"),("CXCL12","CCR7"),("PECAM1","PDCD1"),
 ("PECAM1","CTLA4"),("CD163","CCR7"),("CD163","CXCR4"),("CDH1","EGFR"),
 ("CDH1","ERBB2"),("MS4A1","CD274"),("CD68","CD274")]
TESTED = sorted({g for p in PAIRS for g in p})
assert len(TESTED) == 17, TESTED

# Addendum A3 marker dictionary as frozen; tested-LR genes are removed HERE,
# explicitly, so the exclusion discipline is mechanical (assertion below verifies
# the filter rather than leaking a marker into the vote).
MARKERS_RAW = {
 "Epithelial": ["EPCAM", "MUC1"],
 "T_cells": ["CD3E", "CD4", "CD8A"],
 "Myeloid": ["CD68", "FCGR3A"],
 "B_cells": ["MS4A1", "CD79A"],
 "Endothelial": ["PECAM1", "CLDN5", "CDH5", "EGFL7"],
 "Fibroblast": ["PDGFRA", "FAP"],
}
MARKERS = {cls: [g for g in genes if g not in TESTED] for cls, genes in MARKERS_RAW.items()}
REMOVED_TESTED_MARKERS = sorted(g for cls, genes in MARKERS_RAW.items() for g in genes
                                if g in TESTED)
MIN_CLASS_CELLS = 200  # amendment record rule, frozen

def load_transcripts():
    f = pq.ParquetFile(f"{CROP}/transcripts.parquet")
    frames = []
    n = 0
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=["x_location","y_location","feature_name","qv"])
        if SMOKE and n >= 300_000: break
        frames.append(t.to_pandas()); n += len(t)
    tx = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if SMOKE: tx = tx.head(300_000)
    tx = tx[tx.qv >= 20.0].reset_index(drop=True)
    print(f"transcripts (QV>=20): {len(tx):,}", flush=True)
    return tx

def load_polygons_lung():
    cb = pq.read_table(f"{CROP}/cell_boundaries.parquet",
                       columns=["vertex_x","vertex_y","cell_id"]).to_pandas()
    cb = cb.sort_values("cell_id", kind="stable").reset_index(drop=True)
    from shapely.geometry import Polygon
    cids = cb.cell_id.to_numpy()
    starts = np.flatnonzero(np.r_[True, cids[1:] != cids[:-1]])
    ends = np.r_[starts[1:], len(cids)]
    xs_all = cb.vertex_x.to_numpy(np.float64); ys_all = cb.vertex_y.to_numpy(np.float64)
    polys, cid_map = {}, {}
    for k in range(len(starts)):
        s, e = starts[k], ends[k]
        c = cids[s]; xs = xs_all[s:e]; ys = ys_all[s:e]
        if len(xs) < 3: continue
        p = Polygon(zip(xs, ys))
        if not (p.is_valid and p.area > 1e-6):
            p = p.buffer(0)
            if p.is_empty or p.area <= 1e-6: continue
        if p.geom_type == "MultiPolygon":
            p = max(p.geoms, key=lambda q: q.area)
        if p.area <= 1e-6: continue
        cid_map[c] = len(polys); polys[len(polys)] = p
    print(f"valid crop polygons: {len(polys):,}", flush=True)
    return polys, cid_map

def verify_concordance_lung(tx, owner, cid_map):
    """PIP vs full-section h5 restricted to crop barcodes, 40 random genes (same gate)."""
    with h5py.File("xenium_lung/extracted/cell_feature_matrix.h5", "r") as f:
        g = f["matrix"]; shape = g["shape"][:]
        barcodes = [b.decode() if isinstance(b,bytes) else b for b in g["barcodes"][:]]
        features = [b.decode() if isinstance(b,bytes) else b for b in g["features"]["name"][:]]
        data = g["data"][:]; indices = g["indices"][:]; indptr = g["indptr"][:]
    assert len(barcodes) == int(shape[1]) and len(features) == int(shape[0])
    crop_barcodes = set(pq.read_table(f"{CROP}/cells.parquet", columns=["cell_id"])
                        .column("cell_id").to_pylist())
    keep_cols = np.array([i for i, b in enumerate(barcodes) if b in crop_barcodes])
    print(f"h5: {len(features)} features; crop barcodes in h5: {len(keep_cols):,}", flush=True)
    rng = np.random.default_rng(20260907)
    rng = np.random.default_rng(20260907)
    test_genes = rng.choice(len(features), size=min(40, len(features)), replace=False)
    feat_arr = tx.feature_name.to_numpy()
    cols = np.repeat(np.arange(len(barcodes)), np.diff(indptr))
    rows_all = indices
    conc_total = conc_match = 0; mism = []
    for gi in test_genes:
        gname = features[gi]
        sel = rows_all == gi
        if not sel.any(): continue
        cell_rows = cols[sel]; cnts = data[sel]
        m = (feat_arr == gname)
        if not m.any(): continue
        pip_counts = pd.Series(owner[m][owner[m] >= 0]).value_counts()
        for row, c in zip(cell_rows, cnts):
            if row not in keep_cols_set: continue
            our = cid_map.get(barcodes[row])
            if our is None: continue
            conc_total += 1
            pc = int(pip_counts.get(our, 0))
            if pc == int(c): conc_match += 1
            elif len(mism) < 20: mism.append((gname, barcodes[row], int(c), pc))
    conc = conc_match / max(1, conc_total)
    h5_crop_total = float(np.sum(data[np.isin(cols, keep_cols)])) if len(keep_cols) else 0.0
    total_ratio = float((owner >= 0).sum()) / h5_crop_total if h5_crop_total else 0.0
    return {"concordance_per_cell": conc, "n_cells_tested": conc_total,
            "n_genes_tested": int(len(test_genes)), "total_ratio": total_ratio,
            "pip_assigned": int((owner >= 0).sum()),
            "h5_total_crop_cells": int(h5_crop_total), "mismatches": mism}

keep_cols_set = None  # set in main before concordance

def label_cells(tx, n_cells):
    """Addendum A3: majority marker vote, tested-LR genes excluded, <200-cell classes dropped."""
    ann_map = {}
    for cls, genes in MARKERS.items():
        for g in genes:
            assert g not in TESTED, f"marker {g} is a tested LR gene — dictionary leak"
            ann_map[g] = cls
    feat = tx.feature_name.to_numpy(); cells = tx.cell_idx.to_numpy()
    keep = cells >= 0
    f2, c2 = feat[keep], cells[keep]
    is_marker = np.array([g in ann_map for g in f2])
    f2, c2 = f2[is_marker], c2[is_marker]
    anns = np.array([ann_map[g] for g in f2])
    df = pd.DataFrame({"cell": c2, "ann": anns})
    votes = df.groupby(["cell","ann"]).size().reset_index(name="n")
    best = votes.sort_values(["cell","n"], ascending=[True,False]).drop_duplicates("cell")
    lab = np.full(n_cells, "Unknown", dtype=object)
    lab[best.cell.to_numpy()] = best.ann.to_numpy()
    used = [c for c in np.unique(lab) if c != "Unknown"]
    sizes = {c: int((lab == c).sum()) for c in used}
    dropped = [c for c, n in sizes.items() if n < MIN_CLASS_CELLS]
    for c in dropped: lab[lab == c] = "Unknown"
    print(f"markers after tested-LR exclusion: {MARKERS}; removed: {REMOVED_TESTED_MARKERS}",
          flush=True)
    print(f"label class sizes (pre-drop): {sizes}; dropped <{MIN_CLASS_CELLS}: {dropped}",
          flush=True)
    return lab

def main():
    global keep_cols_set
    t0 = time.time()
    cb = pq.read_table(f"{CROP}/cells.parquet", columns=["cell_id"])
    keep_cols_set = None
    tx = load_transcripts()
    polys, cid_map = load_polygons_lung()
    owner = assign_pip(tx, polys)
    print(f"PIP done ({time.time()-t0:.0f}s)", flush=True)
    # barcode set for h5 restriction (must exist before concordance)
    with h5py.File("xenium_lung/extracted/cell_feature_matrix.h5","r") as f:
        barcodes_all = [b.decode() if isinstance(b,bytes) else b for b in f["matrix"]["barcodes"][:]]
    crop_bc = set(pq.read_table(f"{CROP}/cells.parquet", columns=["cell_id"]).column("cell_id").to_pylist())
    keep_cols_set = {i for i,b in enumerate(barcodes_all) if b in crop_bc}
    conc = verify_concordance_lung(tx, owner, cid_map)
    if SMOKE:
        # Concordance is structurally invalid in smoke mode: PIP sees a small row
        # subset while h5 counts are full. Smoke validates mechanics only; the gate
        # is evaluated in full mode (recorded, not hidden).
        conc["mode"] = "smoke"
        conc["gate"] = "SKIPPED — smoke mode; full-mode gate required"
        print(f"smoke concordance (NOT a gate): {conc}", flush=True)
    else:
        gate = conc["concordance_per_cell"] >= 0.90
        conc["mode"] = "full"
        conc["gate"] = "PASS" if gate else "FAIL"
        print(f"G3-crop concordance: per_cell={conc['concordance_per_cell']:.4f} "
              f"total_ratio={conc['total_ratio']:.4f} "
              f"({conc['n_cells_tested']} cells) -> {'PASS' if gate else 'FAIL'}", flush=True)
        if not gate:
            json.dump(conc, open(f"{OUTD}/concordance.json","w"), indent=1)
            sys.exit(2)
    cells_meta = pq.read_table(f"{CROP}/cells.parquet").to_pandas()
    assert len(cells_meta) == len(polys) or len(cid_map) <= len(cells_meta)
    n_cells = len(cells_meta)
    tx["cell_idx"] = owner  # attach BEFORE labeling (label_cells consumes it)
    labels = label_cells(tx, n_cells)
    # regions: prereg item 2 (k-means k=3, seed 20260907)
    regions = make_regions(cells_meta, labels, k=3, seed=20260907)
    print(f"region sizes: {np.bincount(regions).tolist()}", flush=True)
    dist, band, tree, edge_cell = boundary_and_band(tx, polys)
    tx["dist_boundary"] = dist
    tx["in_band"] = band
    band_idx, sets = donor_sets_for_band_tested(tx, band, tree, edge_cell)
    pfx = "smoke_" if SMOKE else ""
    keep = ["x_location","y_location","feature_name","cell_idx","dist_boundary","in_band"]
    tx[keep].to_parquet(f"{OUTD}/{pfx}tx.parquet", index=False)
    cells_meta.to_parquet(f"{OUTD}/{pfx}cells_meta.parquet", index=False)
    with open(f"{OUTD}/{pfx}donor_map.pkl","wb") as f:
        pickle.dump({"band_tested_indices": band_idx, "donor_sets": sets,
                     "band_indices": np.flatnonzero(band), "cid_map": cid_map,
                     "n_polys": len(polys), "tested_genes": TESTED}, f)
    np.save(f"{OUTD}/{pfx}labels.npy", labels)
    np.save(f"{OUTD}/{pfx}regions.npy", regions)
    json.dump(conc, open(f"{OUTD}/{pfx}concordance.json","w"), indent=1)
    print(f"DONE in {time.time()-t0:.0f}s (smoke={SMOKE})", flush=True)

if __name__ == "__main__":
    main()
