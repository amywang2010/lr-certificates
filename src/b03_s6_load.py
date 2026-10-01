"""B03 Phase 7: breast_s6 crop vendor load — PIP, concordance, band geometry, labels, regions.

Reuses certified machinery BY IMPORT (b03_load: assign_pip, boundary_and_band,
donor_sets_for_band_tested, BAND_UM; b03_certify: make_regions, seed 20260907).
s6-specific I/O per B03_PHASE7_ADDENDUM_A1 + A1-b (frozen before any compute):
  - transcripts: <workspace>/xenium_breast_s6/crop/transcripts.parquet (frozen window
    x[7500,11000] x y[6500,10000], frozen crop rule verbatim on centroid density).
  - polygons:    crop cell_boundaries.parquet as per-cell WKB (simpler than the 10x
    exploded format; same validity filter: area > 1e-6, buffer(0) repair, largest part).
  - concordance: A1-c totals-level gate — vendor obs/transcript_counts (polygon id
    space, raw int64) vs our PIP over is_gene & QV>=20 transcripts; gates r >= 0.99,
    rho >= 0.99, median |diff| <= 5, ratio in [0.90, 1.10]. The packaging ships no
    raw count matrix on the polygon id space (A1-c), so the per-gene gate form of
    tissues 1/2 is structurally unavailable here; this gate was frozen before any
    gene-level computation.
  - labels:      A3 frozen marker dictionary, majority vote, tested-LR genes excluded;
    classes < 200 cells dropped (frozen minimum-class rule).
  - regions:     k-means k=3 on (log10(1+transcript_counts), x/1000, y/1000), seed
    20260907. transcript_counts joined from the zip obs (vendor per-cell total),
    asserted to match the CSR row totals exactly.

Outputs -> the s6 crop data dir: tx.parquet, donor_map.pkl, cells_meta.parquet,
concordance.json, labels.npy, regions.npy  (lung-loader schema, downstream-compatible).
Smoke mode (B03_S6_SMOKE=1): first 300k transcripts, smoke_* artifacts, all
assertions still run; concordance recorded but not gated (structural mismatch of
row subsets, same policy as the lung loader).
"""
import os, sys, json, pickle, time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_load
from b03_load import assign_pip, boundary_and_band, donor_sets_for_band_tested, BAND_UM
from b03_certify import make_regions
from shapely import from_wkb
from shapely.geometry import Polygon

SMOKE = os.environ.get("B03_S6_SMOKE") == "1"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
CROP = os.path.join(WORKSPACE, "xenium_breast_s6", "crop")
TC = os.path.join(WORKSPACE, "xenium_breast_s6", "extracted", "table_cells")
OUTD = f"{CROP}/data"
os.makedirs(OUTD, exist_ok=True)
LOG = os.path.join(WORKSPACE, "xenium_breast_s6", "load_log.json")

def log(d):
    d = {k: (v.item() if hasattr(v, "item") else v) for k, v in d.items()}
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(json.dumps(d) + "\n")

# Frozen tested-gene set = unique genes of the 18 testable pairs (A3 config)
PAIRS = [("CD274","PDCD1"),("CXCL12","CXCR4"),("PECAM1","KDR"),
 ("ERBB2","EGFR"),("ERBB2","PDCD1"),("ESR1","PGR"),("PGR","ESR1"),
 ("CD274","CTLA4"),("CD274","CD8A"),("CXCL12","CCR7"),("PECAM1","PDCD1"),
 ("PECAM1","CTLA4"),("CD163","CCR7"),("CD163","CXCR4"),("CDH1","EGFR"),
 ("CDH1","ERBB2"),("MS4A1","CD274"),("CD68","CD274")]
TESTED = sorted({g for p in PAIRS for g in p})
assert len(TESTED) == 17, TESTED

# A3 frozen marker dictionary; tested-LR genes removed HERE (assertion verifies
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
MIN_CLASS_CELLS = 200  # minimum-class rule, frozen

def decode_zarr_string_dir(d):
    """Decode a zarr v2 string array extracted as chunk files (numcodecs direct)."""
    from numcodecs import blosc, VLenUTF8
    meta = json.load(open(f"{d}/.zarray.json"))
    shape = int(meta["shape"][0]) if isinstance(meta["shape"], (list, tuple)) else int(meta["shape"])
    files = sorted([f for f in os.listdir(d) if f.split(".")[0].isdigit()],
                   key=lambda s: int(s))
    out = []
    for f in files:
        out.extend(VLenUTF8().decode(blosc.decompress(open(f"{d}/{f}", "rb").read())))
    assert len(out) >= shape, (len(out), shape)
    return [None if v == "" else v for v in out[:shape]]

def decode_zarr_int_dir(d):
    from numcodecs import blosc
    meta = json.load(open(f"{d}/.zarray.json"))
    shape = int(meta["shape"][0])
    dt = np.dtype(meta["dtype"])
    files = sorted([f for f in os.listdir(d) if f.split(".")[0].isdigit()],
                   key=lambda s: int(s))
    out = []
    for f in files:
        out.append(np.frombuffer(blosc.decompress(open(f"{d}/{f}", "rb").read()), dtype=dt))
    v = np.concatenate(out)
    assert len(v) >= shape, (len(v), shape)
    return v[:shape]

def load_transcripts():
    f = pq.ParquetFile(f"{CROP}/transcripts.parquet")
    frames = []
    n = 0
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=["x", "y", "feature_name", "qv", "is_gene"]).to_pandas()
        t = t.rename(columns={"x": "x_location", "y": "y_location"})
        if SMOKE and n >= 300_000: break
        frames.append(t); n += len(t)
    tx = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if SMOKE: tx = tx.head(300_000)
    tx = tx[tx.qv >= 20.0].reset_index(drop=True)
    print(f"transcripts (QV>=20): {len(tx):,}", flush=True)
    return tx

def load_polygons_s6():
    """Per-cell WKB -> shapely polygons; same validity filter as tissues 1 and 2."""
    cb = pq.read_table(f"{CROP}/cell_boundaries.parquet", columns=["geometry"]).to_pandas()
    geoms = from_wkb(cb["geometry"].to_numpy())
    polys, cid_map = {}, {}
    cells_meta = pd.read_parquet(f"{CROP}/cells_meta.parquet")
    assert len(geoms) == len(cells_meta), (len(geoms), len(cells_meta))
    n_invalid = 0
    for k in range(len(geoms)):
        p = geoms[k]
        if p is None or p.is_empty:
            n_invalid += 1; continue
        if not (p.is_valid and p.area > 1e-6):
            p = p.buffer(0)
            if p.is_empty or p.area <= 1e-6:
                n_invalid += 1; continue
        if p.geom_type == "MultiPolygon":
            p = max(p.geoms, key=lambda q: q.area)
        if p.area <= 1e-6:
            n_invalid += 1; continue
        cid = cells_meta.cell_id.iloc[k]
        cid_map[cid] = len(polys); polys[len(polys)] = p
    print(f"valid crop polygons: {len(polys):,} (invalid dropped: {n_invalid})", flush=True)
    return polys, cid_map, cells_meta

def verify_totals_s6(tx, owner, cid_map, cm):
    """A1-c totals-level G3 gate (frozen before any gene-level computation).

    Reference: vendor obs/transcript_counts per crop cell (polygon id space).
    Candidate: our PIP over is_gene & QV>=20 crop transcripts.
    Gates: Pearson r >= 0.99 AND Spearman rho >= 0.99 AND median |diff| <= 5
           AND total ratio in [0.90, 1.10].
    Rationale (A1-c): this packaging ships no raw count matrix on the polygon id
    space (table_cells.X is processed; points cell_id is a different segmentation
    layer), so the per-gene gate form of tissues 1/2 is structurally unavailable.
    """
    from scipy.stats import spearmanr
    is_gene = tx["is_gene"].to_numpy().astype(bool)
    o = owner[is_gene]
    cnt = pd.Series(o[o >= 0]).value_counts()
    pip_totals = np.array([int(cnt.get(i, 0)) for i in range(len(cm))])
    obs_ids_all = decode_zarr_string_dir(f"{TC}/obs_cell_id")
    tx_counts_all = decode_zarr_int_dir(f"{TC}/obs_transcript_counts")
    lut = dict(zip(obs_ids_all, tx_counts_all.tolist()))
    obs_totals = np.array([lut[c] for c in cm.cell_id])
    r = float(np.corrcoef(pip_totals, obs_totals)[0, 1])
    rs = float(spearmanr(pip_totals, obs_totals).statistic)
    med_abs = float(np.median(np.abs(pip_totals - obs_totals)))
    ratio = float(pip_totals.sum()) / float(obs_totals.sum())
    gate = (r >= 0.99) and (rs >= 0.99) and (med_abs <= 5) and (0.90 <= ratio <= 1.10)
    return {"gate_form": "A1-c totals-level", "pearson_r": r, "spearman_rho": rs,
            "median_abs_diff": med_abs, "total_ratio": ratio,
            "n_cells": int(len(cm)), "pip_assigned_total": int(pip_totals.sum()),
            "vendor_total": int(obs_totals.sum()),
            "gate": "PASS" if gate else "FAIL"}

def label_cells(tx, n_cells):
    """A3: majority marker vote, tested-LR genes excluded, <200-cell classes dropped."""
    ann_map = {}
    for cls, genes in MARKERS.items():
        for g in genes:
            assert g not in TESTED, f"marker {g} is a tested LR gene - dictionary leak"
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
    t0 = time.time()
    tx = load_transcripts()
    polys, cid_map, cells_meta = load_polygons_s6()
    owner = assign_pip(tx, polys)
    print(f"PIP done ({time.time()-t0:.0f}s)", flush=True)
    log({"step": "pip", "tx": len(tx), "polys": len(polys), "t": time.time()-t0})
    conc = verify_totals_s6(tx, owner, cid_map, cells_meta)
    if SMOKE:
        conc["mode"] = "smoke"
        conc["gate"] = "SKIPPED - smoke mode; full-mode gate required"
        print(f"smoke totals-gate (NOT a gate): {conc}", flush=True)
    else:
        print(f"G3 totals-gate (A1-c): r={conc['pearson_r']:.5f} "
              f"rho={conc['spearman_rho']:.5f} med|d|={conc['median_abs_diff']:.1f} "
              f"ratio={conc['total_ratio']:.5f} -> {conc['gate']}", flush=True)
        json.dump(conc, open(f"{OUTD}/concordance.json", "w"), indent=1)
        if conc["gate"] != "PASS":
            sys.exit(2)
    # vendor per-cell totals: join from the zip obs by cell-id string, then assert
    # the join is total on the crop (every crop cell found exactly once).
    obs_ids_all = decode_zarr_string_dir(f"{TC}/obs_cell_id")
    tx_counts_all = decode_zarr_int_dir(f"{TC}/obs_transcript_counts")
    lut = {}
    for r_, b in enumerate(obs_ids_all):
        if b is not None and b in cid_map:
            lut[b] = int(tx_counts_all[r_])
    assert set(cells_meta.cell_id) <= set(lut), "crop cell missing from vendor obs"
    cells_meta = cells_meta.copy()
    cells_meta["transcript_counts"] = [lut[c] for c in cells_meta.cell_id]
    print(f"vendor totals joined: {len(cells_meta):,} cells", flush=True)
    n_cells = len(cells_meta)
    tx["cell_idx"] = owner  # attach BEFORE labeling
    labels = label_cells(tx, n_cells)
    regions = make_regions(cells_meta, labels, k=3, seed=20260907)
    print(f"region sizes: {np.bincount(regions).tolist()}", flush=True)
    dist, band, tree, edge_cell = boundary_and_band(tx, polys)
    tx["dist_boundary"] = dist
    tx["in_band"] = band
    # A3 registered mover filter: donor sets cover band transcripts of the 17
    # tested genes ONLY. b03_load's module default is tissue-1's 25-pair superset;
    # the lung leg was unaffected because its panel lacks every superset-only gene
    # (verified), but s6's panel carries CCL5/CD3D. Override explicitly and assert.
    b03_load.TESTED_GENES = list(TESTED)
    assert set(b03_load.TESTED_GENES) == set(TESTED) and len(TESTED) == 17
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
    log({"step": "done", "tx": len(tx), "cells": n_cells, "t": time.time()-t0})
    print(f"DONE in {time.time()-t0:.0f}s (smoke={SMOKE})", flush=True)

if __name__ == "__main__":
    main()
