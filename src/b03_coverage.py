"""B03 Step 5: held-out coverage check (pre-registered) + nucleus sensitivity report.

Design decisions (frozen before expression access):

1. COVERAGE VARIANT = 2 um-eroded cell masks (design freeze: "Proseg or erosion variant as
   held-out check"). Rationale: U consists of reassignments of the SAME molecules among
   the SAME cell inventory. The eroded-mask assignment preserves the inventory (every
   cell persists; only its transcript set changes) and is therefore a member of U by
   construction -- the check certifies that the certificate machinery covers a concrete,
   independently implemented alternative assignment.
   - Base polygons are imported from b03_load.load_polygons() (guarantees polygon-index
     alignment with frozen labels); the ASSIGNMENT is a fresh code path: PIP of all
     transcripts onto the eroded geometry. No donor/band/certificate machinery is used.
   - Proseg is NOT a valid in-U variant (it merges/splits cells -> changes the cell
     inventory), so it cannot serve this gate; it is deferred to the cross-platform
     phase of the full campaign (requires a Rust toolchain install (one-time environment setup)).

2. NUCLEUS-ONLY VARIANT = sensitivity report only (OUTSIDE U, not gated): nucleus PIP
   removes cytoplasmic transcripts wholesale; T(A_nucleus) quantifies how much of each
   contrast survives when cells are shrunk to nuclei. Reported alongside [T_lo, T_hi].
   Cells without a valid nucleus are given an empty polygon (all their transcripts
   unassigned). Fixed cell-count convention (inventory preserved) disclosed; conservative.

Coverage criterion (design freeze): >= 90% of the 75 (pair, region) rows must satisfy
T_lo <= T(eroded) <= T_hi.
"""
import sys
import numpy as np
import pandas as pd
import shapely
from shapely import STRtree
import json, time, math

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from b03_load import load_polygons  # exact same polygon construction/order as A0

EPS = 0.5
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")
ERODE_UM = 2.0
NUC_FILE = f"{DATA}/Xenium_FFPE_Human_Breast_Cancer_Rep1_nucleus_boundaries.parquet"

POSITIVE_CONTROLS = [("CD274", "PDCD1"), ("CXCL12", "CXCR4"), ("CCL5", "CCR7"),
                     ("PECAM1", "KDR"), ("CD163", "CD3D")]
DISPUTED = [("ERBB2", "EGFR"), ("ERBB2", "PDCD1"), ("ESR1", "PGR"), ("PGR", "ESR1"),
            ("CD274", "CTLA4"), ("CD274", "CD8A"), ("CXCL12", "CCR7"), ("CCL5", "CXCR4"),
            ("PECAM1", "PDCD1"), ("PECAM1", "CTLA4"), ("CD163", "CCR7"), ("CD163", "CXCR4"),
            ("ACTA2", "EGFR"), ("ACTA2", "ERBB2"), ("KRT8", "ERBB2"), ("CDH1", "EGFR"),
            ("CDH1", "ERBB2"), ("CD3D", "ERBB2"), ("MS4A1", "CD274"), ("CD68", "CD274")]
ALL_PAIRS = POSITIVE_CONTROLS + DISPUTED
TYPE_MAP = {"CD274": "Breast cancer", "CXCL12": "Fibroblasts", "CCL5": "Macrophages",
            "PECAM1": "Endothelial cells", "CD163": "Macrophages", "ERBB2": "Breast cancer",
            "ESR1": "Breast cancer", "PGR": "Breast cancer", "ACTA2": "Smooth muscle cells",
            "KRT8": "Breast glandular cells", "CDH1": "Breast glandular cells",
            "CD3D": "T cells", "MS4A1": "B cells", "CD68": "Macrophages"}
RECV_MAP = {"PDCD1": "T cells", "CXCR4": "T cells", "CCR7": "T cells",
            "KDR": "Endothelial cells", "CD3D": "T cells", "EGFR": "Breast cancer",
            "PGR": "Breast cancer", "ESR1": "Breast cancer", "CTLA4": "T cells",
            "CD8A": "T cells", "ERBB2": "Breast cancer", "CD274": "Breast cancer"}


def t_log2(n_l, n_r, n_l_cells, n_r_cells):
    return (math.log2(max(n_l, 0) / n_l_cells + EPS) +
            math.log2(max(n_r, 0) / n_r_cells + EPS))


def assign_pip_fresh(tx_xy, poly_list):
    """Fresh-code-path PIP (independent of the loader's assign_pip): deterministic
    first-polygon (smallest index) assignment for multiply-contained transcripts."""
    tree = STRtree(np.array(poly_list, dtype=object))
    pts = shapely.points(tx_xy[:, 0], tx_xy[:, 1])
    pairs = tree.query(pts, predicate="intersects")
    tx_idx, poly_idx = pairs[0], pairs[1]
    ok = (tx_idx >= 0) & (tx_idx < len(tx_xy)) & (poly_idx >= 0) & (poly_idx < len(poly_list))
    tx_idx, poly_idx = tx_idx[ok], poly_idx[ok]
    order = np.lexsort((poly_idx, tx_idx))
    tx2, po2 = tx_idx[order], poly_idx[order]
    first = np.ones(len(tx2), dtype=bool)
    first[1:] = tx2[1:] != tx2[:-1]
    cell_of = np.full(len(tx_xy), -1, dtype=np.int64)
    cell_of[tx2[first]] = po2[first]
    return cell_of


def load_nucleus_polys_aligned(cid_map, n_out):
    """Nucleus polygons aligned to CELL polygon index space (via cid_map)."""
    import pyarrow.parquet as pq
    nb = pq.read_table(NUC_FILE, columns=["vertex_x", "vertex_y", "cell_id"]).to_pandas()
    nb = nb.sort_values("cell_id", kind="stable").reset_index(drop=True)
    nids = nb.cell_id.to_numpy()
    starts = np.flatnonzero(np.r_[True, nids[1:] != nids[:-1]])
    ends = np.r_[starts[1:], len(nids)]
    xs_all = nb.vertex_x.to_numpy(np.float64)
    ys_all = nb.vertex_y.to_numpy(np.float64)
    nuc_by_cid = {}
    for k in range(len(starts)):
        s, e = starts[k], ends[k]
        c = nids[s]
        xs, ys = xs_all[s:e], ys_all[s:e]
        if len(xs) < 3:
            continue
        p = shapely.Polygon(np.stack([xs, ys], axis=1)) if hasattr(shapely, "polygons") \
            else None
        from shapely.geometry import Polygon
        p = Polygon(zip(xs, ys))
        if not (p.is_valid and p.area > 1e-6):
            p = p.buffer(0)
            if p.is_empty or p.area <= 1e-6:
                continue
        if p.geom_type == "MultiPolygon":
            p = max(p.geoms, key=lambda q: q.area)
        nuc_by_cid[c] = p
    empty = shapely.Polygon()
    out = []
    for c, i in cid_map.items():
        while len(out) < i:
            out.append(empty)
        out.append(nuc_by_cid.get(c, empty))
    while len(out) < n_out:
        out.append(empty)
    return out[:n_out]


def main():
    t0 = time.time()
    tx = pd.read_parquet(f"{DATA}/tx.parquet")
    cells_meta = pd.read_parquet(f"{DATA}/cells_meta.parquet")
    labels = np.load(f"{DATA}/labels.npy", allow_pickle=True)
    region_of_cell = np.load(f"{DATA}/regions.npy")
    res = pd.read_csv(f"{OUT}/scout_bounds.csv")
    assert len(cells_meta) == len(labels) == len(region_of_cell)
    tested = sorted({g for p in ALL_PAIRS for g in p})

    xy = np.c_[tx.x_location.to_numpy(np.float64), tx.y_location.to_numpy(np.float64)]
    feat = tx.feature_name.to_numpy()
    print(f"[{time.time()-t0:.0f}s] transcripts {len(tx):,}", flush=True)

    # ---- base polygons (loader semantics -> index alignment with labels) ----
    polys_dict, cid_map = load_polygons()
    n_cells = len(cells_meta)
    assert len(polys_dict) == n_cells, \
        f"valid polygons {len(polys_dict)} != cells {n_cells}"
    poly_list = [polys_dict[i] for i in range(len(polys_dict))]

    # ---- variant 1: eroded cell masks (in-U coverage check) ----
    eroded = shapely.buffer(np.array(poly_list, dtype=object), -ERODE_UM)
    eroded = shapely.make_valid(eroded)
    cell_idx_er = assign_pip_fresh(xy, list(eroded))
    print(f"[{time.time()-t0:.0f}s] eroded PIP done: assigned "
          f"{int((cell_idx_er >= 0).sum()):,}", flush=True)

    # ---- variant 2: nucleus masks (outside-U sensitivity) ----
    nuc_list = load_nucleus_polys_aligned(cid_map, n_cells)
    cell_idx_nu = assign_pip_fresh(xy, nuc_list)
    print(f"[{time.time()-t0:.0f}s] nucleus PIP done: assigned "
          f"{int((cell_idx_nu >= 0).sum()):,}", flush=True)

    n_slots = len(set(TYPE_MAP.values()) | set(RECV_MAP.values())) + 1
    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    type_shift = {t: i + 1 for i, t in enumerate(types_needed)}
    cid = region_of_cell * n_slots + np.array([type_shift.get(l, 0) for l in labels],
                                              dtype=np.int64)

    rows = []
    for pair in ALL_PAIRS:
        Lg, Rg = pair
        S, R = TYPE_MAP[Lg], RECV_MAP[Rg]
        for region_id in range(3):
            rr = res[(res.pair == f"{Lg}->{Rg}") & (res.region == region_id)]
            if rr.empty:
                continue
            rr = rr.iloc[0]
            nLc, nRc = int(rr.nl_cells), int(rr.nr_cells)
            cid_S = region_id * n_slots + type_shift[S]
            cid_R = region_id * n_slots + type_shift[R]
            out = dict(pair=f"{Lg}->{Rg}", region=region_id,
                       control=pair in POSITIVE_CONTROLS)
            for tag, cell_idx in (("eroded", cell_idx_er), ("nucleus", cell_idx_nu)):
                valid = cell_idx >= 0
                cid_of_tx = np.where(valid, cid[np.clip(cell_idx, 0, None)], -1)
                maskL = (feat == Lg) & (cid_of_tx == cid_S)
                maskR = (feat == Rg) & (cid_of_tx == cid_R)
                N_L, N_R = int(maskL.sum()), int(maskR.sum())
                out[f"N_L_{tag}"] = N_L
                out[f"N_R_{tag}"] = N_R
                out[f"T_{tag}"] = t_log2(N_L, N_R, nLc, nRc)
            out["T_lo"] = float(rr.T_lo)
            out["T_hi"] = float(rr.T_hi)
            out["T0"] = float(rr.T0)
            out["covered_eroded"] = bool(rr.T_lo - 1e-9 <= out["T_eroded"] <=
                                         rr.T_hi + 1e-9)
            out["nucleus_outside_U"] = bool((out["T_nucleus"] < rr.T_lo - 1e-9) or
                                            (out["T_nucleus"] > rr.T_hi + 1e-9))
            rows.append(out)
        print(f"[{time.time()-t0:.0f}s] {Lg}->{Rg} done", flush=True)

    df = pd.DataFrame(rows)
    cov = float(df.covered_eroded.mean())
    df.to_csv(f"{OUT}/coverage_check.csv", index=False)
    summary = {
        "coverage_variant": f"cell masks eroded by {ERODE_UM} um (in-U, fresh assignment path)",
        "coverage_fraction": cov,
        "coverage_gate_90pct": bool(cov >= 0.90),
        "n_rows": len(df),
        "n_covered": int(df.covered_eroded.sum()),
        "nucleus_variant": "sensitivity only (outside U)",
        "nucleus_outside_U_rows": int(df.nucleus_outside_U.sum()),
    }
    with open(f"{OUT}/coverage_summary.json", "w") as fj:
        json.dump(summary, fj, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"DONE in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
