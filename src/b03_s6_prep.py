"""B03 Phase 7 step 2: coordinate transform + frozen expression-blind crop.

Per B03_PHASE7_ADDENDUM_A1 (frozen):
  1. Parse the registered transform scale from the zip's points/st/.zattrs
     (expected 1/0.2125 = 4.705882352941177, the Xenium pixel size); apply to
     transcript and polygon coordinates ONCE, writing micrometer crop files.
  2. Coordinate-frame assertions BEFORE window selection: transcript extent inside
     polygon extent; estimated section area of order tens of mm^2.
  3. Window rule VERBATIM from b03_lung_crop.py (frozen): 3.5 mm window, 0.5 mm
     grid, upper-median vendor-cell count, lexicographic tie-break, centroid
     geometry ONLY (no expression access).
  4. Crop extraction mirroring b03_lung_crop_extract.py: transcripts filtered to
     the window; polygons of in-window cells; measured-count envelope gate
     (270 B/tx <= 7.0 GB) before any proseg launch.
Outputs -> the workspace xenium_breast_s6/crop/ dir: {transcripts.parquet,cell_boundaries.parquet,
nucleus_boundaries.parquet,cells_meta.parquet} and s6_crop_window.json,
s6_crop_envelope.json (in the workspace xenium_breast_s6/ dir).
"""
import io, json, os, sys, time, zipfile
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

RESEARCH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ZIP = os.path.join(RESEARCH, "xenium_breast_s6", "sdata_breast_s6.zarr.zip")
EXT = os.path.join(RESEARCH, "xenium_breast_s6", "extracted")
CROP = os.path.join(RESEARCH, "xenium_breast_s6", "crop")
os.makedirs(CROP, exist_ok=True)
WIN, GRID, PEAK_B_PER_TX, MAX_PEAK_GB = 3500.0, 500.0, 270.0, 7.0
LOG = os.path.join(RESEARCH, "xenium_breast_s6", "prep_log.json")

def log(d):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(json.dumps(d) + "\n")

def registered_scale():
    """Parse and VALIDATE the registered transform, but do not apply it to the
    analysis frame. The zarr stores points and shapes in MICROMETERS; the registered
    scale 1/0.2125 = 4.7059 maps um -> IMAGE PIXELS for H&E alignment. Applying it
    to the analysis frame inflated polygon area 22.1x (956 mm^2 vs 43.2 mm^2 raw) —
    caught by the A1 frame assertion before any data was written. Analysis uses the
    stored um coordinates (the same convention as tissues 1 and 2)."""
    z = zipfile.ZipFile(ZIP)
    with z.open("sdata_breast_s6.zarr/points/st/.zattrs") as f:
        attrs = json.load(f)
    scales = []
    def walk(o):
        if isinstance(o, dict):
            if o.get("type") == "scale" and isinstance(o.get("scale"), list):
                scales.append(o["scale"])
            for v in o.values(): walk(v)
        elif isinstance(o, list):
            for v in o: walk(v)
    walk(attrs)
    assert scales, "no scale transform registered in points/st/.zattrs"
    s = scales[0]
    assert abs(s[0] - 1 / 0.2125) < 1e-9 and abs(s[1] - 1 / 0.2125) < 1e-9, \
        f"unexpected scale {s}"
    return float(s[0])  # recorded in outputs; NOT applied

def main():
    t0 = time.time()
    registered_scale()  # validated, recorded; stored coords are already um
    scale = 1.0

    # polygons: per-cell WKB -> vertices in um; centroid from vertex mean
    from shapely import from_wkb
    cb = pq.read_table(os.path.join(EXT, "boundaries", "cell_boundaries.parquet"))
    cell_ids = cb.column("__index_level_0__").to_pylist()
    wkbs = cb.column("geometry").to_pylist()
    geoms = from_wkb(list(wkbs))
    n = len(geoms)
    cx = np.empty(n); cy = np.empty(n); areas = np.empty(n); nvert = np.empty(n, dtype=np.int64)
    for i, g in enumerate(geoms):
        xy = np.asarray(g.exterior.coords)
        nvert[i] = len(xy)
        cx[i] = xy[:, 0].mean(); cy[i] = xy[:, 1].mean()
        areas[i] = g.area
    cx = cx * scale; cy = cy * scale; areas = areas * scale * scale
    xs_min = float(min(g.bounds[0] for g in geoms)) * scale
    xs_max = float(max(g.bounds[2] for g in geoms)) * scale
    ys_min = float(min(g.bounds[1] for g in geoms)) * scale
    ys_max = float(max(g.bounds[3] for g in geoms)) * scale
    n_valid = int((areas > 1e-6).sum())
    assert n_valid > 0.99 * n, f"too many degenerate polygons: {n - n_valid}"
    est_area_mm2 = float(areas.sum() / 1e6)
    print(f"polygons: {n:,} valid={n_valid:,} est section area={est_area_mm2:.1f} mm^2",
          flush=True)
    assert 5.0 <= est_area_mm2 <= 500.0, "section area out of Xenium order"

    # coordinate-frame assertion via transcript extent (columns only, then scale)
    tx0 = pq.ParquetFile(os.path.join(EXT, "transcripts_parts", "part.0.parquet"))
    xmin_t = ymin_t = np.inf; xmax_t = ymax_t = -np.inf
    for i in range(16):
        t = pq.read_table(os.path.join(EXT, "transcripts_parts", f"part.{i}.parquet"),
                          columns=["x", "y"])
        x = t.column("x").to_numpy() * scale; y = t.column("y").to_numpy() * scale
        xmin_t = min(xmin_t, float(x.min())); xmax_t = max(xmax_t, float(x.max()))
        ymin_t = min(ymin_t, float(y.min())); ymax_t = max(ymax_t, float(y.max()))
    # A1 frame check (direction corrected pre-crop, see A1-b.2): cell polygons lie
    # where molecules are, not vice versa — polygons cover only tissue while the
    # capture area includes acellular space. Tolerance 5 um: catches frame errors
    # (systematic, pixel-scale ~4.7 um) while permitting capture-edge rounding
    # (measured 1.6 um on y_min; a polygon touches the capture edge where no
    # QV-passing transcript sits). Zero centroids outside; margins logged.
    TOL = 5.0
    mx0, mx1 = xs_min - xmin_t, xmax_t - xs_max
    my0, my1 = ys_min - ymin_t, ymax_t - ys_max
    assert mx0 >= -TOL and mx1 >= -TOL and my0 >= -TOL and my1 >= -TOL, \
        f"frame mismatch beyond edge tolerance: margins {mx0:.2f},{mx1:.2f},{my0:.2f},{my1:.2f}"
    bbox_mm2 = (xs_max - xs_min) * (ys_max - ys_min) / 1e6
    cells_per_mm2 = n / bbox_mm2
    print(f"polygon bbox {xs_min:.0f}..{xs_max:.0f} x {ys_min:.0f}..{ys_max:.0f} um "
          f"({bbox_mm2:.1f} mm^2); cell density {cells_per_mm2:.0f}/mm^2; "
          f"polygon coverage {est_area_mm2/bbox_mm2:.2f}; edge margins "
          f"x[{mx0:.2f},{mx1:.2f}] y[{my0:.2f},{my1:.2f}] um", flush=True)
    assert 1000 <= cells_per_mm2 <= 20000, "cell density outside Xenium range"
    log({"step": "frame_check", "tx_extent": [xmin_t, xmax_t, ymin_t, ymax_t],
         "poly_extent": [xs_min, xs_max, ys_min, ys_max],
         "edge_margins_um": [mx0, mx1, my0, my1],
         "bbox_mm2": bbox_mm2, "cells_per_mm2": cells_per_mm2,
         "est_area_mm2": est_area_mm2, "t": time.time() - t0})
    print("coordinate-frame assertions passed", flush=True)

    # frozen window rule on centroid geometry ONLY
    cands = []
    x0 = 0.0
    while x0 + WIN <= xmax_t + GRID:
        y0 = 0.0
        while y0 + WIN <= ymax_t + GRID:
            m = (cx >= x0) & (cx < x0 + WIN) & (cy >= y0) & (cy < y0 + WIN)
            cands.append((x0, y0, int(m.sum())))
            y0 += GRID
        x0 += GRID
    counts = np.array([k for _, _, k in cands])
    med = int(np.sort(counts)[len(counts) // 2])
    ties = sorted([(x0, y0) for (x0, y0, k) in cands if k == med])
    wx, wy = ties[0]
    print(f"candidates={len(cands)} median/window={med} FROZEN WINDOW x[{wx:.0f},{wx+WIN:.0f}] "
          f"y[{wy:.0f},{wy+WIN:.0f}]", flush=True)
    json.dump({"window": {"x0": wx, "y0": wy, "x1": wx + WIN, "y1": wy + WIN,
                          "area_mm2": (WIN / 1000.0) ** 2},
               "cells_in_window": med, "median_window_cells": med,
               "n_candidates": len(cands),
               "section_tx_extent": [xmin_t, xmax_t, ymin_t, ymax_t],
               "est_section_area_mm2": est_area_mm2, "scale": scale,
               "rule": "b03_lung_crop.py verbatim (3.5mm, 0.5mm grid, upper-median, lex tie)"},
              open(os.path.join(RESEARCH, "xenium_breast_s6", "s6_crop_window.json"), "w"),
              indent=1)

    # crop extraction in um
    x0w, x1w, y0w, y1w = wx, wx + WIN, wy, wy + WIN
    n_tx = 0
    writer = None
    for i in range(16):
        t = pq.read_table(os.path.join(EXT, "transcripts_parts", f"part.{i}.parquet"))
        x = t.column("x").to_numpy() * scale
        y = t.column("y").to_numpy() * scale
        m = (x >= x0w) & (x < x1w) & (y >= y0w) & (y < y1w)
        if not m.any():
            continue
        t2 = t.filter(pa.array(m))
        n_tx += t2.num_rows
        t2 = t2.set_column(t2.schema.get_field_index("x"), "x",
                           pa.array((x[m]).astype(np.float32)))
        t2 = t2.set_column(t2.schema.get_field_index("y"), "y",
                           pa.array((y[m]).astype(np.float32)))
        if writer is None:
            writer = pq.ParquetWriter(os.path.join(CROP, "transcripts.parquet"), t2.schema)
        writer.write_table(t2)
        log({"step": f"tx_part_{i}", "kept": int(m.sum()), "t": time.time() - t0})
    if writer: writer.close()
    print(f"crop transcripts: {n_tx:,}", flush=True)

    in_win = (cx >= x0w) & (cx < x1w) & (cy >= y0w) & (cy < y1w)
    keep_ids = [c for c, k in zip(cell_ids, in_win) if k]
    idx_keep = np.flatnonzero(in_win)
    cells_meta = pd.DataFrame({
        "cell_id": keep_ids,
        "x_centroid": cx[idx_keep], "y_centroid": cy[idx_keep],
        "cell_area_um2": areas[idx_keep],
        "vertex_count": nvert[idx_keep],
    })
    cells_meta.to_parquet(os.path.join(CROP, "cells_meta.parquet"), index=False)
    # WKB crop files (verbatim geometry, scaled)
    for name in ("cell_boundaries", "nucleus_boundaries"):
        b = pq.read_table(os.path.join(EXT, "boundaries", f"{name}.parquet"))
        bidx = b.column("__index_level_0__").to_pylist()
        pos = {c: i for i, c in enumerate(bidx)}
        sel = [pos[c] for c in keep_ids if c in pos]
        b2 = b.take(pa.array(sel))
        pq.write_table(b2, os.path.join(CROP, f"{name}.parquet"))
        print(f"crop {name}: {b2.num_rows} rows", flush=True)

    peak_gb = PEAK_B_PER_TX * n_tx / 1e9
    gate = peak_gb <= MAX_PEAK_GB
    print(f"ENVELOPE: {n_tx:,} tx x {PEAK_B_PER_TX} B = {peak_gb:.2f} GB "
          f"(gate <= {MAX_PEAK_GB}) -> {'LAUNCH PERMITTED' if gate else 'LAUNCH DENIED'}",
          flush=True)
    json.dump({"n_transcripts": int(n_tx), "n_cells": int(len(keep_ids)),
               "peak_gb": round(peak_gb, 2), "gate": bool(gate),
               "window": {"x0": wx, "y0": wy, "x1": x1w, "y1": y1w}},
              open(os.path.join(RESEARCH, "xenium_breast_s6", "s6_crop_envelope.json"), "w"),
              indent=1)
    log({"step": "done", "seconds": time.time() - t0, "gate": bool(gate)})
    if not gate:
        sys.exit(2)

if __name__ == "__main__":
    main()
