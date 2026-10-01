"""B03 Phase 5 addendum A2 step 2: extract the frozen crop, count transcripts, envelope gate.

Reads the frozen window from results/lung_crop_window.json. Writes:
  <workspace>/xenium_lung/crop/transcripts.parquet  (rows with window <= x,y < window)
  <workspace>/xenium_lung/crop/cell_boundaries.parquet, nucleus_boundaries.parquet, cells.parquet
  (vendor cells whose centroid is inside the window; boundaries matched by cell key)
Envelope gate (frozen): peak = 270 B/tx x n_tx must be <= 7.0 GB else no proseg launch.
"""
import json, os, sys
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
CROP = os.path.join(WORKSPACE, "xenium_lung", "crop")
EXT = os.path.join(WORKSPACE, "xenium_lung", "extracted")
env = json.load(open(os.path.join(ROOT, "results", "lung_crop_window.json")))
w = env["window"]
x0, x1, y0, y1 = w["x0"], w["x1"], w["y0"], w["y1"]
os.makedirs(CROP, exist_ok=True)

# --- transcripts: stream row groups, keep in-window ---
src = pq.ParquetFile(os.path.join(EXT, "transcripts.parquet"))
writer = None
n_tx = 0
for i in range(src.num_row_groups):
    t = src.read_row_group(i)
    x = pc.and_(pc.greater_equal(t.column("x_location"), x0),
                pc.less(t.column("x_location"), x1))
    y = pc.and_(pc.greater_equal(t.column("y_location"), y0),
                pc.less(t.column("y_location"), y1))
    m = pc.and_(x, y)
    if pc.all(m).as_py():          # whole group inside: copy as-is
        keep = t
    else:
        keep = t.filter(m)
    if keep.num_rows == 0:
        continue
    n_tx += keep.num_rows
    if writer is None:
        writer = pq.ParquetWriter(os.path.join(CROP, "transcripts.parquet"), keep.schema)
    writer.write_table(keep)
if writer: writer.close()
print(f"crop transcripts: {n_tx:,}")

# --- cells metadata: filter to window ---
c = pq.read_table(os.path.join(EXT, "cells.parquet"))
cm = pc.and_(pc.and_(pc.greater_equal(c.column("x_centroid"), x0),
                     pc.less(c.column("x_centroid"), x1)),
             pc.and_(pc.greater_equal(c.column("y_centroid"), y0),
                     pc.less(c.column("y_centroid"), y1)))
c_crop = c.filter(cm)
pq.write_table(c_crop, os.path.join(CROP, "cells.parquet"))
print("crop cells:", c_crop.num_rows)

# --- boundaries: keep polygons of in-window cells ---
cell_ids = set(c_crop.column("cell_id").to_pylist())
for name in ("cell_boundaries", "nucleus_boundaries"):
    b = pq.read_table(os.path.join(EXT, f"{name}.parquet"))
    bcol = b.column("cell_id") if "cell_id" in b.column_names else b.column(0)
    mask = pa.array([v in cell_ids for v in bcol.to_pylist()])
    b_crop = b.filter(mask)
    pq.write_table(b_crop, os.path.join(CROP, f"{name}.parquet"))
    print(f"crop {name}: {b_crop.num_rows} rows")

# --- envelope gate ---
peak_gb = 270.0 * n_tx / 1e9
print(f"ENVELOPE: predicted proseg peak = {peak_gb:.2f} GB (gate <= {7.0} GB) -> "
      + ("LAUNCH PERMITTED" if peak_gb <= 7.0 else "LAUNCH DENIED"))
json.dump(dict(n_transcripts=int(n_tx), n_cells=c_crop.num_rows,
               peak_gb=round(peak_gb, 2), gate=peak_gb <= 7.0, window=w),
          open(os.path.join(ROOT, "results", "lung_crop_envelope.json"), "w"), indent=1)
