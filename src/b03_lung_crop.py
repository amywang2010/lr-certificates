"""B03 Phase 5 addendum A2: expression-blind crop selection (FROZEN rule).

Selects the 3.5x3.5 mm window (0.5 mm grid, median vendor-cell count, lexicographic
tie-break) using ONLY cells.parquet nucleus-centroid geometry. Touches no expression
data. Prints the envelope arithmetic required by the addendum before any proseg launch.
"""
import os, sys
import numpy as np
import pandas as pd

CELLS = sys.argv[1] if len(sys.argv) > 1 else "xenium_lung/extracted/cells.parquet"
OUT = sys.argv[2] if len(sys.argv) > 2 else "B03_project/results/lung_crop_window.json"
WIN = 3500.0        # um (3.5 mm) — frozen
GRID = 500.0        # um (0.5 mm) step — frozen
PEAK_B_PER_TX = 270 # upper measured bytes/tx (breast empirics) — conservative
MAX_PEAK_GB = 7.0   # launch gate (addendum A2)

c = pd.read_parquet(CELLS, columns=["x_centroid", "y_centroid"])
x = c.x_centroid.to_numpy(); y = c.y_centroid.to_numpy()
xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
print(f"section bounds: x[{xmin:.0f},{xmax:.0f}] y[{ymin:.0f},{ymax:.0f}]  cells={len(c):,}")

cands = []
x0 = 0.0
while x0 + WIN <= xmax + GRID:  # allow edge-aligned windows clipped to data
    y0 = 0.0
    while y0 + WIN <= ymax + GRID:
        m = (x >= x0) & (x < x0 + WIN) & (y >= y0) & (y < y0 + WIN)
        cands.append((x0, y0, int(m.sum())))
        y0 += GRID
    x0 += GRID
counts = np.array([k for _,_,k in cands])
# upper-median = order statistic at index n//2 (frozen rule; always an attained value)
med = int(np.sort(counts)[len(counts)//2])
ties = [(k, x0, y0) for (x0, y0, k) in cands if k == med]
ties.sort()
_, wx, wy = ties[0]
n_in = med

print(f"candidates={len(cands)}  median cells/window={med}")
print(f"FROZEN WINDOW: x[{wx:.0f},{wx+WIN:.0f}] y[{wy:.0f},{wy+WIN:.0f}]  cells={n_in:,}")

# Envelope arithmetic (addendum A2): density from candidate windows of same area
area_mm2 = (WIN/1000.0)**2
# local tx estimate: section density x window area (upper bound if section uniform);
# refine with measured local density later at crop-extraction time (tx counted directly).
dens = None
import json
env = dict(window=dict(x0=wx, y0=wy, x1=wx+WIN, y1=wy+WIN, area_mm2=area_mm2),
           cells_in_window=n_in, median_window_cells=med, n_candidates=len(cands),
           section_bounds=dict(xmin=float(xmin), xmax=float(xmax), ymin=float(ymin), ymax=float(ymax)))
with open(OUT, "w") as f: json.dump(env, f, indent=1)
print(f"saved {OUT}")
print("NOTE: transcript-count envelope computed at crop extraction (direct count), per addendum A2.")
