"""Pseudocount sensitivity of the interval signs (S5).

The score uses eps = 0.5. This script recomputes T_lo and T_hi from the stored
count endpoints at eps in {0.1, 0.5, 1.0} and counts positive, negative and
zero-crossing intervals for each section. The 0.5 row is checked against the
stored endpoints before the table is written.

Output: sensitivity/eps_sensitivity.csv
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EPS_GRID = (0.1, 0.5, 1.0)
SECTIONS = (
    ("breast", os.path.join(ROOT, "results", "breast", "scout_final_matched.csv")),
    ("lung", os.path.join(ROOT, "results", "lung", "scout_final_lung.csv")),
    ("breast_s6", os.path.join(ROOT, "results", "breast_s6", "scout_final_s6.csv")),
)


def endpoints(df, eps):
    lo = (np.log2((df.nL + df.dL_min) / df.nl_cells + eps)
          + np.log2((df.nR + df.dR_min) / df.nr_cells + eps))
    hi = (np.log2((df.nL + df.dL_max) / df.nl_cells + eps)
          + np.log2((df.nR + df.dR_max) / df.nr_cells + eps))
    return lo.to_numpy(float), hi.to_numpy(float)


def main():
    rows = []
    for name, path in SECTIONS:
        df = pd.read_csv(path)
        for eps in EPS_GRID:
            lo, hi = endpoints(df, eps)
            if eps == 0.5:
                dev = max(float(np.abs(lo - df.T_lo).max()),
                          float(np.abs(hi - df.T_hi).max()))
                assert dev < 1e-12, (name, dev)
            rows.append(dict(
                section=name, eps=eps,
                positive=int((lo > 0).sum()),
                negative=int((hi < 0).sum()),
                crossing=int(((lo <= 0) & (hi >= 0)).sum()),
                rows=len(df)))
    out = pd.DataFrame(rows)
    os.makedirs(os.path.join(ROOT, "sensitivity"), exist_ok=True)
    out.to_csv(os.path.join(ROOT, "sensitivity", "eps_sensitivity.csv"), index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
