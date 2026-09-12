"""B03 Phase 6: lung B=10,000 matched-null driver (B03_PHASE6_PREREG.md, frozen 2026-09-08).

Mirrors the sealed Phase 5 driver (b03_lung_null.py) exactly: same seed
(20260907), same data dir, same A3-frozen pair list and type maps from the
shared config module, same B03_NULL_CAP=1. Differences per prereg: B=10000 and
a fresh OUT directory so the sealed B=1000 artifacts remain untouched.
Per-perm seeds 20260907+100003+p: perms 0..999 are bitwise identical to the
B=1000 run; the chain script enforces the exact-match determinism gate against
null_v3_checkpoint.npz before any output is used.
"""
import os
import sys

os.environ["B03_NULL_SEED"] = "20260907"
os.environ["B03_NULL_DATA"] = "xenium_lung/crop/data"
os.environ["B03_NULL_OUT"] = "xenium_lung/crop/data/b10k"
os.environ["B03_NULL_CAP"] = "1"
os.environ["B03_NULL_B"] = "10000"
os.environ["B03_NULL_WORKERS"] = "4"
os.environ["B03_NULL_CKPT"] = "50"

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import pandas as pd
import b03_lung_config as cfg
import b03_robust_null_v3 as v3


def main():
    # label-count verification file (v3 asserts exact match)
    labels = np.load(f"{cfg.D}/labels.npy", allow_pickle=True)
    vc = pd.Series(labels).value_counts().to_dict()
    json_path = f"{cfg.D}/label_counts.json"
    import json
    json.dump({k: int(v) for k, v in vc.items()}, open(json_path, "w"))
    # parent-side configuration override (workers need only env-driven globals)
    v3.ALL_PAIRS = cfg.PAIRS
    v3.TYPE_MAP = cfg.TYPE_MAP
    v3.RECV_MAP = cfg.RECV_MAP
    v3.main()


if __name__ == "__main__":
    main()
