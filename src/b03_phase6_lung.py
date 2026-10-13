"""B03 Phase 6: lung B=10,000 matched-null driver (the high-B configuration).

Mirrors the B=1,000 driver (b03_lung_null.py): same seed (20260907), same data
dir, same pair list and type maps from the shared config module, same
B03_NULL_CAP=1. Differences: B=10000 and a fresh OUT directory, so the B=1,000
artifacts are left untouched. Per-perm seeds 20260907+100003+p: perms 0..999 are
bitwise identical to the B=1,000 run, and the chain script enforces the
exact-match determinism gate against the B=1,000 checkpoints before any output
is used.
"""
import os
import sys

os.environ["B03_NULL_SEED"] = "20260907"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
os.environ["B03_NULL_DATA"] = os.path.join(WORKSPACE, "xenium_lung", "crop", "data")
os.environ["B03_NULL_OUT"] = os.path.join(WORKSPACE, "xenium_lung", "crop", "data", "b10k")
os.environ["B03_NULL_CAP"] = "1"
os.environ["B03_NULL_B"] = "10000"
os.environ["B03_NULL_WORKERS"] = "4"
os.environ["B03_NULL_CKPT"] = "50"

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import pandas as pd
import b03_lung_config as cfg
import b03_matched_null as v3


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
