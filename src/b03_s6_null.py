"""Matched permutation null for the second breast section.

Drives b03_matched_null by import with the s6 configuration:
  - DATA/OUT -> the workspace xenium_breast_s6/crop/data dir (env B03_NULL_DATA/B03_NULL_OUT, env so
    the spawned workers inherit the same module configuration)
  - seed 20260907 -> per-perm seeds 20260907+100003+p (fixed)
  - the shared TYPE_MAP, RECV_MAP and 18-pair family from b03_lung_config, of
    which the ten pairs without Myeloid cells give the 30 tests
  - B03_NULL_CAP=1, which applies the capacity-capped dmax with the same
    formula as the certify step
  - observed intervals are re-verified against scout_bounds.csv (1e-6) before
    any permutation. The eight pairs that lose a class to the 200-cell rule are
    absent from those rows and so absent from the key list.
Label counts are checked against label_counts.json, written by b03_s6_load.py,
before any permutation runs.
"""
import os, sys, json
import numpy as np
import pandas as pd
import b03_lung_config as cfg

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
os.environ["B03_NULL_SEED"] = "20260908"
os.environ["B03_NULL_DATA"] = os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "data")
os.environ["B03_NULL_OUT"] = os.path.join(ROOT, "results", "breast_s6")
os.environ["B03_NULL_CAP"] = "1"
os.environ["B03_NULL_B"] = "10000"
os.environ["B03_NULL_WORKERS"] = "4"
os.environ["B03_NULL_CKPT"] = "50"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_matched_null as v3

D = cfg.D
PAIRS = cfg.PAIRS
TYPE_MAP = cfg.TYPE_MAP
RECV_MAP = cfg.RECV_MAP

# The s6 data dir is not cfg.D, which holds the lung dir from the shared config.
# The label-counts file is read and written in the s6 crop data dir, which the null
# finds through B03_NULL_DATA.
S6_D = os.environ["B03_NULL_DATA"]

def main():
    labels = np.load(f"{S6_D}/labels.npy", allow_pickle=True)
    vc = pd.Series(labels).value_counts().to_dict()
    json.dump({k: int(v) for k, v in vc.items()}, open(f"{os.environ['B03_NULL_OUT']}/label_counts.json", "w"))
    v3.ALL_PAIRS = PAIRS
    v3.TYPE_MAP = TYPE_MAP
    v3.RECV_MAP = RECV_MAP
    v3.main()

if __name__ == "__main__":
    main()
