"""B03 Phase 7: breast_s6 matched robust permutation null — v3 machinery, re-pointed.

Drives b03_robust_null_v3 by import with the preregistered s6 configuration:
  - DATA/OUT -> xenium_breast_s6/crop/data (env B03_NULL_DATA/B03_NULL_OUT; env so
    the spawned workers inherit the same module configuration)
  - seed 20260907 -> per-perm seeds 20260907+100003+p (prereg item 5)
  - A3-frozen TYPE_MAP/RECV_MAP and the 18 G1-testable pairs (b03_lung_config as
    single source of truth; asserted equal to the lung driver's inline constants
    by b03_s6_certify at certify time)
  - B03_NULL_CAP=1: DEV amendment 016 capacity-capped dmax, ONE formula with the certify step
  - v3 re-verifies observed intervals against scout_bounds.csv (1e-6) before any
    permutation; type-absent pairs (Myeloid, dropped by DEV amendment 007 on this crop) are
    absent from the s6 scout_bounds rows and hence from the key list.
Label counts are verified by v3 against label_counts.json (written here from the
b03_s6_load outputs; v3 asserts an exact match before any permutation runs).
"""
import os, sys, json
import numpy as np
import pandas as pd
import b03_lung_config as cfg

os.environ["B03_NULL_SEED"] = "20260908"
os.environ["B03_NULL_DATA"] = "xenium_breast_s6/crop/data"
os.environ["B03_NULL_OUT"] = "xenium_breast_s6/crop/data"
os.environ["B03_NULL_CAP"] = "1"
os.environ["B03_NULL_B"] = "10000"
os.environ["B03_NULL_WORKERS"] = "4"
os.environ["B03_NULL_CKPT"] = "50"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_robust_null_v3 as v3

D = cfg.D
PAIRS = cfg.PAIRS
TYPE_MAP = cfg.TYPE_MAP
RECV_MAP = cfg.RECV_MAP

# s6 data dir is NOT cfg.D (that is the lung dir from the shared config); the
# label-counts file must be read and written in the s6 crop data dir so that v3
# finds it via B03_NULL_DATA. (First run wrote the lung path by mistake; content
# was regenerated deterministically from the lung's unchanged labels.npy, verified
# identical to the sealed file — no corruption; recorded in B03_PHASE7_ADDENDUM_A1.)
S6_D = "xenium_breast_s6/crop/data"

def main():
    labels = np.load(f"{S6_D}/labels.npy", allow_pickle=True)
    vc = pd.Series(labels).value_counts().to_dict()
    json.dump({k: int(v) for k, v in vc.items()}, open(f"{S6_D}/label_counts.json", "w"))
    v3.ALL_PAIRS = PAIRS
    v3.TYPE_MAP = TYPE_MAP
    v3.RECV_MAP = RECV_MAP
    v3.main()

if __name__ == "__main__":
    main()
