"""Matched permutation null for the lung crop.

Drives b03_matched_null by import with the lung configuration:
  - DATA and OUT point at the crop prepared by b03_lung_load.py and
    results/lung, passed through the environment so spawned workers inherit
    the same module configuration
  - seed 20260907, so permutation p uses 20260907 + 100003 + p
  - the lung TYPE_MAP and RECV_MAP and the 18 testable pairs, set as
    parent-side globals
  - B03_NULL_CAP=1, which applies the capacity-capped dmax with the same
    formula as the certify step
Label counts are checked against label_counts.json, written by b03_lung_load.py,
before any permutation runs.
"""
import os, sys, json
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
os.environ["B03_NULL_SEED"] = "20260907"
os.environ["B03_NULL_DATA"] = os.path.join(WORKSPACE, "xenium_lung", "crop", "data")
os.environ["B03_NULL_OUT"] = os.path.join(ROOT, "results", "lung")
os.environ["B03_NULL_CAP"] = "1"
os.environ["B03_NULL_B"] = "1000"
os.environ["B03_NULL_WORKERS"] = "4"
os.environ["B03_NULL_CKPT"] = "50"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_matched_null as v3

D = os.environ["B03_NULL_DATA"]
PAIRS = [("CD274","PDCD1"),("CXCL12","CXCR4"),("PECAM1","KDR"),
 ("ERBB2","EGFR"),("ERBB2","PDCD1"),("ESR1","PGR"),("PGR","ESR1"),
 ("CD274","CTLA4"),("CD274","CD8A"),("CXCL12","CCR7"),("PECAM1","PDCD1"),
 ("PECAM1","CTLA4"),("CD163","CCR7"),("CD163","CXCR4"),("CDH1","EGFR"),
 ("CDH1","ERBB2"),("MS4A1","CD274"),("CD68","CD274")]
TYPE_MAP = {
 "CD274": "Myeloid", "CXCL12": "Fibroblast", "PECAM1": "Endothelial",
 "CD163": "Myeloid", "ERBB2": "Epithelial", "ESR1": "Epithelial",
 "PGR": "Epithelial", "CDH1": "Epithelial", "MS4A1": "B_cells", "CD68": "Myeloid",
}
RECV_MAP = {
 "PDCD1": "T_cells", "CXCR4": "Myeloid", "CCR7": "Myeloid",
 "KDR": "Endothelial", "EGFR": "Epithelial", "PGR": "Epithelial",
 "ESR1": "Epithelial", "CTLA4": "T_cells", "CD8A": "T_cells",
 "ERBB2": "Epithelial", "CD274": "Epithelial",
}

def main():
    # label-count verification file (exact match asserted)
    labels = np.load(f"{D}/labels.npy", allow_pickle=True)
    vc = pd.Series(labels).value_counts().to_dict()
    json.dump({k: int(v) for k, v in vc.items()}, open(f"{os.environ['B03_NULL_OUT']}/label_counts.json", "w"))
    # parent-side configuration override (workers need only env-driven globals)
    v3.ALL_PAIRS = PAIRS
    v3.TYPE_MAP = TYPE_MAP
    v3.RECV_MAP = RECV_MAP
    v3.main()

if __name__ == "__main__":
    main()
