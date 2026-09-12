"""B03 Phase 5 (step C): E3 proseg-side matched null — phase3 machinery, re-pointed.

Drives b03_phase3_proseg_null by import with the tissue-2 configuration:
env-driven DATA/RES/rows; b03_scout globals patched BEFORE phase3 import so its
module-level `from b03_scout import ALL_PAIRS, TYPE_MAP, RECV_MAP` binds the lung
values. Seeds: 20260907+100003+p (prereg item 5; phase3's SEED is already 20260907).
Gate E3: >= 85% of certified rows significant in the certified direction under the
proseg-side null.
"""
import os, sys

os.environ["B03_P3_DATA"] = "xenium_lung/crop/data"
os.environ["B03_P3_RES"] = "xenium_lung/crop/data"
os.environ["B03_P3_ROWS"] = "54"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_scout

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
b03_scout.ALL_PAIRS = PAIRS
b03_scout.TYPE_MAP = TYPE_MAP
b03_scout.RECV_MAP = RECV_MAP

import b03_phase3_proseg_null as p3

if __name__ == "__main__":
    p3.main()
