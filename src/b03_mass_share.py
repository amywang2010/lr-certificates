"""Tested-gene share of transcript mass on the two panel-density sections (paper Section 5.4).

The like-for-like comparison in Section 5.4 reads the same tissue and the same ten pairs
on a 313-plex section and on a 5,000-plex section, and attributes the collapse in
certification to the tested genes' share of transcript mass. This script computes that
share directly from the vendor molecule files, before any assignment is made:

    share = (QV >= 20 transcripts whose feature name is one of the tested genes)
            / (QV >= 20 transcripts in the section)

The tested genes are the union of the ligand and receptor genes of the ten pairs that the
200-cell class floor leaves evaluable on the second section, imported from the shared
configuration so the two legs cannot drift apart.

Inputs (vendor molecule files, not redistributed; see the Data availability section of the
README):
  - data/tx.parquet, the Rep1 molecule table written by b03_load.py (already QV >= 20)
  - <workspace>/xenium_breast_s6/crop/transcripts.parquet, the frozen second-section crop
    written by b03_s6_prep.py

Output: results/tested_mass_share.json
"""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_lung_config import PAIRS  # noqa: E402  (shared, frozen pair family)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
OUT = os.path.join(ROOT, "results")

EVALUABLE_PAIRS = ["CDH1->EGFR", "CDH1->ERBB2", "ERBB2->EGFR", "ERBB2->PDCD1", "ESR1->PGR",
                   "MS4A1->CD274", "PECAM1->CTLA4", "PECAM1->KDR", "PECAM1->PDCD1",
                   "PGR->ESR1"]
TESTS = [
    ("rep1", "Xenium v1.0.1 FFPE Human Breast Cancer Rep1",
     os.path.join(ROOT, "data", "tx.parquet")),
    ("s6", "second breast section, Xenium Prime 5K panel",
     os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "transcripts.parquet")),
]


def tested_genes(pairs):
    """Genes of a pair family, given as (ligand, receptor) tuples or "L->R" strings."""
    genes = set()
    for p in pairs:
        genes.update(p.split("->") if isinstance(p, str) else p)
    return sorted(genes)


def share_of(path, genes):
    """Fraction in percent of QV >= 20 transcripts carrying a tested gene, and its counts."""
    tx = pd.read_parquet(path)
    col = "feature_name" if "feature_name" in tx.columns else "gene"
    qcol = next((c for c in ("qv", "quality_score") if c in tx.columns), None)
    sub = tx if qcol is None else tx[tx[qcol] >= 20]
    mask = sub[col].isin(genes).to_numpy()
    return 100.0 * mask.mean(), int(len(sub)), int(mask.sum())


def main():
    genes = tested_genes(EVALUABLE_PAIRS)
    assert len(genes) == 11, f"expected 11 tested genes, found {len(genes)}"
    assert set(genes) <= set(tested_genes(PAIRS)), \
        "the evaluable pairs are not a subfamily of the frozen pair list"

    out = {
        "design": "tested-gene share of transcript mass on the two panel-density sections",
        "definition": "fraction of QV >= 20 transcripts in the section whose feature name is "
                      "one of the genes in the ten evaluable pairs, computed from the vendor "
                      "molecule file before any assignment is made",
        "genes": genes,
        "n_genes": len(genes),
        "pairs": EVALUABLE_PAIRS,
    }
    for key, section, path in TESTS:
        pct, n_all, n_in = share_of(path, genes)
        out[key] = {"section": section, "source": os.path.basename(path),
                    "n_transcripts_qv20": n_all, "n_in_tested_genes": n_in,
                    "share_pct": round(pct, 4)}
        print(f"{key}: {pct:.4f}% of {n_all:,} QV>=20 transcripts ({n_in:,} in tested genes)",
              flush=True)
        assert 0.0 < out[key]["share_pct"] < 100.0

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "tested_mass_share.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, sort_keys=True)
        f.write("\n")
    print("wrote results/tested_mass_share.json", flush=True)


if __name__ == "__main__":
    main()
