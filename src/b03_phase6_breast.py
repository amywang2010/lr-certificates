"""B03 Phase 6: breast B=10,000 matched-null driver (design frozen 2026-09-08).

Configuration of the shipped Phase 1 run (seed 20260905, repository data/ directory,
CAP unset: the capacity cap never bound on breast and must stay off so the run
reproduces the shipped configuration), except B_PERM=10000 and a fresh OUT
directory. Per-perm seeds 20260905+100003+p: perms 0..999 are bitwise identical
to the B=1000 run; the chain script verifies this determinism check against
null_v3_checkpoint.npz before any output is used. A check failure on either
tissue kills Phase 6 and triggers a code audit (registered anomaly path).
"""
import os
import sys

os.environ["B03_NULL_SEED"] = "20260905"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ["B03_NULL_DATA"] = os.path.join(ROOT, "data")
os.environ["B03_NULL_OUT"] = os.path.join(ROOT, "results", "b10k_breast")
os.environ["B03_NULL_B"] = "10000"
os.environ["B03_NULL_WORKERS"] = "4"
os.environ["B03_NULL_CKPT"] = "50"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_robust_null_v3 as v3

if __name__ == "__main__":
    v3.main()
