"""B03 Phase 6: breast B=10,000 matched-null driver (B03_PHASE6_PREREG.md, frozen 2026-09-08).

Configuration of the sealed Phase 1 run (seed 20260905, B03_project/data,
CAP unset: the amendment record cap never bound on breast and must stay off so the run
reproduces the sealed configuration), except B_PERM=10000 and a fresh OUT
directory. Per-perm seeds 20260905+100003+p: perms 0..999 are bitwise identical
to the B=1000 run; the chain script verifies this determinism check against
null_v3_checkpoint.npz before any output is used. A check failure on either
tissue kills Phase 6 and triggers a code audit (registered anomaly path).
"""
import os
import sys

os.environ["B03_NULL_SEED"] = "20260905"
os.environ["B03_NULL_DATA"] = "B03_project/data"
os.environ["B03_NULL_OUT"] = "B03_project/results/b10k_breast"
os.environ["B03_NULL_B"] = "10000"
os.environ["B03_NULL_WORKERS"] = "4"
os.environ["B03_NULL_CKPT"] = "50"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_robust_null_v3 as v3

if __name__ == "__main__":
    v3.main()
