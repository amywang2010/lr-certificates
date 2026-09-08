"""B03 Phase 5 (step A): dataset arrival audit — per B03_PHASE5_PREREG.md gates G1/G2.

Runs the moment the zip lands. Zero statistics, zero assumptions:
  1. Zip CRC integrity on all members.
  2. G1 panel audit: which of the 25 preregistered pairs have BOTH genes in
     gene_panel.json -> the mechanical testable-pair subset (every exclusion
     printed with its missing gene).
  3. G2 schema audit: file formats, parquet columns, transcript header — dumped
     as facts.

Importable (main(zip_path, out_path)); executable directly with defaults.
Exit code 1 on integrity failure (zip corrupt, needed member missing).
"""
import gzip
import io
import json
import sys
import zipfile

sys.path.insert(0, "B03_project/src")
from b03_scout import ALL_PAIRS

ZIP_DEFAULT = "xenium_lung/Xenium_Prime_Human_Lung_Cancer_FFPE_outs.zip"
OUT_DEFAULT = "B03_project/logs/lung_audit.json"

NEED = [
    "gene_panel.json",
    "metrics_summary.csv",
    "cells.parquet",
    "cell_feature_matrix.h5",
    "cell_boundaries.parquet",
    "nucleus_boundaries.parquet",
]
TX_MEMBERS = ["transcripts.csv.gz", "transcripts.parquet"]  # at least one must exist


def main(zip_path=ZIP_DEFAULT, out_path=OUT_DEFAULT):
    tested_genes = sorted({g for p in ALL_PAIRS for g in p})
    zf = zipfile.ZipFile(zip_path)
    names = zf.namelist()
    print(f"zip members: {len(names)}")

    def find(member_candidates):
        return {c: next((m for m in names if m.split("/")[-1] == c), None)
                for c in member_candidates}

    needed_hits = find(NEED + TX_MEMBERS)
    missing = [c for c in NEED if not needed_hits[c]]
    if not any(needed_hits[t] for t in TX_MEMBERS):
        missing.append(f"transcripts.* (neither {TX_MEMBERS})")
    print("member map:", json.dumps(needed_hits, indent=1))
    if missing:
        print(f"MISSING NEEDED MEMBERS: {missing}")
        sys.exit(1)

    bad = zf.testzip()
    if bad is not None:
        print(f"ZIP CRC FAILURE at first bad member: {bad}")
        sys.exit(1)
    print("zip CRC integrity: all members OK")

    # ---- G1 panel audit ----
    panel = json.loads(zf.read(needed_hits["gene_panel.json"]))

    def collect_genes(obj, out):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k in ("gene", "gene_name", "symbol", "name") and isinstance(v, str):
                    out.add(v)
                else:
                    collect_genes(v, out)
        elif isinstance(obj, list):
            for it in obj:
                collect_genes(it, out)

    raw = set()
    collect_genes(panel, raw)
    panel_genes = {str(g).upper() for g in raw if isinstance(g, str)}
    print(f"panel genes found: {len(panel_genes)} (top-level keys: {list(panel.keys())[:6]})")

    present = [g for g in tested_genes if g in panel_genes]
    absent = [g for g in tested_genes if g not in panel_genes]
    testable_pairs = [f"{L}->{R}" for L, R in ALL_PAIRS if L in present and R in present]
    excluded_pairs = [(f"{L}->{R}", [g for g in (L, R) if g not in present])
                      for L, R in ALL_PAIRS if not (L in present and R in present)]
    print(f"G1: tested genes present {len(present)}/{len(tested_genes)}; "
          f"testable pairs {len(testable_pairs)}/25")
    if absent:
        print("  absent genes:", absent)
    for p, miss in excluded_pairs:
        print(f"  excluded pair {p}: missing {miss}")

    # ---- G2 schema audit ----
    import pandas as pd

    audit = dict(n_members=len(names), needed=needed_hits,
                 panel_n_genes=len(panel_genes),
                 tested_present=present, tested_absent=absent,
                 testable_pairs=testable_pairs, excluded_pairs=excluded_pairs,
                 panel_schema_keys=list(panel.keys()))

    with zf.open(needed_hits["cells.parquet"]) as fh:
        cells_df = pd.read_parquet(io.BytesIO(fh.read()))
    audit["cells_columns"] = list(cells_df.columns)
    audit["n_cells"] = len(cells_df)
    print(f"cells.parquet: {len(cells_df)} rows, cols {list(cells_df.columns)}")

    if needed_hits.get("transcripts.parquet"):
        with zf.open(needed_hits["transcripts.parquet"]) as fh:
            t0 = pd.read_parquet(io.BytesIO(fh.read()))
        audit["transcripts_columns"] = list(t0.columns)
        audit["n_transcripts"] = len(t0)
        print(f"transcripts.parquet: {len(t0)} rows, cols {list(t0.columns)}")
    else:
        with zf.open(needed_hits["transcripts.csv.gz"]) as fh:
            with gzip.open(fh, "rt") as g:
                head = [next(g) for _ in range(2)]
        audit["transcripts_csv_header"] = head[0].strip()
        audit["transcripts_csv_sample"] = head[1].strip()
        print("transcripts.csv.gz header:", head[0].strip())

    with open(out_path, "w") as f:
        json.dump(audit, f, indent=1, default=str)
    print(f"AUDIT COMPLETE -> {out_path}")
    return audit


if __name__ == "__main__":
    main()
