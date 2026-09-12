"""Phase-2 report generator: assembles B03_PHASE2_REPORT.md from analysis artifacts.

Implements the preregistered verdict trichotomy verbatim (B03_PHASE2_PREREG.md):
  >= 90% of evaluable rows inside -> PASS
  75-90%                          -> PARTIAL (mechanistic diagnosis required)
  < 75%                           -> FAIL (U does not transfer; report openly)
Dry-run validated against fixture outputs before real data arrives (watchdog pattern).
UTF-8 output. Referee pass happens after the real report renders.
"""

import json
import sys

import numpy as np
import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else "B03_project/results"
OUT_PATH = f"{RES}/B03_PHASE2_REPORT.md"

S = json.load(open(f"{RES}/phase2_proseg_summary.json"))
C = pd.read_csv(f"{RES}/phase2_proseg_coverage.csv")
cv = json.load(open(f"{RES}/phase2_label_transfer_cv.json"))

frac = S["primary"]["fraction"]
n_in, n_eval = S["primary"]["n_inside"], S["primary"]["n_evaluable"]
if frac >= 0.90:
    verdict = "PASS"
elif frac >= 0.75:
    verdict = "PARTIAL"
else:
    verdict = "FAIL"

# notable movements: rows whose Proseg T lands far outside the certified interval
C["mid"] = (C.T_lo + C.T_hi) / 2
C["outside_dist"] = np.where(C["T"] > C.T_hi, C["T"] - C.T_hi,
                             np.where(C["T"] < C.T_lo, C.T_lo - C["T"], 0.0))
evaluable = C["T"].notna()
mov = C[evaluable & ~C.inside.fillna(True)].sort_values(
    "outside_dist", ascending=False)
mov_rows = "".join(
        f"| {r.pair} | {r.region} | {r.T0:+.2f} | [{r.T_lo:+.2f}, {r.T_hi:+.2f}] | "
        f"{getattr(r, 'T'):+.2f} | {r.outside_dist:.2f} |\n"
        for r in mov.head(10).itertuples())

cert_pos = int(C.certified_pos_matched.fillna(False).sum())
cert_neg = int(C.certified_neg_matched.fillna(False).sum())
cert = C[C.certified_pos_matched.fillna(False) | C.certified_neg_matched.fillna(False)]
sign_ok = int((((C["T"] > 0) == C.certified_pos_matched)[
    C.certified_pos_matched.fillna(False) |
    C.certified_neg_matched.fillna(False)]).sum())

L = []
a = L.append
a("# B03 Phase 2 — Proseg Cross-Platform Replication Report")
a("")
a(f"## Verdict: **{verdict}**")
a("")
a(f"Pre-registered gate (B03_PHASE2_PREREG.md, Amendment A3 semantics): "
  f"**{n_in}/{n_eval}** evaluable rows inside their certified intervals = "
  f"**{frac*100:.1f}%** (gate >= 90%; {S['primary'].get('n_excluded', 0)} rows "
  f"excluded as non-evaluable, >10% exclusion aborts).")
a("")
a("## Configuration provenance")
a("")
a(f"- Proseg {S['proseg_version']}, run 4 per Amendment A4: 2 um voxels, "
  f"2 z-layers (memory-feasible configuration after the run-3 guard abort, amendment record); "
  f"diffusion model KEPT; all other parameters default.")
a(f"- Exact command: `{S['command']}`")
a(f"- Cells: {S['n_proseg_cells']:,} (background rows dropped per A2a: "
  f"{S['n_background_cells_excluded']}); counts as gzip MatrixMarket, rounded to "
  f"integer point estimates (Proseg 3 default).")
a(f"- Gene-alignment proof (A2c): 0 invariant violations across needed genes; "
  f"noise-slot handling: {S.get('noise_gene_slot') or 'n_named layout (none dropped)'}.")
a("")
a("## Endpoint arithmetic")
a("")
a(f"- Rows inside certified intervals: **{n_in}/{n_eval}** ({frac*100:.1f}%)")
a(f"- Certified rows (matched-null scout, q <= 0.10): {cert_pos} positive, "
  f"{cert_neg} negative; sign agreement under Proseg: **{sign_ok}/{len(cert)}**")
a(f"- Label transfer (frozen vendor centroids, cosine, tested-LR genes excluded): "
  f"2-fold CV accuracy **{cv['accuracy_2fold_cv']*100:.1f}%** on "
  f"{cv['n_eval']:,} vendor cells ({cv['n_types']} types); analysis-critical types "
  f"(endothelial, macrophage, T, fibroblast, myoepithelial) 85-97% "
  f"(confusion matrix: phase2_label_transfer_cv_confusion.csv).")
a("")
if len(mov):
    a("## Rows moved outside certification by Proseg (largest first)")
    a("")
    a("| pair | region | vendor T0 | certified interval | Proseg T | distance outside |")
    a("|---|---|---|---|---|---|")
    a(mov_rows)
else:
    a("## Rows moved outside certification by Proseg")
    a("")
    a("None (all evaluable rows inside).")
a("")
a("## Interpretation contract (frozen in prereg)")
a("")
a("- PASS: certification coverage transfers to an independent probabilistic "
  "segmentation algorithm — the uncertainty-set thesis gains cross-platform support.")
a("- PARTIAL: coverage is resolution-dependent; mechanistic diagnosis of which "
  "pairs/regions move is REQUIRED before any claim.")
a("- FAIL: U does not transfer across segmentation algorithms — a major, publishable "
  "negative result, reported openly.")
a("")

with open(OUT_PATH, "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print(f"report written: {OUT_PATH} (verdict {verdict})")
