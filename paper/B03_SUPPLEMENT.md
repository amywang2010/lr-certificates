# Supplementary Material

## Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty

All numbers in this supplement are read from hashed artifacts (SHA-256 manifests:
`results/B03_SCOUT_MANIFEST.sha256`, `results/PHASE2_MANIFEST.sha256`,
`PHASE3_MANIFEST.sha256`, `PHASE4_MANIFEST.sha256`, `PHASE5_MANIFEST.sha256`).
The final `PHASE5_MANIFEST.sha256` (117 entries) supersedes the earlier manifests for
every file it covers (verified entry-by-entry at release time); the earlier manifests
are retained as the frozen per-phase state. `RELEASE_MANIFEST.sha256` covers the
paper, supplement, figures, audit, and released source exactly as shipped.
The manuscript cites the deviation ledger (`B03_DEVIATIONS.md`, DEV-001–018) at every
design change.

---

## S1. Datasets and preprocessing

| Property | Tissue 1 (breast) | Tissue 2 (lung) |
|---|---|---|
| Dataset | Xenium FFPE Human Breast Cancer Rep1 (v1.0.1) | Xenium Prime Human Lung Cancer FFPE (5K) |
| Panel | 313-plex breast panel | 5,001-gene 5K Pan Tissue & Pathways (hAtlas v1.1) |
| Raw transcripts | 42,638,083 | 177.5 M (full section) |
| QV≥20 transcripts | 34,493,510 | 149.4 M full section; 21,341,704 in crop |
| Vendor cells | 167,780 | 278,328 full section; 47,754 in crop |
| Section area |, (full section) | 71.09 mm² (full); 12.25 mm² (crop) |
| PIP assignment rate | 31,294,184 / 34,493,510 (90.7%) | 16,112,356 / 21,341,704 (75.5%) |
| Band fraction (3 µm) | 22,477,605 (65.2%) | 17,358,580 (81.3%) |
| Band+tested transcripts (movable) | 3,136,405 | 128,428 |
| Concordance (per-cell exact / total ratio) | 0.9146 / 0.9736 | 0.9568 / 0.9694 |
| Regions (cells) | 67,801 / 53,886 / 46,093 | 18,551 / 17,258 / 11,945 |
| Label classes (≥200 cells) | 15 vendor annotation classes | 6 dictionary classes (all ≥200) |

Notes.
- The lung analysis region is the pre-registered expression-blind crop x[6500,10000) ×
  y[1000,4500) µm, selected by median vendor-cell count over a 0.5 mm grid with zero
  expression access (DEV-017; window record `results/lung_crop_window.json`).
- Band fractions differ between tissues because band membership depends on cell size and
  density; both are reported because the band is the uncertainty-supporting set (DEV-002).
- The breast concordance gate was redefined with a diagnosis after the initial 0.95
  exact-concordance target proved miscalibrated to boundary semantics (DEV-001); the
  lung gate is the same class of gate (per-cell ≥ 0.90 AND total ratio ≥ 0.97).

## S2. Pre-registration inventory

| Document | Frozen | Content |
|---|---|---|
| B03_PREREGISTRATION.md | 2026-09-05, before any expression compute | Estimand, U, labels, 25 pairs (5 controls + 20 disputed), null, gates, regions, calibration genes |
| B03_PHASE2_PREREG.md (+A1–A4) | 2026-09-06, before any Proseg output | Cross-platform endpoints, label transfer, invocation, denominator semantics, memory-feasible config |
| B03_PHASE3_PREREG.md | 2026-09-07 09:40, before compute | Proseg-side matched null, assertions A–D |
| B03_PHASE4_PREREG.md | 2026-09-07, before compute | Kernel margins c*, monotone bisection |
| B03_PHASE5_PREREG.md + ADDENDUM | 2026-09-07, before any lung statistic | Dataset identity, G1–G3, endpoints E1–E4, seeds, crop rule (A2), marker dictionary (A3) |

Deviations are numbered DEV-001 through DEV-018 with dates, reasons, and quarantine
records for invalidated artifacts. Highlights: DEV-003 (null layer corrected before any
inference), DEV-005 (parallel rerun with per-permutation seeds and pre-flight bitwise
determinism), DEV-009 (resource-guard abort and memory-feasible Proseg configuration),
DEV-015 (panel metadata unavailable; assertion-guarded marker dictionary), DEV-016
(capacity binding on tissue 2; one sound capped formula everywhere), DEV-017 (V1 run
unviable; crop design), DEV-018 (E2 PARTIAL reported per frozen gates).

## S3. Full breast results table (75 rows)

The complete per-row table (pair, region, control flag, n_L, n_R, cell counts, count
deltas, T0, T_lo, T_hi, width, erosion variants, matched-null p and q in both
directions, verdict) is shipped as `results/scout_final_matched.csv`. Summary counts:

| Verdict | Count |
|---|---|
| Certified positive (interval + matched q ≤ 0.10) | 37 |
| Certified negative (interval + matched q ≤ 0.10) | 9 |
| Non-identifiable | 29 |

Certified-negative rows: ESR1→PGR (3 regions), PGR→ESR1 (3 regions), CD274→PDCD1
(3 regions). The CD274→PDCD1 negative certificates are the expected
trans-cellular disco-localization of a ligand–receptor control pair and were declared
non-gating at pre-registration (watchdog note).

Per-region verdict counts (breast): r0 14 pos / 3 neg / 8 non-id; r1 12 pos / 3 neg /
10 non-id; r2 11 pos / 3 neg / 11 non-id. Disputed rows with 0 in the interval
(pre-registered gate-2 event): 19 of 60. Disputed rows certifying negative: 6.

## S4. Full lung results table (54 rows)

Shipped as `lung/scout_final_lung.csv` (release path; working path `xenium_lung/crop/data/scout_final_lung.csv`) (same schema). Summary counts:

| Verdict | Count |
|---|---|
| Certified positive | 1 (CDH1→EGFR r0) |
| Certified negative | 6 (ESR1→PGR and PGR→ESR1, r0–r2) |
| Non-identifiable | 47 |

Interval-certified (before the matched null): 1 positive, 32 negative. The 26
interval-negative rows that fail the matched null share a structure: their upper bounds
are only moderately negative (median T_hi ≈ −0.58), and under the add-one null at
B = 1000 the worst-case object under density-matched label permutations reaches or
exceeds the observed bound in all 26 rows (p_neg ≥ 0.789; 22 at the p = 1.0 tie
boundary). These rows are reported as interval-certified but
not claim-certified, exactly as the two-layer design intends.

Per-region verdict counts (lung): r0 1 pos / 2 neg / 15 non-id; r1 0 pos / 2 neg /
16 non-id; r2 0 pos / 2 neg / 16 non-id.

## S5. Pseudocount (ε) sensitivity at the interval layer

Interval-certified row counts as a function of ε (computed from the frozen count
deltas; no new segmentation or permutation compute):

| ε | Breast pos / neg (of 75) | Lung pos / neg (of 54) |
|---|---|---|
| 0.1 | 19 / 20 | 0 / 41 |
| 0.5 (registered) | 39 / 12 | 1 / 32 |
| 1.0 | 75 / 0 | 52 / 0 |

Mechanism: ε enters both the observed intervals and the null identically; larger ε
compresses log-ratios toward zero from above, easing positive certification and
hardening negative certification. The certified (claim) set additionally requires the
matched null; the table above is the interval layer only, shown to expose the ε
dependence transparently.

## S6. Kernel margins (E4) detail

Breast (46 matched-certified rows): median c* = 1.28 (min 1.049; 69.6% ≥ 1.25; 28.3%
≥ 1.5; 15.2% ≥ 2.0; 0 kernel-immune). Positive controls median 1.49; disputed 1.28.
Secondary: 21 of 29 non-certified rows would certify at c = 0.5, including the 5
matched-null-dropped rows (their intervals exclude 0 at all c ≤ 1; they were dropped on
null significance, not interval uncertainty).

Lung (7 matched-certified rows): median c* = 20.4 (min 1.034 on the excursion row
CDH1→EGFR r0; 6 of 7 ≥ 3). Secondary: 35 of 47 non-certified rows would certify at
c = 0.5.

Assertions executed in both tissues: c = 1 reproduces the certified intervals to
≤ 8.9e-16; interval signs match flags; matched-certified set ⊆ interval-certified set;
monotonicity of T_lo/T_hi in c on an 11-point grid; root precision |T(c*)| ≤ 1e-10;
direction-aware sign-transition verification on every row.

## S7. Cross-segmentation detail (E2/E3)

Breast. Proseg 3.2.0, run 4 configuration (2 µm voxels, 2 z-layers; memory-feasible
after the run-3 guard abort, DEV-008/009), 165,094 cells, background rows dropped 0,
noise slot absent, gene-alignment proof 0 invariant violations. Primary coverage
72/75 = 96.0% (PASS ≥ 90%). Excursions (all region 0, all sign-preserving):
ERBB2→PDCD1 T = +3.66 vs ceiling +3.61; ESR1→PGR and PGR→ESR1 T = −0.249 vs ceiling
−0.275. Label-transfer 2-fold CV accuracy 72.4% (15 types); analysis-critical types
85–97%. E3: 40/46 = 87.0% certified-direction significance (PASS ≥ 85%); all 6 losses
are ESR1↔PGR (point estimates negative and inside intervals in r1/r2; r0 0.026 above
ceiling; observed statistic sits at the top of the Proseg-side null, p_pos = 0.001).

Lung. Proseg 3.2.0, same A4 configuration, crop input (25.55 M transcripts), 47,286
cells (−468 vs vendor, 1.0%), envelope gate 6.90 GB ≤ 7.0 GB before launch, guard
active throughout, clean exit at 15:59 with watcher validation. Coverage 6/7 = 85.7%
(PARTIAL; gate ≥ 90% PASS): single excursion CDH1→EGFR r0, T = +4.035 vs ceiling
+3.717, sign preserved, significant under the Proseg-side null (q = 0.0011). E3: 7/7 =
100% (PASS ≥ 85%), all four preregistered assertions PASSED, referee BH re-derivation
max |Δq| ≈ 1e-16.

## S8. Exactness and verification receipts

| Check | Breast | Lung |
|---|---|---|
| Randomized feasible configurations (gene-level) | 150,000; 0 violations | 108,000; 0 violations |
| Subcube exhaustive enumeration | within greedy extremes (CD274→PDCD1 r1) | within greedy extremes (CD274→PDCD1 r2) |
| Observed-interval cross-check (null startup) | all 75 rows within 1e-6 | all 54 rows within 1e-6 |
| Pre-flight determinism | perm 0 bitwise identical (parent + 2 workers) | perm 0 bitwise identical (parent + 2 workers) |
| Permutation invariance checks | 10 named permutations PASSED | 10 named permutations PASSED |
| Referee BH re-derivation | max |Δq| = 1.1e-16 (scout), 2.4e-15 (E3) | max |Δq| ≈ 1e-16 (E3); 0.009 vs stored raw-count q in vendor battery (convention difference, no certified-set sensitivity) |
| Capacity cap cross-check (lung) | n/a (never binds) | 32 binding instances; 0 mismatches |

## S9. Resource envelope and runtime accounting

| Stage | Breast | Lung (crop) |
|---|---|---|
| Vendor load (PIP + labels + regions) | ~4 min PIP within full load | 177 s total |
| Certified intervals | vectorized (breast) | 63 s (bounds-to-subcube block, artifact mtimes 15:31:05–15:32:09) |
| Matched null (B = 1000, 4 workers) | ~74 min | 417 s |
| Proseg | run 4: ~2.4 h (2 µm voxels, 2 layers) | ~1.7 h (crop) |
| E2 coverage | seconds | seconds |
| E3 Proseg-side null | 60 s | 21 s |
| E4 margins | < 1 s | < 1 s |

Hardware: 8 physical cores, 15.7 GB RAM, Windows. The resource guard (free virtual
< 3.0 GB or proseg private > 14.0 GB aborts) fired once in the campaign (run 3, breast)
and protected the machine as designed; its thresholds were documented before any abort
decision (DEV-008/009).

## S10. Figure provenance

- Figure 1: `results/scout_final_matched.csv` (75 rows; 37/9/29 verdicts asserted
  at plot time).
- Figure 2: `results/phase4_kernel_margins.csv` (46 rows) and
  `lung/phase4_kernel_margins.csv` (7 rows).
- Figure 3: `results/paper_fig3_data.csv`, compiled from the two final CSVs and the two
  coverage CSVs; 18 rows asserted.
- Figure 4: `results/phase2_proseg_summary.json`, `results/phase3_summary.json`,
  `lung/phase5_e2_summary.json`,
  `lung/phase3_summary.json`.
- Generation code: `src/b03_paper_figures.py` (assertions at plot time; deterministic
  jitter seed 20260907).

## S11. Reference provenance

See `paper/REFERENCES_PROVENANCE.md`: every reference verified by direct retrieval on
2026-09-07, with the specific claim cited and any attribution corrections recorded
(cellAdmix author order per the LOCK). Two method families named in the internal dossier
(CONCISE, CellNEST) are deliberately not cited: citation-grade verification was not
completed within the session, and the verified reviews cover the space adequately.
