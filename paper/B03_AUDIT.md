# B03 Strict-Reviewer Audit (pre-release, 2026-09-07)

Scope: every line of `B03_MANUSCRIPT.md`, `B03_SUPPLEMENT.md`, `REFERENCES_PROVENANCE.md`,
the four figures and their generation code, re-read against the hashed artifacts without
reliance on session memory. Each item below states the check, the evidence, and the
disposition.

## 1. Number-by-number verification (manuscript vs artifacts)

| Claim in text | Artifact | Result |
|---|---|---|
| Breast 75 rows; 37 pos / 9 neg / 29 NI | `results/scout_final_matched.csv` | exact |
| Disputed rows 60; 19 collapse; 6 certify negative; max disputed point estimate +2.24 | same (non-control subset) | exact |
| Controls: 7 of 15 control rows certify positive | same | exact |
| ERBB2→EGFR r0: T0 +3.71, interval [+1.55, +4.58], q_pos = 1.0 | same | exact |
| ESR1→PGR r0: T0 −0.71, T_hi −0.27, width 1.23 | same | exact |
| CD68→CD274 r0: T_lo +0.07 | same | exact |
| Certified-positive width median 3.34, IQR 2.40–4.58 | same | **fixed** (was 3.21 / 2.23–3.74) |
| Interval layer: 39 pos / 12 neg (ε = 0.5) | same | exact |
| Breast E2: 72/75 = 96.0%; 3 excursions r0, 0.03–0.05, sign-preserving | `results/phase2_proseg_summary.json`, `phase2_proseg_coverage.csv` | exact (excursions 0.052/0.025/0.025) |
| Breast E3: 40/46 = 87.0%; losses = ESR1↔PGR ×3 regions ×2 directions | `results/phase3_summary.json`, `phase3_proseg_null.csv` | exact |
| Breast E4: median c* 1.28; controls 1.49; 69.6% ≥ 1.25; 15.2% ≥ 2.0; 21 of 29 at c = 0.5 | `results/phase4_summary.json` | exact |
| ε-sensitivity 19/20 (0.1), 39/12 (0.5), 75/0 (1.0) breast; 0/41, 1/32, 52/0 lung | re-derived independently from count deltas | exact |
| d = 2 sensitivity: 75/75 breast, 54/54 lung inside | `T_erode2` columns, both final CSVs | exact |
| Nucleus outside-U: 7 breast rows, all CXCL12/CD68 senders | `results/coverage_check.csv` | exact |
| Lung 54 rows; 1 pos (CDH1→EGFR r0 [+0.21, +3.72], T0 +3.39); 6 neg (ESR1↔PGR, T_hi ∈ [−1.88, −1.79], q_neg ≤ 0.027) | `xenium_lung/crop/data/scout_final_lung.csv` | exact (**q bound fixed** from 0.036) |
| Lung certified width median 0.18, range 0.11–0.20 | same | exact |
| Lung E2: 6/7 = 85.7% PARTIAL; excursion CDH1→EGFR r0 +4.04 vs +3.72; q = 0.0011 | `phase5_e2_summary.json`, `phase5_e2_coverage.csv`, `phase3_proseg_null.csv` | exact |
| Lung E3: 7/7 = 100% | `phase3_summary.json` (crop) | exact |
| Lung E4: median c* 20.4; min 1.03 (excursion row); 6/7 ≥ 3; 35 of 47 at c = 0.5 | `phase4_summary.json` (crop), `phase4_kernel_margins.csv` | exact |
| Concordance: breast per-cell 0.9146 / ratio 0.9736 / per-gene r median 0.991 min 0.958; lung 0.9568 / 0.9694; lung full-section 20/20 bitwise | `data/concordance*.json`, `xenium_lung/crop/data/concordance.json`, `logs/lung_concordance.log` | exact |
| Band fractions 65.2% (22,477,605) breast; 81.3% (17,358,580) lung; movable 3,136,405 / 128,428 | `logs/load.log`, `logs/finish.log`, `xenium_lung/crop/lung_load.log` | exact |
| PIP rates 90.7% / 75.5% | same logs | exact |
| Region sizes 67,801/53,886/46,093 and 18,551/17,258/11,945 | scout report; `label_counts`/regions | exact |
| Capacity: 197 cells bind, 41 zero-count artifacts, 32 cross-check instances | `xenium_lung/crop/data/capacity_check.json` | exact |
| 15 label-transfer types (breast), CV 72.4%, critical types 85–97% | `phase2_label_transfer_cv.json`, `b03_phase2_vendor_prep.py` (≥50-cell rule) | exact; manuscript label-class wording **fixed** (20 panel classes → 15 transfer targets) |
| Lung section 149.4M QV≥20; 71.09 mm²; 177.5M all-decoded; crop 25.55M / 21.34M QV≥20; 47,754 cells | `xenium_lung/extracted/metrics_summary.csv`; `B03_PHASE5_*` | exact; basis-mismatch **fixed** (was "177.5M … 5.1×"; now quotes QV≥20 both tissues, 4.3×) |
| V1 run: 541 features, missing 5 tested genes | `B03_PHASE5_ADDENDUM.md` (audit) | exact |
| Proseg cells 165,094 (breast) / 47,286 (lung crop) | phase summaries, `proseg_done.json` | exact |
| Runtimes: null 74 min breast / 417 s lung; proseg ~2.4 h breast / ~1.7 h lung; E3 60 s / 21 s | `B03_DEVIATIONS.md`, guard logs, artifact mtimes | proseg-lung **fixed** from "~1 h" to ~1.7 h (guard log 14:16:50 → 15:59:12) |

## 2. Logic and overclaiming audit

- The phrase "certified cross-tissue sign flip" (internal phase reports) was examined and
  deliberately **not** used in the manuscript: both tissues certify the negative direction
  for ESR1↔PGR; the cross-tissue contrast is fragility vs robustness, not a sign change.
  Section 7.3 states this explicitly.
- Every "PASS/PARTIAL" verdict is tied to a gate number frozen in a prereg document;
  no post hoc gate softening exists in the text.
- Non-identifiable rows are described as non-identifiability, never as absence of
  interaction (no equivalence-claim fallacy).
- The abstract's claims are scoped to "boundary-local reassignment, mask erosion/dilation,
  per-gene leakage intervals"; no claim of coverage over merge/split or label uncertainty
  is made anywhere outside the limitations.
- The five matched-null drops are reported as the framework refusing to overclaim, with
  the mechanism (worst-case layer not extreme under matched permutations).

## 3. Reference audit

- All 20 listed references are cited in the text (Salas et al. was initially uncited;
  now cited in Section 7.2). Reference 19/20 are dataset citations.
- Provenance ledger records verification mode per reference; two unverifiable method
  families (CONCISE, CellNEST) are deliberately not cited, and the manuscript makes no
  claims that depend on them.
- The ESR1/PGR co-expression claim cites Perou et al. (2000); the add-one convention
  cites Phipson & Smyth (2010) and North et al. (2002); integrality cites Ahuja et al. (1993).

## 4. Style audit (anti-AI-slop)

- Zero em dashes in all three paper documents after the sweep; parenthetical and
  appositive constructions were recast where a bare comma splice would have read poorly.
- No signposting transitions ("Moreover", "Furthermore", "Additionally" as connectives),
  no "delve/landscape/leverage/underscore" vocabulary, no "It is important to note".
- Sentences vary in length; paragraphs state claims first and evidence second; the
  discussion concedes limitations in numbered, specific terms.
- Figures: four figures, each generated by assertion-guarded code from sealed artifacts
  only; captions name the artifact and hash location; layout checked for label/point
  overlap (Fig 1 label gutter at x = −7.6 vs min T_lo −1.98; Fig 3 row spacing 0.6).

## 5. Code-quality audit (release scripts)

- `b03_lung_certify.py`, `b03_lung_null.py`, `b03_lung_load.py`, `b03_lung_concordance.py`,
  `b03_lung_crop*.py`, `b03_paper_figures.py`: single-source-of-truth imports from the
  tissue-1 machinery (`b03_scout`, `b03_robust_null_v3`, `b03_certify`); env-var
  parameterization with tissue-1 defaults unchanged; assertions at every seam
  (alignment, label counts, tested-gene leaks, cap cross-check, exactness, subcube).
- The DEV-016 cap formula appears in exactly one form, used by certify, null, and the
  observed cross-check; an invalid assertion found during development was removed rather
  than silently kept.
- Spawn-bootstrap guard (`if __name__ == "__main__"`) present in multiprocessing drivers.
- Known and documented (not hidden): the null driver's startup banner prints "all 75"
  where 54 lung rows are evaluated; the check itself is `len(rows)` and passed on 54
  (verified in `lung_null.log` line 6 wording is cosmetic only).

## 6. Issues found and fixed in this audit

1. Certified-positive width stats (3.21/2.23–3.74 → 3.34/2.40–4.58).
2. Lung q_neg bound (0.036 → 0.027) in two places.
3. Lung/breast transcript-count basis (5.1× → 4.3× with explicit basis statement).
4. Figure 3 caption (omitted the PECAM1→KDR control block present in the plotted data).
5. Breast label-class description (19 → 20 panel classes, 15 transfer targets).
6. Proseg lung runtime (~1 h → ~1.7 h).
7. Salas et al. reference now cited in text.
8. Em-dash sweep and appositive recasts.

## 7. Residual known limitations (disclosed in the paper, not defects)

- Single section per tissue; certificates are per-section statements.
- One alternative segmentation engine (Proseg); the gates accept more without modification.
- The lung analysis is a pre-registered crop, not the full section.
- B = 1000 floors p-values at 9.99e-4; no verdict changes at larger B.
- The certified statistic is a co-expression product, not a physical-interaction measure.
