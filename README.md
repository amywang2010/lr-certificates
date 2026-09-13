# Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty

Code, pre-registrations, and SHA-256-hashed result artifacts for the study:

> **Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty in spatial transcriptomics**

The study computes, for each ligand–receptor (LR) contrast in a spatial transcriptomics
section, an **interval valid over an explicit uncertainty set of molecule-to-cell
assignments** (boundary-local reassignment with capacity constraints, mask
erosion/dilation, per-gene leakage intervals). A biological direction is reported only
when both bounds exclude zero **and** a permutation null matched to the worst-case layer
of the statistic confirms it (within-region Benjamini–Hochberg FDR). Everything else is
reported as non-identifiable. The manuscript and supplement are provided with the
journal submission; this repository holds the analysis code and the artifacts every
number traces to.

## What the study finds

| Section | Dataset | Panel | Result |
|---|---|---|---|
| Breast (Rep1) | Xenium v1.0.1 FFPE Human Breast Cancer Rep1 (10x Genomics) | 313-plex | 37 positive / 9 negative / 29 non-identifiable of 75 contrasts |
| Lung (crop) | Xenium Prime 5K FFPE Human Lung Cancer (10x Genomics; pre-registered expression-blind crop) | 5,001-gene | 1 positive (CDH1→EGFR) / 6 negative (ESR1↔PGR mutual exclusion) of 54 contrasts |
| Breast (second section) | FFPE breast cancer, Xenium Prime 5K panel (public EBI archive repackage S-BIAD2146; pre-registered expression-blind crop) | 5,000-gene | 0 certified of 30 evaluable contrasts — all reported non-identifiable, with Proseg interval coverage 30/30 |

The third section isolates a mechanism the two-tissue comparison could not: on the same
tissue and the same 10 pairs, certification collapses (22/30 → 0/30) while certified
intervals get *narrower* (median width 2.60 → 1.60 log2 units). What governs
identifiability is the tested genes' share of transcript mass (8.1% of the 313-plex
section's transcripts, 0.4% of the 5,000-plex section's), which compresses effect sizes
below the certification floor.

Cross-segmentation validation uses Proseg as an independent probe on all three sections:
coverage of certified intervals (E2), cross-segmentation significance (E3), and
per-certificate kernel margins (E4). All design choices were frozen in pre-registration
documents before any expression computation; changes made after the initial
pre-registration are documented in dated addenda, and the per-phase reports are
deposited with the journal submission.

## Repository map

| Path | Content |
|---|---|
| `preregistration/` | All frozen pre-registrations: breast design (`B03_PREREGISTRATION.md`), cross-segmentation through the third-segmentation leg (`B03_PHASE2..10_PREREG*.md`, `B03_PHASE5_ADDENDUM.md`, `B03_PHASE7_ADDENDUM_A1.md`), plus `B03_THEORY.md` (bound construction and exactness argument) and `B03_DATA_AUDIT.md` (dataset identity and audit gates) |
| `src/` | Analysis code, one script per step: load and concordance (`b03_load.py`, `b03_scout.py`), certified bounds (inline in `b03_scout.py`), matched null (`b03_robust_null_v3.py`), B=10,000 chain (`b03_phase6_*.py`), cross-segmentation and kernel margins (`b03_phase2..5_*.py`), lung pipeline (`b03_lung_*.py`), second-breast-section pipeline (`b03_s6_*.py`), neighbor effect (`b03_phase9_neighbor.py`), figures (`b03_paper_figures.py`) |
| `results/` | Breast-section artifacts: `scout_final_matched.csv` (75 rows, the source of Figures 1 and 3), `null_robust_matched*.csv`, `b10k_breast/` and `b10k_lung/` (B=10,000 reruns), `phase6/` (determinism gate), `exactness_verification.csv`, `sensitivity/` (independent re-derivation of all intervals), per-phase manifests `*.sha256` |
| `lung/` | Lung-crop artifacts, same schema: `scout_final_lung.csv` (54 rows), null outputs, E2/E3/E4 summaries |
| `breast_s6/` | Second-breast-section artifacts, same schema: `scout_final_s6.csv` (30 rows), matched-null output at B=10,000 (`null_robust_matched_final.csv`), E2 coverage (`phase7_e2_coverage.csv`, `phase7_e2_summary.json`), corrected-rule re-derivation (`s6_tight_bounds_sensitivity.csv`, `s6_tight_bounds_summary.json`), verification receipts (`exactness_verification.csv`, `subcube_check.json`, `capacity_check.json`, `type_absent_skips.json`, `concordance.json`) |
| `RELEASE_MANIFEST.sha256` | Hashes this repository exactly as shipped; verifies byte-for-byte on any platform |

## Provenance

Every number in the manuscript traces to a SHA-256-hashed artifact in `results/`,
`lung/`, or `breast_s6/` (per-phase manifests record artifact and source hashes as they
existed at analysis time). Loaders reproduce each vendor assignment from molecules and
polygons and verify it before any statistic is computed: per-gene bitwise or
totals-level gates per package (breast Rep1 exact-per-cell ≥ 0.90 with per-gene r ≥ 0.958
over 40 audited genes; lung 20/20 audited genes bitwise-exact; second breast section
totals-level gates, r = 0.999, per the packaging addendum). The released code differs
from the analysis-time files only in comment text (internal record identifiers removed);
byte-identical analysis-time code and the dated addenda are deposited with the journal as
supplementary material. The interval construction was independently re-derived after the
initial analysis; the re-derivation reproduces every certified verdict and shows the
reported intervals are conservative supersets (`results/sensitivity/`, and for the third
section `breast_s6/s6_tight_bounds_sensitivity.csv`).

## Reproduction

Hardware used: 8-core workstation, 15.7 GB RAM, Windows. Full-section Proseg on the lung
and 5,000-plex datasets is memory-infeasible at this scale; both use the pre-registered
expression-blind crop (selection rule in `preregistration/B03_PHASE5_ADDENDUM.md`,
window records in `results/lung_crop_window.json` and the Phase 7 addendum).

1. Download the three public datasets (exact URLs and byte sizes are recorded in the
   pre-registration documents and the Phase 7 addendum).
2. Breast pipeline: `src/b03_load.py` → `src/b03_scout.py` → `src/b03_robust_null_v3.py`
   → `src/b03_finish.py`.
3. Cross-segmentation: install Proseg 3.2.0, run the command registered in
   `preregistration/B03_PHASE2_PREREG.md` Amendment A4; then
   `src/b03_phase2_vendor_prep.py`, `src/b03_phase2_analysis.py`,
   `src/b03_phase3_proseg_null.py`, `src/b03_phase4_kernel.py`.
4. Lung pipeline: `src/b03_lung_audit.py` → `src/b03_lung_crop.py` →
   `src/b03_lung_crop_extract.py` → `src/b03_lung_load.py` → `src/b03_lung_certify.py`
   → `src/b03_lung_null.py`, then the E2/E3/E4 scripts under `src/b03_phase5_*.py`.
5. Second breast section (5,000-plex): `src/b03_s6_extract.py` → `src/b03_s6_prep.py`
   (frozen expression-blind crop) → `src/b03_s6_load.py` → `src/b03_s6_certify.py` →
   `src/b03_s6_null.py` → `src/b03_s6_finalize.py` → `src/b03_s6_vendor_prep.py` →
   `src/b03_s6_e2.py` → `src/b03_s6_tight.py` (the corrected-rule re-derivation). Proseg
   input construction is documented in the Phase 7 addendum.
6. Figures: `python src/b03_paper_figures.py` (reads only sealed artifacts; assertions
   verify row counts before plotting).

The permutation nulls checkpoint every 50 permutations and resume deterministically.
Expected runtimes: certified bounds ~2 min; B=10,000 matched null ~12 h on 4 workers
(second section: ~12 min; the 30-row surface is smaller); cross-segmentation battery
< 1 h per section.

## Statistical conventions

- Permutation p-values use the add-one convention (Phipson & Smyth 2010; North et al. 2002).
- BH false-discovery control is applied within region over the pre-specified pair list.
- The pseudocount ε = 0.5 was frozen pre-analysis; interval-level ε-sensitivity is in the
  supplement.
- The null recomputes the same worst-case objects as the observed statistic under each
  permutation; testing worst-case bounds against point-statistics nulls would be
  trivially conservative.

## Data availability

All three datasets are public. Breast Rep1 and the lung section are 10x Genomics Xenium
bundles (Xenium FFPE Human Breast Cancer Rep1, v1.0.1; Xenium Prime 5K FFPE Human Lung
Cancer). The second breast-cancer section is distributed as the EBI's public zarr
repackage (accession S-BIAD2146); its packaging differences from the 10x bundles are
arbitrated item by item in `preregistration/B03_PHASE7_ADDENDUM_A1.md`. The lung and
second-breast analysis regions are pre-registered expression-blind spatial crops; crop
boundaries and selection rules are in the addenda and the window records above.
