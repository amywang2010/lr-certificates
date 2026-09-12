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

| Tissue | Dataset | Panel | Result |
|---|---|---|---|
| Breast | Xenium v1.0.1 FFPE Human Breast Cancer Rep1 | 313-plex | 37 positive / 9 negative / 29 non-identifiable of 75 contrasts |
| Lung | Xenium Prime 5K FFPE Human Lung Cancer (pre-registered expression-blind crop) | 5,001-gene | 1 positive (CDH1→EGFR) / 6 negative (ESR1↔PGR mutual exclusion) of 54 contrasts |

Cross-segmentation validation uses Proseg as an independent probe: coverage of certified
intervals (E2), cross-segmentation significance (E3), and per-certificate kernel margins
(E4). All design choices were frozen in pre-registration documents before any expression
computation; changes made after the initial pre-registration are documented in dated
amendments and per-phase reports (available on request, deposited with the journal
submission).

## Repository map

| Path | Content |
|---|---|
| `preregistration/` | All frozen pre-registrations: breast design (`B03_PREREGISTRATION.md`), cross-segmentation through neighbor-effect and third-segmentation legs (`B03_PHASE2..10_PREREG*.md`, `B03_PHASE5_ADDENDUM.md`), plus `B03_THEORY.md` (bound construction and exactness argument) and `B03_DATA_AUDIT.md` (dataset identity and audit gates) |
| `src/` | Analysis code, one script per step: load and concordance (`b03_load.py`, `b03_scout.py`), certified bounds (inline in `b03_scout.py`), matched null (`b03_robust_null_v3.py`), B=10,000 chain (`b03_phase6_*.py`), cross-segmentation and kernel margins (`b03_phase2..5_*.py`), lung pipeline (`b03_lung_*.py`), figures (`b03_paper_figures.py`) |
| `results/` | Breast-section artifacts: `scout_final_matched.csv` (75 rows, the source of Figures 1 and 3), `null_robust_matched*.csv`, `b10k_breast/` (B=10,000 rerun + checkpoints), `exactness_verification.csv`, `sensitivity/` (independent re-derivation of all intervals), per-phase manifests `*.sha256` |
| `lung/` | Lung-crop artifacts, same schema: `scout_final_lung.csv` (54 rows), null outputs, E2/E3/E4 summaries, manifests |
| `RELEASE_MANIFEST.sha256` | Hashes this repository exactly as shipped; verifies byte-for-byte on any platform |

## Provenance

Every number in the manuscript traces to a SHA-256-hashed artifact in `results/` or
`lung/` (per-phase manifests record artifact and source hashes as they existed at
analysis time). The released code differs from the analysis-time files only in comment
text; byte-identical analysis-time code and the dated amendment records are deposited
with the journal as supplementary material. The interval construction was independently
re-derived after the initial analysis; the re-derivation reproduces every certified
verdict and shows the reported intervals are conservative supersets
(`results/sensitivity/`).

## Reproduction

Hardware used: 8-core workstation, 15.7 GB RAM, Windows. Full-section Proseg on the lung
dataset is memory-infeasible at this scale; the lung analysis uses the pre-registered
expression-blind crop (selection rule in `preregistration/B03_PHASE5_ADDENDUM.md`,
window record in `results/lung_crop_window.json`).

1. Download the two public datasets (exact URLs and byte sizes are recorded in the
   pre-registration documents).
2. Breast pipeline: `src/b03_load.py` → `src/b03_scout.py` → `src/b03_robust_null_v3.py`
   → `src/b03_finish.py`.
3. Cross-segmentation: install Proseg 3.2.0, run the command registered in
   `preregistration/B03_PHASE2_PREREG.md` Amendment A4; then
   `src/b03_phase2_vendor_prep.py`, `src/b03_phase2_analysis.py`,
   `src/b03_phase3_proseg_null.py`, `src/b03_phase4_kernel.py`.
4. Lung pipeline: `src/b03_lung_audit.py` → `src/b03_lung_crop.py` →
   `src/b03_lung_crop_extract.py` → `src/b03_lung_load.py` → `src/b03_lung_certify.py`
   → `src/b03_lung_null.py`, then the E2/E3/E4 scripts under
   `src/b03_phase5_*.py`.
5. Figures: `python src/b03_paper_figures.py` (reads only sealed artifacts; assertions
   verify row counts before plotting).

The permutation nulls checkpoint every 50 permutations and resume deterministically.
Expected runtimes: certified bounds ~2 min; B=10,000 matched null ~12 h on 4 workers;
cross-segmentation battery < 1 h.

## Statistical conventions

- Permutation p-values use the add-one convention (Phipson & Smyth 2010; North et al. 2002).
- BH false-discovery control is applied within region over the pre-specified pair list.
- The pseudocount ε = 0.5 was frozen pre-analysis; interval-level ε-sensitivity is in the
  supplement.
- The null recomputes the same worst-case objects as the observed statistic under each
  permutation; testing worst-case bounds against point-statistics nulls would be
  trivially conservative.

## Data availability

Both datasets are public (10x Genomics): Xenium FFPE Human Breast Cancer Rep1
(v1.0.1) and Xenium Prime 5K FFPE Human Lung Cancer. The lung analysis region is a
pre-registered expression-blind spatial crop; crop boundaries and the selection rule are
in `preregistration/B03_PHASE5_ADDENDUM.md` and `results/lung_crop_window.json`.
