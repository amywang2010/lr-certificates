# Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty

Code, pre-registrations, and result artifacts for:

> **Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty in spatial transcriptomics**

The framework computes, for each ligand–receptor (LR) contrast in a spatial
transcriptomics section, an **interval valid over an explicit uncertainty set of
molecule-to-cell assignments** (boundary-local reassignment with capacity constraints,
mask erosion/dilation, per-gene leakage intervals). A direction is reported only when
both bounds exclude zero **and** a permutation null matched to the worst-case layer of
the statistic confirms it (within-region Benjamini–Hochberg FDR). Everything else is
reported as non-identifiable.

Two public Xenium sections are analyzed end-to-end with identical machinery:

| Tissue | Dataset | Panel | Result |
|---|---|---|---|
| Breast | Xenium v1.0.1 FFPE Human Breast Cancer Rep1 | 313-plex | 37 positive / 9 negative / 29 non-identifiable of 75 contrasts |
| Lung | Xenium Prime 5K FFPE Human Lung Cancer (pre-registered crop) | 5,001-gene | 1 positive (CDH1→EGFR) / 6 negative (ESR1↔PGR mutual exclusion) of 54 contrasts |

Cross-segmentation validation uses Proseg as an independent probe: coverage of certified
intervals (E2), cross-segmentation significance (E3), and per-certificate kernel margins
(E4). All design choices were frozen in pre-registration documents before any expression
computation; changes made after the initial pre-registration are documented in dated
amendments and per-phase reports (available on request, to be deposited with the journal
submission).

## Repository layout

```
paper/
  B03_MANUSCRIPT.md            full manuscript
  B03_SUPPLEMENT.md            supplementary methods and tables
  REFERENCES_PROVENANCE.md     per-reference verification record
  figures/                     Figure 1–4 (PNG 600 dpi + PDF), generated only from
                               sealed artifacts by src/b03_paper_figures.py
src/                           analysis code (single source of truth per step)
results/                       breast-section artifacts (CSV/JSON) + manifests
xenium_lung/crop/data/         lung-crop artifacts (CSV/JSON) + labels/regions
B03_THEORY.md                  exactness theorem, uncertainty-set construction
B03_PREREGISTRATION.md         frozen pre-registration (breast)
B03_PHASE2..10_PREREG*.md      frozen pre-registrations (cross-segmentation, null,
                               margins, lung replication, B=10,000 rerun, neighbor
                               effect, third segmentation leg) + addenda
(per-phase reports with receipts are deposited with the journal submission)
B03_DATA_AUDIT.md              dataset design audit
B03_SCOUT_MANIFEST.sha256      SHA-256 manifest, breast scout artifacts
results/PHASE2_MANIFEST.sha256 SHA-256 manifest, phase-2 artifacts
PHASE3/4/5_MANIFEST.sha256     SHA-256 manifests, later phases (PHASE5 supersedes
                               earlier manifests for files it covers)
```

## Reproduction

Hardware used: 8-core workstation, 15.7 GB RAM, Windows. Full-section Proseg on the lung
dataset is memory-infeasible at this scale; the lung analysis uses the pre-registered
expression-blind crop (selection rule in `B03_PHASE5_ADDENDUM.md` and
`results/lung_crop_window.json`).

1. Download the two public datasets (exact URLs and byte sizes are recorded in the
   pre-registration addenda).
2. Breast pipeline: `src/b03_load.py` → `src/b03_scout.py` → `src/b03_robust_null_v3.py`
   → `src/b03_finish.py` (labels, regions, certified intervals, matched null, referee pass).
3. Cross-segmentation: install Proseg 3.2.0, run the exact command in
   `B03_PHASE2_PREREG.md` Amendment A4; then `src/b03_phase2_vendor_prep.py`,
   `src/b03_phase2_analysis.py`, `src/b03_phase3_proseg_null.py`, `src/b03_phase4_kernel.py`.
4. Lung pipeline: `src/b03_lung_audit.py`, `src/b03_lung_crop.py`,
   `src/b03_lung_crop_extract.py`, `src/b03_lung_load.py`, `src/b03_lung_certify.py`,
   `src/b03_lung_null.py`, then `src/b03_phase5_e2.py`, `src/b03_phase5_e3.py`,
   `src/b03_phase5_vendor_prep.py` and the phase-4 driver with lung environment variables.
5. Figures: `python src/b03_paper_figures.py` (reads only sealed artifacts; assertions
   verify row counts before plotting).

Every result artifact is covered by a SHA-256 manifest; `RELEASE_MANIFEST.sha256` covers
this repository exactly as shipped.

## Statistical conventions

- Permutation p-values use the add-one convention (Phipson & Smyth 2010; North et al. 2002).
- BH false-discovery control is applied within region over the pre-specified pair list.
- The pseudocount ε = 0.5 was frozen pre-analysis; interval-level ε-sensitivity is in the
  supplement (S5).
- The worst-case layer of the statistic is what the null permutes; testing worst-case
  bounds against point-statistics nulls is trivially conservative, and the registered
  null recomputes the same worst-case objects under each permutation.

## Data availability

Both datasets are public (10x Genomics). The lung analysis region is a pre-registered
expression-blind spatial crop of the full section; crop boundaries and the selection rule
are in `B03_PHASE5_ADDENDUM.md` and `results/lung_crop_window.json`.
