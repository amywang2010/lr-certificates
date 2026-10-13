# Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty

Code and SHA-256-hashed result artifacts for the study:

> **Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty in spatial transcriptomics**

For each ligand–receptor contrast in a spatial transcriptomics section the study computes an
interval valid over an explicit uncertainty set of molecule-to-cell assignments (boundary-local
reassignment with per-cell capacity constraints, mask erosion and dilation, and a parameterised
donor radius). A biological direction is reported only when both bounds exclude zero and a
permutation null matched to the worst-case layer confirms it, with within-region
Benjamini–Hochberg control. Everything else is reported as non-identifiable. This repository
holds the reproduction materials: the analysis code, the frozen configuration, and the hashed
result artifact behind every reported number.

## What the study finds

| Section | Dataset | Panel | Result |
|---|---|---|---|
| Breast (Rep1) | Xenium v1.0.1 FFPE Human Breast Cancer Rep1 (10x Genomics) | 313-plex | 37 positive, 9 negative, 29 non-identifiable of 75 contrasts |
| Lung (crop) | Xenium Prime 5K FFPE Human Lung Cancer (10x Genomics, expression-blind crop) | 5,001 genes | 1 positive (CDH1→EGFR) and 6 negative (ESR1↔PGR) of 54 contrasts |
| Breast (second section) | FFPE breast cancer on the Xenium Prime 5K panel (public EBI repackage, BioStudies S-BIAD2146, expression-blind crop) | 5,000 genes | 0 certified of 30 evaluable contrasts, all reported non-identifiable |

The third section isolates a mechanism the two-tissue comparison cannot. On the same tissue and
the same ten pairs, certification collapses from 22 of 30 rows to 0 of 30 while the certified
intervals get narrower, with median width falling from 2.60 to 1.60 log2 units. The tested genes
carry 8.08% of the 313-plex section's transcripts and 0.36% of the 5,000-plex section's, which
compresses the effect sizes below the certification floor. The released computation is
`results/tested_mass_share.json`, written by `src/b03_mass_share.py` from the vendor molecule
files.

A 20-run controlled intervention study on the first section (paper Section 5.6) tests that
mechanism under intervention. Injected-pair certification transitions between 1% and 4% ambient
tested-gene mass share and is seed-replicable; the smallest injection that certifies is 5.3
thousand transcripts on a 4%-mass substrate against 87.6 thousand on the full 13.84%-mass
section. Per-run artifacts and the fleet summary are under `results/synthstudy/`, and
`src/b03_synthstudy_analysis.py` recomputes every Section 5.6 number from them, including the
dose and power curves (`results/synthstudy/dose_curve.csv`, `power_curve.csv`, `fcr_table.csv`).

The conventional co-expression permutation pipeline of the CellChat/CellPhoneDB family, run on
the same breast section, calls 73 of 75 contrasts interactions. The certificates confirm 39 of
those calls in the same direction, render 28 non-identifiable, and reverse 6, all of them
ESR1–PGR rows (`results/baseline_comparison/`).

Cross-segmentation validation uses Proseg as an independent probe on all three sections:
interval coverage of certified rows (`E2`), certified-direction significance under a Proseg-side
matched null (`E3`), and per-certificate kernel margins (`E4`). Regions, crop windows, gates and
thresholds are fixed before the corresponding statistic is computed and are listed in the paper
and in this README.

## Repository map

| Path | Content |
|---|---|
| `src/` | Analysis code, one script per step; the reproduction chain below lists the order |
| `results/` | Breast-section artifacts: the 75-row matched table and its null, the B=10,000 reruns, the determinism gate, the feasibility and subcube receipts, the independent re-derivation under `sensitivity/`, the conventional-pipeline comparison under `baseline_comparison/`, the intervention study under `synthstudy/`, and the tested-gene mass share |
| `lung/` | Lung-crop artifacts under the same schema: the 54-row table, its matched null, the feasibility, subcube and capacity receipts, and the E2/E3/E4 summaries |
| `breast_s6/` | Second-section artifacts under the same schema: the 30-row table, its matched null at B=10,000, the feasibility, subcube and capacity receipts, the E2 coverage summary, and the corrected-rule re-derivation |
| `figures/` | The three figures in PNG and PDF, drawn at the LNCS text width |
| `RELEASE_MANIFEST.sha256` | SHA-256 of every other file in the repository, LF-terminated and sorted by path. Verify with `sha256sum -c RELEASE_MANIFEST.sha256` on Linux, or `shasum -a 256 -c RELEASE_MANIFEST.sha256` where available |

## Figures

`python src/b03_paper_figures.py` writes all three figures into `figures/` from released
artifacts only, asserting row counts before plotting.

| Figure | File | Panels | Source artifacts |
|---|---|---|---|
| 1 | `fig1_main.pdf` | (a) certified intervals for the 25 breast pairs across three regions, (b) cross-tissue certificates for the shared pairs | the 75-row matched table, and the 18-row cross-tissue table `results/paper_fig3_data.csv` |
| 2 | `fig2_support.pdf` | (a) cross-segmentation interval coverage against the gate, (b) kernel margins c* per certified row | the coverage summaries of the breast and lung legs, and the two kernel-margin tables (46 and 7 rows) |
| S1 | `fig0_workflow.pdf` | the construction: both inputs enter one interval machinery | drawn from no data, by `src/b03_schematic.py` |

The schematic script measures every label's rendered bounding box against its shape and checks
all shape pairs for overlap, aborting instead of writing a figure that fails either check.

## Reproduction

All scripts resolve paths from the repository root, independently of the working directory.
Breast-section scripts read the vendor bundle from `data/` and write intermediates to `data/` and
`results/`. The lung and second-section pipelines operate on a workspace directory next to the
repository clone, its parent by default, overridable with `B03_WORKSPACE`: the lung workspace is
`xenium_lung/` (the vendor zip, an `extracted/` subdirectory, and a `crop/` subdirectory), and
the second-section workspace is `xenium_breast_s6/` (the S-BIAD2146 zarr zip, `extracted/`,
`crop/`, `chunks/`). Vendor bundles are not redistributed; everything they produce that a
reported number depends on is released under `results/`, `lung/` and `breast_s6/`.

1. Download the three public datasets listed under Data availability and place the breast bundle
   files in `data/`, the lung bundle in `xenium_lung/`, and the second-section zarr in
   `xenium_breast_s6/`.
2. Breast certified bounds and matched null: `src/b03_load.py` → `src/b03_scout.py` →
   `src/b03_matched_null.py` → `src/b03_finalize.py`. The last step writes the two tables the
   breast results are reported from, `results/scout_final_matched.csv` (75 rows) and
   `results/null_robust_matched_final.csv`, and records its consistency battery in
   `finalize_checks.md`.
3. Verification and sensitivity on the breast leg: `src/b03_coverage.py` (held-out coverage on a
   second assignment path, and the nucleus-only sensitivity), and `src/b03_phase9_neighbor.py`
   (the independent re-derivation of all 159 intervals under the corrected per-gene-slot rule,
   writing `results/sensitivity/`, plus the neighbour-effect check).
4. Cross-segmentation on the breast leg: install Proseg 3.2.0 and run
   `proseg transcripts.parquet` on the identical molecule file, then
   `src/b03_phase2_vendor_prep.py`, `src/b03_phase2_label_validation.py`,
   `src/b03_phase2_analysis.py`, `src/b03_phase3_proseg_null.py`, `src/b03_phase4_kernel.py`.
5. Lung leg: `src/b03_lung_audit.py` → `src/b03_lung_crop.py` →
   `src/b03_lung_crop_extract.py` → `src/b03_lung_load.py` → `src/b03_lung_certify.py` →
   `src/b03_lung_null.py` → `src/b03_lung_concordance.py` → `src/b03_phase5_vendor_prep.py` →
   `src/b03_phase5_e2.py` → `src/b03_phase5_e3.py`. The configuration shared with the
   second-section leg lives in `src/b03_lung_config.py`.
6. High-permutation reruns and the determinism gate: `src/b03_phase6_lung.py`,
   `src/b03_phase6_breast.py`, then `src/b03_phase6_analysis.py`.
7. Second section on the 5,000-plex panel: `src/b03_s6_extract.py` → `src/b03_s6_prep.py` →
   `src/b03_s6_load.py` → `src/b03_s6_certify.py` → `src/b03_s6_null.py` →
   `src/b03_s6_finalize.py` → `src/b03_s6_vendor_prep.py` → `src/b03_s6_e2.py` →
   `src/b03_s6_tight.py`.
8. Conventional-pipeline comparison: `src/b03_baseline_conventional_null.py` →
   `src/b03_baseline_comparison.py`.
9. Tested-gene mass share: `src/b03_mass_share.py` (needs the breast intermediates from step 2
   and the second-section crop from step 7).
10. Intervention study: `python src/b03_synthstudy.py --mode verify` (invariants),
    `--mode run` for each of the 20 runs (arm, seed and measured injection are recorded in each
    `run_manifest.json`), then `python src/b03_synthstudy_analysis.py` (writes the dose and power
    curves and the fleet summary). Fleet wall time is about 2.5 h on the reference workstation,
    and results are deterministic given the recorded seeds.
11. Figures: `python src/b03_paper_figures.py`.

The permutation nulls checkpoint every 50 permutations and resume deterministically. Expected
runtimes on the reference workstation: certified bounds about 2 minutes, a 75-minute matched null
at B=1,000 over 4 workers, about 12 hours at B=10,000 on the breast leg and 12 minutes on the
second section, and under an hour per section for the cross-segmentation battery. The vendor load
takes about 4 minutes for breast, 177 s for lung and 135 s for the second section. Proseg takes
about 2.4 h on the full breast section, about 1.7 h on the lung crop and 2 h 07 m on the
second-section crop, at a peak of 6.9 GB private memory on the lung crop.

Hardware: 8 physical cores, 15.7 GB RAM, Windows. Full-section Proseg on the lung and
second-section datasets is memory-infeasible at this scale, so both legs use an expression-blind
crop: 3.5 × 3.5 mm candidate windows on a 0.5 mm grid over the full section bounds, ranked by
vendor-cell count with no access to expression data, selecting the median-density window with
ties broken lexicographically. The frozen window records are `results/lung_crop_window.json` and
`breast_s6/s6_crop_window.json`.

## Verification and provenance

Every number in the manuscript is backed by a hashed artifact, and the checks below are the
receipts for the verification claims made in the paper.

| Claim | Receipt |
|---|---|
| The reproduced assignment matches the vendor matrix | `results/concordance.json`, `results/concordance_diagnosis.json` (per-cell exact 0.9146, per-gene Pearson median 0.9912 and minimum 0.9578 over 40 audited genes), `lung/concordance.json` (0.9568, 20 audited genes bitwise exact), `breast_s6/concordance.json` (totals-level gate, r = 0.9991) |
| Bounds are never violated | `exactness_verification.csv` in each of the three legs: 150,000 + 108,000 + 60,000 randomized feasible configurations, 0 violations |
| The interval bounds are attained | `results/subcube_check.json` and the two leg copies: exhaustive enumeration of a subcube matches the greedy extremes |
| Both interval code paths agree | `results/null_path_validation.csv`: bitwise equality of the shipped and fast paths at every null startup |
| The null layers agree with the observed bounds | `scout_summary_matched.json` in each leg: all rows within 1e-6 of the released bounds |
| All 159 intervals re-derive under the corrected rule | `results/sensitivity/` (75 breast and 54 lung rows, each with a `.sha256`), `breast_s6/s6_tight_bounds_sensitivity.csv` and `s6_tight_bounds_summary.json`: 0 verdict flips, median width tightening 9.2%, 14.1% and 10.4%. `src/b03_phase9_neighbor.py` recomputes the corrected bounds on the breast and lung legs with its own aggregation path and asserts bit-exact agreement with the released tables on every row and gene, and `src/b03_s6_tight.py` writes the second-section tables |
| The capacity rule binds where the paper says it does | `lung/capacity_check.json` (197 cells flagged) and `breast_s6/capacity_check.json` (1,350 cells flagged), with the capped formula applied inside the interval construction |
| The higher-permutation rerun changes no verdict | `results/phase6/phase6_summary.json` and `results/phase6/determinism_gate.json`: 0 flips across 75 breast and 54 lung rows at B=10,000 |
| Section 5.6 numbers | the 20 `results/synthstudy/*/run_manifest.json` files and the five tables written by `src/b03_synthstudy_analysis.py`, receipted by `results/synthstudy/synthstudy_analysis.sha256` |
| The tested-gene mass share | `results/tested_mass_share.json`, written by `src/b03_mass_share.py` from the vendor molecule files |

Two conventions are worth stating because the released tables expose them. First, the interval
layer and the certified layer are different: on the breast leg 39 of 75 rows are interval-positive
and 12 interval-negative, and the matched null then leaves 37 and 9. Second, the p-value floor at
B permutations is 1/(B+1), which is 9.99e-4 at B=1,000 and 9.999e-5 at B=10,000; rows at that
floor are reported there.

## Statistical conventions

- Permutation p-values use the add-one convention (Phipson and Smyth 2010, North et al. 2002).
- Benjamini–Hochberg is applied within region over the pair list. Benjamini–Yekutieli, valid
  under arbitrary dependence, changes 1 of 159 verdicts.
- The pseudocount is 0.5 throughout, and the interval-level sensitivity to it is tabulated in
  the appendix.
- The null recomputes the same worst-case objects as the observed statistic under each
  permutation. Testing worst-case bounds against a null of point statistics would be trivially
  conservative, since every permutation's worst-case minimum is at most its point value.

## Data availability

All three datasets are public. Breast Rep1 and the lung section are 10x Genomics Xenium bundles
(Xenium FFPE Human Breast Cancer Rep1 at v1.0.1, and Xenium Prime 5K FFPE Human Lung Cancer from
the post-Xenium technical note). The second breast-cancer section is distributed as the EBI's
public zarr repackage under BioStudies accession S-BIAD2146, from the STHELAR release. Its
packaging differences from the 10x bundles, namely the zarr layout and the deprecated codeword
categories, are handled explicitly in `src/b03_s6_extract.py` and `src/b03_s6_load.py`. The
lung and second-section analysis regions are expression-blind crops whose boundaries are the
window records named above.
