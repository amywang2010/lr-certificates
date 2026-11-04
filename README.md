# Joint Attainability of Worst-Case Co-expression Bounds via Capacitated Flows

Code and hashed result artifacts for the RECOMB 2027 paper of the same title
(Amy Wang and Manolis Kellis, MIT).

The repository computes an interval for a ligand–receptor co-expression score over a
declared set of plausible transcript-to-cell assignments, tests whether one assignment
attains both per-gene optima of that interval, and certifies a direction with a
permutation null that recomputes the same worst-case object under every label
permutation.

## What the code computes

**Interval ends.** The score separates into a ligand term and a receptor term. Each
upper endpoint is a maximum flow on a network whose nodes are source cells, donor-set
signatures and counted cells. Each lower endpoint has a closed form. `src/b03_exact_end.py`
solves the 150 breast gene-ends this way and reads every flow back into an assignment.
With source-loss budgets inactive, a Hall condition on the counted cells decides when a
per-transcript count already attains the flow value, and the shortfall otherwise equals
the maximum deficiency.

**Joint attainment.** Two further flows on a merged network decide whether a single
assignment reaches both per-gene ends. `src/b03_coupled_end.py` runs the test on all three
sections and reports the gap when it fails.

**Inference.** `src/b03_matched_null.py` recomputes the interval under each
within-region, decile-matched label permutation and writes add-one Monte Carlo p-values.
The claim layer reads those stored values and doubles the p-value selected by the
interval sign before Benjamini–Yekutieli correction within region
(`src/b03_claim_layer.py`). A test carries a claim when the interval excludes zero and
the adjusted q-value is at most 0.10.

## Results in the released artifacts

| Section | Tests | Interval sign (+/−/0) | Claims (+) | Claims (−) |
|---|---|---|---|---|
| Breast Rep1, 313-plex | 75 | 39 / 12 / 24 | 36 | 9 |
| Lung crop, 5K panel | 54 | 1 / 32 / 21 | 1 | 4 |
| Breast, second section, 5K panel | 30 | 0 / 17 / 13 | 0 | 0 |

The joint-attainment test covers 157 of the 159 tests. Two lung tests have no movable
transcript for either gene and no joint-flow instance. The lower end is attained on all
157, the upper end on 149.

## Contents

| Path | Content |
|---|---|
| `src/` | Analysis scripts, one per pipeline step |
| `results/` | Breast artifacts, exact-end tables, conventional-null comparison, high-permutation reruns, intervention study, tested-gene mass share |
| `lung/` | Lung-crop artifacts under the same column schema |
| `breast_s6/` | Second-section artifacts under the same column schema |
| `figures/` | Figure 1, Figure S2 and the workflow schematic, as PDF and PNG |
| `tests/` | Release checks that run on the shipped artifacts |
| `requirements.txt` | Python environment |
| `CITATION.cff` | Citation metadata |
| `RELEASE_MANIFEST.sha256` | SHA-256 of every other file in the repository |
| `LICENSE` | MIT |

Each section directory holds a 75-, 54- or 30-row table (`scout_final_matched.csv`), its
matched null, the claim layer (`claim_layer.csv`), feasibility and subcube receipts, and
the cross-segmentation summaries.

## Installation

```
python -m pip install -r requirements.txt
```

The analysis is CPU-only. Proseg 3.2.0 is a separate tool, invoked on the molecule
files for the cross-segmentation checks.

## Quick start from the shipped artifacts

These steps need no vendor data and finish in seconds to minutes.

```
python src/b03_claim_layer.py      # claim columns, claim_layer.csv per section
python tests/test_release.py       # release checks
python -m unittest discover -s tests
```

Verify every shipped file against the manifest:

```
sha256sum -c RELEASE_MANIFEST.sha256
```

## Full reproduction

Scripts resolve paths from the repository root. Breast scripts read the vendor bundle
from `data/`. The lung and second-section pipelines read a workspace directory next to
the clone, set with `B03_WORKSPACE` (default: the parent directory), containing
`xenium_lung/` and `xenium_breast_s6/` with the vendor archive, an `extracted/`
subdirectory and a `crop/` subdirectory.

1. Breast intervals and matched null: `src/b03_load.py`, `src/b03_scout.py`,
   `src/b03_certify.py`, `src/b03_matched_null.py`, `src/b03_finalize.py`.
2. Interval verification and sensitivity: `src/b03_coverage.py`,
   `src/b03_phase9_neighbor.py`, `src/b03_exact_end.py`. The exhaustive check of the
   flow ends and the attainment criteria on 3,000 small instances is
   `src/b03_enumerate_check.py`, which writes `results/enumeration_check.json`.
3. Claim layer: `src/b03_claim_layer.py`.
4. Cross-segmentation on breast: run Proseg 3.2.0 on the molecule file, then
   `src/b03_phase2_vendor_prep.py`, `src/b03_phase2_label_validation.py`,
   `src/b03_phase2_analysis.py`, `src/b03_phase3_proseg_null.py`,
   `src/b03_phase4_kernel.py`.
5. Lung crop: `src/b03_lung_audit.py`, `src/b03_lung_crop.py`,
   `src/b03_lung_crop_extract.py`, `src/b03_lung_load.py`, `src/b03_lung_certify.py`,
   `src/b03_lung_null.py`, `src/b03_lung_concordance.py`, `src/b03_phase5_vendor_prep.py`,
   `src/b03_phase5_e2.py`, `src/b03_phase5_e3.py`. Shared settings are in
   `src/b03_lung_config.py`.
6. High-permutation reruns: `src/b03_phase6_breast.py`, `src/b03_phase6_lung.py`,
   `src/b03_phase6_analysis.py`.
7. Second section: `src/b03_s6_extract.py`, `src/b03_s6_prep.py`, `src/b03_s6_load.py`,
   `src/b03_s6_certify.py`, `src/b03_s6_null.py`, `src/b03_s6_finalize.py`,
   `src/b03_s6_vendor_prep.py`, `src/b03_s6_e2.py`, `src/b03_s6_tight.py`.
8. Conventional-null comparison: `src/b03_baseline_conventional_null.py`,
   `src/b03_baseline_comparison.py`.
9. Tested-gene mass share: `src/b03_mass_share.py`.
10. Intervention study: `python src/b03_synthstudy.py --mode verify`, then `--mode run`
    for each of the 20 runs, then `python src/b03_synthstudy_analysis.py`.
11. Figures: `python src/b03_paper_figures.py` and `python src/b03_schematic.py`.

Permutation nulls checkpoint every 50 permutations and resume deterministically. On the
reference workstation (8 cores, 15.7 GB RAM, Windows) the exact-end table takes 540 s,
a matched null at B = 1000 takes 74 min on the breast section, and the B = 10,000 rerun
takes about 12 h.

## Data

All three sections are public and are not redistributed here.

| Section | Source |
|---|---|
| Breast Rep1 | 10x Genomics Xenium FFPE Human Breast Cancer Rep1, v1.0.1 |
| Lung crop | 10x Genomics Xenium Prime 5K FFPE Human Lung Cancer, technical note release |
| Breast, second section | EBI zarr repackage of the STHELAR release, BioStudies S-BIAD2146 |

The lung and second-section analyses use expression-blind crops. Window coordinates are
recorded in `results/lung_crop_window.json` and `breast_s6/s6_crop_window.json`.

## Tests

`tests/test_release.py` runs without vendor data. It checks the released row tables
against the counts reported in the paper, recomputes the claim layer from the stored
permutation values, and checks the Benjamini–Hochberg and Benjamini–Yekutieli routines
on fixed inputs.

```
python -m unittest discover -s tests -v
```

## Citation

Cite the paper and this repository. Metadata is in `CITATION.cff`. Code is released
under the MIT licence.
