# Joint Attainability of Worst-Case Co-expression Bounds via Capacitated Flows

Code and hashed result artifacts for the RECOMB 2027 paper by Amy Wang and
Manolis Kellis (MIT). The analysis bounds a ligand–receptor co-expression score
over a declared set of transcript-to-cell assignments, tests whether one
assignment attains both per-gene optima, and certifies a direction with a
permutation null that recomputes the same bound.

## Layout

| Path | Contents |
|---|---|
| `src/` | One script per analysis step |
| `preprocessing/` | Crop windows, point-in-polygon concordance, concordance diagnosis |
| `results/breast/` | 75-test breast table, claim layer, exact-flow endpoints, joint attainment |
| `results/lung/` | 54-test lung-crop table and joint attainment |
| `results/breast_s6/` | 30-test second-section table and joint attainment |
| `results/coupled_end_summary.json` | Joint-attainment counts across the three sections |
| `coverage/` | Proseg coverage, Proseg-side null, label transfer, nucleus-only mask |
| `margins/` | Critical multiplier `c*` for certified tests |
| `sensitivity/` | Corrected-bound and pseudocount sensitivity tables |
| `provenance/` | B = 10,000 stability run, determinism gate, per-run tables |
| `intervention/` | 20-run injection study, arm tables, slot identity audit |
| `baseline_comparison/` | Conventional point-null comparison |
| `figures/` | Figure 1, Figures S1 and S2, and the panel B data |
| `tests/` | Release checks that run on the shipped artifacts |

## Reproducing the published numbers

Install the environment once. Every step below runs from the repository root.

```
python -m pip install -r requirements.txt
```

The shipped tables support the paper without vendor data. These commands
recompute the derived tables and check them against the stored values.

```
python src/b03_claim_layer.py        # claim columns in all three section tables
python src/b03_eps_sensitivity.py    # Table S5
python -m unittest discover -s tests -v
sha256sum -c RELEASE_MANIFEST.sha256
```

The heavy pipelines need the public sections. Breast scripts read the vendor
bundle from `data/`. The lung and second-section scripts read a workspace next
to the clone, set by `B03_WORKSPACE` (default: the parent directory), holding
`xenium_lung/` and `xenium_breast_s6/` with the vendor archive, an `extracted/`
subdirectory and a `crop/` subdirectory.

Intervals, matched null and claims:

```
python src/b03_load.py
python src/b03_scout.py
python src/b03_certify.py
python src/b03_matched_null.py
python src/b03_finalize.py
python src/b03_claim_layer.py
```

Endpoints and joint attainment:

```
python src/b03_exact_end.py          # 150 breast gene-ends, 540 s
python src/b03_coupled_end.py        # three sections, writes results/coupled_end_summary.json
python src/b03_enumerate_check.py    # 3,000 small instances
```

Evidence beyond the declared uncertainty set:

```
python src/b03_coverage.py           # nucleus-only mask, coverage/coverage_summary.json
python src/b03_phase2_vendor_prep.py
python src/b03_phase2_label_validation.py
python src/b03_phase2_analysis.py    # coverage/phase2_proseg_*
python src/b03_phase3_proseg_null.py # coverage/phase3_*
python src/b03_phase4_kernel.py      # margins/phase4_*
python src/b03_phase9_neighbor.py    # sensitivity/breast_tight_bounds_sensitivity.csv
```

Lung crop and second section:

```
python src/b03_lung_audit.py
python src/b03_lung_crop.py
python src/b03_lung_crop_extract.py
python src/b03_lung_load.py
python src/b03_lung_certify.py
python src/b03_lung_null.py
python src/b03_s6_extract.py
python src/b03_s6_prep.py
python src/b03_s6_load.py
python src/b03_s6_certify.py
python src/b03_s6_null.py
python src/b03_s6_finalize.py
python src/b03_phase5_vendor_prep.py
python src/b03_phase5_e2.py
python src/b03_phase5_e3.py
python src/b03_s6_vendor_prep.py
python src/b03_s6_e2.py
python src/b03_s6_tight.py
```

Remaining analyses:

```
python src/b03_phase6_breast.py      # B = 10,000, provenance/b10k_breast
python src/b03_phase6_lung.py        # B = 10,000, provenance/b10k_lung
python src/b03_phase6_analysis.py    # provenance/high_permutation_summary.json
python src/b03_baseline_conventional_null.py
python src/b03_baseline_comparison.py
python src/b03_mass_share.py
python src/b03_synthstudy.py --mode verify
python src/b03_synthstudy_analysis.py
python src/b03_paper_figures.py
python src/b03_schematic.py
```

`src/b03_synthstudy.py --mode run` reproduces the 20 intervention runs, one run
at a time, into `intervention/runs/<run_id>/`. `--mode fleet` runs the batch.

## Reading a section table

`results/breast/scout_final_matched.csv` holds 75 rows, one per pair and region.
The columns record the pair and control flag, the counted transcript totals
`nL` and `nR`, the cell denominators, the departure and arrival bounds `dL_min`,
`dL_max`, `dR_min`, `dR_max`, the point score `T0`, the interval `T_lo`, `T_hi`
and its width, the erosion variants, the matched-null `p_pos`, `p_neg`, `q_pos`
and `q_neg`, the undoubled certificate flags `certified_pos_matched` and
`certified_neg_matched`, the selection-adjusted `p_adj`, `q_adj_pos` and
`q_adj_neg`, and the claim flags `claim_pos` and `claim_neg`.

`results/lung/scout_final_lung.csv` and `results/breast_s6/scout_final_s6.csv`
carry the same schema for 54 and 30 rows. `claim_layer.csv` in each directory
holds the claim computation on its own, with the selected one-sided value and
the undoubled Benjamini–Hochberg values beside the reported ones.

The interval endpoints reproduce from the count columns:

```
T_lo = log2((nL + dL_min)/nl_cells + 0.5) + log2((nR + dR_min)/nr_cells + 0.5)
T_hi = log2((nL + dL_max)/nl_cells + 0.5) + log2((nR + dR_max)/nr_cells + 0.5)
```

`src/b03_eps_sensitivity.py` applies the same expression at ε = 0.1, 0.5 and
1.0 and writes `sensitivity/eps_sensitivity.csv`, which reproduces Table S5.

## What the artifacts show

| Section | Tests | Interval sign (+/−/0) | Claims (+) | Claims (−) |
|---|---|---|---|---|
| Breast Rep1, 313-plex | 75 | 39 / 12 / 24 | 36 | 9 |
| Lung crop, 5K panel | 54 | 1 / 32 / 21 | 1 | 4 |
| Breast, second section, 5K panel | 30 | 0 / 17 / 13 | 0 | 0 |

The joint-attainment test covers 157 of the 159 tests. Two lung tests have no
movable transcript for either gene. The lower endpoint is attained on all 157,
the upper endpoint on 149. `results/coupled_end_summary.json` gives the counts,
the largest shared-loss collision and the largest slack.

The reported interval uses the donor-slot count at Layer A, the layer the matched
null recomputes. A claim additionally requires the null to support the direction
of the interval sign at q <= 0.10 within the region.

Proseg keeps 42 of 45 breast claims inside their intervals (93.3%) and 4 of 5
lung claims (80.0%). The lung fraction falls below the 90% gate fixed in S2, and
the paper reports that outcome.

## Tests

`tests/test_release.py` needs no vendor data. It checks row counts, interval
signs, claim counts, the claim arithmetic, the exact-end summary, the joint
counts, the coverage and margin summaries, and the artifact layout.

```
python -m unittest discover -s tests -v
```

## Data

All three sections are public and are not redistributed here.

| Section | Source |
|---|---|
| Breast Rep1 | 10x Genomics Xenium FFPE Human Breast Cancer Rep1, v1.0.1 |
| Lung crop | 10x Genomics Xenium Prime 5K FFPE Human Lung Cancer, technical note release |
| Breast, second section | EBI zarr repackage of the STHELAR release, BioStudies S-BIAD2146 |

The crops are expression-blind. Their windows are recorded in
`preprocessing/lung_crop_window.json` and
`preprocessing/breast_s6_crop_window.json`.

## Citation

Metadata is in `CITATION.cff`. Code is released under the MIT licence.
