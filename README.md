# Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty

Code and SHA-256-hashed result artifacts for the study:

> **Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty in spatial transcriptomics**

The study computes, for each ligand–receptor (LR) contrast in a spatial transcriptomics
section, an **interval valid over an explicit uncertainty set of molecule-to-cell
assignments** (boundary-local reassignment with capacity constraints, mask
erosion/dilation, per-gene leakage intervals). A biological direction is reported only
when both bounds exclude zero **and** a permutation null matched to the worst-case layer
of the statistic confirms it (within-region Benjamini–Hochberg FDR). Everything else is
reported as non-identifiable. This repository holds the reproduction materials only: the analysis code, the
frozen run configuration, and the result artifacts every number traces to.

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

A controlled 20-run intervention study on the first section (paper Section 5.6) tests
this mechanism directly: injected-pair certification transitions between 1% and 4%
ambient tested-gene mass share (seed-replicable), and certificate soundness under
gene-name relabeling is witnessed bitwise (87/87 class-preserving slots match the
real section's counts exactly). Per-run artifacts and the fleet receipt are under
`results/synthstudy/`; `src/b03_synthstudy_analysis.py` recomputes every Section 5.6
number from them.

The conventional co-expression permutation pipeline (the CellChat/CellPhoneDB-family
null), run on the same breast section under the registered permutation scheme, calls
73 of 75 contrasts interactions; the certificates confirm 39, render 28
non-identifiable, and reverse the six ESR1-PGR rows the conventional test calls
positive (`results/baseline_comparison/`, supplement S8.2).

Cross-segmentation validation uses Proseg as an independent probe on all three sections:
coverage of certified intervals (E2), cross-segmentation significance (E3), and
per-certificate kernel margins (E4). All analysis-free design choices (regions, crop
windows, gates, thresholds) were fixed before any expression computation and are stated
in the paper and this README.

## Repository map

| Path | Content |
|---|---|
| `results/synthstudy/` | Controlled intervention study (paper Section 5.6): one directory per run with `run_manifest.json` (arm, seed, measured injected mass, verdict counts), `scout_bounds.csv`, `null_robust_matched.csv`, `scout_final_matched.csv` (75 rows each), plus the fleet receipt `fleet_summary.csv`; 20 runs = 6 name-permutation replicates (N01..N06), 6 base-null + ring-injection runs (E*), 8 thinning + injection runs (M*) |
| `src/` | Analysis code, one script per step: load and concordance (`b03_load.py`, `b03_scout.py`), certified bounds (inline in `b03_scout.py`), matched null (`b03_robust_null_v3.py`), B=10,000 chain (`b03_phase6_*.py`), cross-segmentation and kernel margins (`b03_phase2..5_*.py`), lung pipeline (`b03_lung_*.py`), second-breast-section pipeline (`b03_s6_*.py`), neighbor effect (`b03_phase9_neighbor.py`), figures (`b03_paper_figures.py`), and the
intervention study (`b03_synthstudy.py`, `b03_synthstudy_analysis.py`) |
| `results/` | Breast-section artifacts: `scout_final_matched.csv` (75 rows, the source of Figures 1 and 3), `null_robust_matched*.csv`, `b10k_breast/` and `b10k_lung/` (B=10,000 reruns), `phase6/` (determinism gate), `exactness_verification.csv`, `sensitivity/` (independent re-derivation of all intervals), per-directory manifests `*.sha256`; `baseline_comparison/` (conventional co-expression permutation baseline vs certificates, the source of the paper's Table 1); `results/synthstudy/` holds the 20 per-run directories of the intervention study plus `fleet_summary.csv` (two B=10 smoke checks retained alongside the 20 fleet runs) |
| `lung/` | Lung-crop artifacts, same schema: `scout_final_lung.csv` (54 rows), null outputs, E2/E3/E4 summaries |
| `breast_s6/` | Second-breast-section artifacts, same schema: `scout_final_s6.csv` (30 rows), matched-null output at B=10,000 (`null_robust_matched_final.csv`), E2 coverage (`phase7_e2_coverage.csv`, `phase7_e2_summary.json`), corrected-rule re-derivation (`s6_tight_bounds_sensitivity.csv`, `s6_tight_bounds_summary.json`), verification receipts (`exactness_verification.csv`, `subcube_check.json`, `capacity_check.json`, `type_absent_skips.json`, `concordance.json`) |
| `RELEASE_MANIFEST.sha256` | Hashes this repository exactly as shipped; verifies byte-for-byte on any platform |

## Provenance

Every number in the manuscript traces to a SHA-256-hashed artifact in `results/`,
`lung/`, or `breast_s6/` (per-directory manifests record artifact and source hashes as they
existed at analysis time). Loaders reproduce each vendor assignment from molecules and
polygons and verify it before any statistic is computed: per-gene bitwise or
totals-level gates per package (breast Rep1 exact-per-cell ≥ 0.90 with per-gene r ≥ 0.958
over 40 audited genes; lung 20/20 audited genes bitwise-exact; second breast section
totals-level gates, r = 0.999, per the packaging addendum). The interval construction was independently re-derived after the
initial analysis; the re-derivation reproduces every certified verdict and shows the
reported intervals are conservative supersets (`results/sensitivity/`, and for the third
section `breast_s6/s6_tight_bounds_sensitivity.csv`).

## Reproduction

All scripts resolve paths from the repository root (independently of the working
directory). Breast-section scripts read the vendor bundle from `data/` and write
intermediates to `data/` and `results/`. The lung and second-breast pipelines operate on
a small workspace directory next to the repository clone (its parent directory by
default; override with the `B03_WORKSPACE` environment variable): the lung workspace is
`xenium_lung/` (`Xenium_Prime_Human_Lung_Cancer_FFPE_outs.zip`, an `extracted/`
subdirectory, and a `crop/` subdirectory), the second-breast workspace is
`xenium_breast_s6/` (the S-BIAD2146 zarr zip, `extracted/`, `crop/`, `chunks/`). Vendor
bundle sizes make these directories unsuitable for version control; everything they
produce that a number depends on is released under `results/`, `lung/`, and
`breast_s6/`.

Hardware used: 8-core workstation, 15.7 GB RAM, Windows. Full-section Proseg on the lung
and 5,000-plex datasets is memory-infeasible at this scale; both use an expression-blind
crop: 3.5 x 3.5 mm candidate windows on a 0.5 mm grid over the full section bounds,
ranked by vendor-cell nucleus count with zero access to expression data, selecting the
median-density window (ties broken lexicographically). The frozen windows are
`results/lung_crop_window.json` and the crop record in `breast_s6/`.

1. Download the three public datasets (sources under Data availability; the lung
   bundle is 38,795,733,696 bytes, SHA-256 verified per chunk at download) and place
   the breast bundle files in `data/` and the other two bundles in the workspace
   directories described above.
2. Breast pipeline: `src/b03_load.py` → `src/b03_scout.py` → `src/b03_robust_null_v3.py`
   → `src/b03_finish.py` (writes `data/tx.parquet`, `data/cells_meta.parquet`,
   `data/donor_map.pkl`, `data/labels.npy`, `data/regions.npy`, and the panel TSV
   alongside the vendor files in `data/`; also copy
   `Xenium_FFPE_Human_Breast_Cancer_Rep1_panel.tsv` into `data/`).
3. Cross-segmentation: install Proseg 3.2.0 and run `proseg transcripts.parquet` on the
   identical molecule file; then
   `src/b03_phase2_vendor_prep.py`, `src/b03_phase2_analysis.py`,
   `src/b03_phase3_proseg_null.py`, `src/b03_phase4_kernel.py`.
4. Lung pipeline: `src/b03_lung_audit.py` → `src/b03_lung_crop.py` →
   `src/b03_lung_crop_extract.py` → `src/b03_lung_load.py` → `src/b03_lung_certify.py`
   → `src/b03_lung_null.py`, then the E2/E3/E4 scripts under `src/b03_phase5_*.py`.
5. Second breast section (5,000-plex): `src/b03_s6_extract.py` → `src/b03_s6_prep.py`
   (frozen expression-blind crop) → `src/b03_s6_load.py` → `src/b03_s6_certify.py` →
   `src/b03_s6_null.py` → `src/b03_s6_finalize.py` → `src/b03_s6_vendor_prep.py` →
   `src/b03_s6_e2.py` → `src/b03_s6_tight.py` (the corrected-rule re-derivation). Proseg
   input construction is in `src/b03_s6_vendor_prep.py`.
6. Controlled intervention study (paper Section 5.6): requires the breast-section
   intermediate artifacts from step 2 (`data/tx.parquet`, `cells_meta.parquet`,
   `donor_map.pkl`, `labels.npy`, and the panel TSV). Then `python
   src/b03_synthstudy.py --mode verify` (invariants), `--mode run` for each of the
   20 runs (arm, seed, and measured injection recorded in each `run_manifest.json`),
   and `--mode analyze` (writes `fleet_summary.csv`). Fleet wall time ~2.5 h on the
   reference workstation; results are deterministic given the recorded seeds.
7. Figures: `python src/b03_paper_figures.py` (reads only released artifacts from
   `results/` and `lung/`; assertions verify row counts before plotting; writes PNG and
   PDF into `figures/`).

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
repackage (accession S-BIAD2146); its packaging differences from the 10x bundles
(zarr layout, deprecated codeword categories) are handled explicitly in
`src/b03_s6_extract.py` and `src/b03_s6_load.py`. The lung and second-breast analysis
regions are expression-blind spatial crops; boundaries are the window records above,
and the selection rule is stated under Reproduction.
