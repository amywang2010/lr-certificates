# B03 Phase 2 — Proseg Cross-Platform Replication Report

## Verdict: **PASS**

Pre-registered gate (B03_PHASE2_PREREG.md, Amendment A3 semantics): **72/75** evaluable rows inside their certified intervals = **96.0%** (gate >= 90%; 0 rows excluded as non-evaluable, >10% exclusion aborts).

## Configuration provenance

- Proseg 3.2.0, run 4 per Amendment A4: 2 um voxels, 2 z-layers (memory-feasible configuration after the run-3 guard abort, DEV-009); diffusion model KEPT; all other parameters default.
- Exact command: `proseg B03_project/data/Xenium_FFPE_Human_Breast_Cancer_Rep1_transcripts.csv.gz --xenium --burnin-voxel-size 4 --voxel-size 2 --voxel-layers 2 --output-counts counts --output-cell-metadata cells --output-gene-metadata genes --output-cell-metadata-fmt parquet --output-gene-metadata-fmt parquet --output-path B03_project/data/proseg_out/`
- Cells: 165,094 (background rows dropped per A2a: 0); counts as gzip MatrixMarket, rounded to integer point estimates (Proseg 3 default).
- Gene-alignment proof (A2c): 0 invariant violations across needed genes; noise-slot handling: n_named layout (none dropped).

## Endpoint arithmetic

- Rows inside certified intervals: **72/75** (96.0%)
- Certified rows (matched-null scout, q <= 0.10): 37 positive, 9 negative; sign agreement under Proseg: **46/46**
- Label transfer (frozen vendor centroids, cosine, tested-LR genes excluded): 2-fold CV accuracy **72.4%** on 166,157 vendor cells (15 types); analysis-critical types (endothelial, macrophage, T, fibroblast, myoepithelial) 85-97% (confusion matrix: phase2_label_transfer_cv_confusion.csv).

## Rows moved outside certification by Proseg (largest first)

| pair | region | vendor T0 | certified interval | Proseg T | distance outside |
|---|---|---|---|---|---|
| ERBB2->PDCD1 | 0 | +3.21 | [+1.39, +3.61] | +3.66 | 0.05 |
| ESR1->PGR | 0 | -0.71 | [-1.51, -0.27] | -0.25 | 0.03 |
| PGR->ESR1 | 0 | -0.71 | [-1.51, -0.27] | -0.25 | 0.03 |


## Interpretation contract (frozen in prereg)

- PASS: certification coverage transfers to an independent probabilistic segmentation algorithm — the uncertainty-set thesis gains cross-platform support.
- PARTIAL: coverage is resolution-dependent; mechanistic diagnosis of which pairs/regions move is REQUIRED before any claim.
- FAIL: U does not transfer across segmentation algorithms — a major, publishable negative result, reported openly.
