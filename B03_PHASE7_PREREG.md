# B03 Phase 7 Pre-Registration — breast_s6 leg (FROZEN 2026-09-08, before any download completes)

## Motivation (registered)

The paper's central cross-tissue claim (ESR1↔PGR: fragile partial negative in
breast-313plex vs near-point-identified exclusion in lung-5K) is confounded between
tissue and panel. One public section unconfounds both axes at once:

- **breast_s6** (10x "Xenium Prime FFPE Human Breast Cancer", Oct 24 2024, XOA v3.0.0,
  5K Pan Tissue & Pathways + 100 custom genes): same tissue as tissue 1, same panel
  family as tissue 2.

Two registered comparisons, decided before data access:
1. breast_s6 vs breast Rep1: same tissue, panel 16x -> isolates the PANEL axis.
2. breast_s6 vs lung crop: same panel family, different tissue -> isolates the
   TISSUE axis for the ESR1↔PGR asymmetry.

Both outcomes are informative and reportable; neither is gated as success/failure
of the framework. The claim in the manuscript will be stated as whatever the data
show, with the asymmetry attributed accordingly.

## Section provenance

Source: the ORIGINAL 10x dataset bundle (transcripts.parquet, cell_boundaries.parquet,
cell_feature_matrix.h5, gene_panel.json, metrics_summary.csv), not the STHHELAR
zarr repackage (S-BIAD2146 documents its provenance: 10x CC-BY-4.0). Download
verified by HTTP 200 + exact byte size per file (DEV-013 discipline), URLs recorded
in the ledger before extraction.

## Gates (unchanged from the frozen template)

- **G0 (panel, pre-expression):** the 18-pair lung testable set requires 17 unique
  genes; all 17 verified present in the lung 5K panel file (2026-09-08, local).
  breast_s6 G1: parse breast_s6's own gene_panel.json inside the zip; the tested
  set is the mechanical subset with both genes present; exclusions recorded with
  missing genes. If any of the 17 is missing, the affected pairs are dropped
  mechanically and the comparison set shrinks (registered, not fatal).
- **G3 (assignment layer):** PIP reproduction of the vendor matrix on the full
  section; per-cell concordance >= 0.90, total ratio >= 0.97.
- **Resource envelope before any proseg launch:** transcripts x 230 B vs 15.7 GB
  physical; if the full section exceeds the envelope, the expression-blind crop
  rule (DEV-017 recipe: density-median window over a 0.5 mm grid, zero expression
  access, envelope gate <= 7.0 GB predicted) applies unchanged.
- **Labels:** the frozen 6-class tissue-generic dictionary (Epithelial, Myeloid,
  Fibroblast, Endothelial, T_cells, B_cells) with the assertion-guarded tested-gene
  leak filter; >= 200 cells per class floor; classes below the floor are dropped
  from the transfer and recorded.
- **Regions:** k-means k=3 on (log10(1+counts), x, y), frozen once from A0.
- **Certification:** identical machinery (interval + matched null). Matched null at
  **B = 10,000** (Phase 6 convention; per-perm seeds registered as 20260908 +
  100003 + p), BH within region q <= 0.10. DEV-016 cap formula active from the start.
- **E2/E3/E4:** Proseg on the same (full or crop) molecule set; gates 90% / 85% /
  descriptive, interpretation clause identical (excursions reported with mechanism).

## Outputs

`xenium_breast_s6/` tree mirroring xenium_lung/crop; scout_bounds, scout_final,
matched-null CSVs, coverage/margins summaries; one report (B03_PHASE7_REPORT.md);
supplement S-table row set; manuscript Section 6.7 (one page) + limitations update
(7.4.1 single-section limitation addressed; 7.4.4 third-segmenter remains for the
multimodal leg, tracked as Phase 8).

## Kill / descope clauses (registered)

- Dataset unavailable or hash-unverifiable: phase aborts, no partial claims.
- G3 concordance below gate: abort before any biology.
- Any tested-gene leak in the dictionary: mechanical filter + DEV entry; abort if
  any class loses its floor as a result.

## Estimated effort

Download ~20-40 GB overnight; G1+G3 ~30 min; proseg ~1-2 h if within envelope;
null at B=10k ~2-6 h depending on crop size; E2-E4 < 30 min. Total: ~1 day of
machine time, overlapped with Phase 6's breast null (this phase starts after the
Phase 6 chain completes; the machine runs one null at a time by envelope design).
