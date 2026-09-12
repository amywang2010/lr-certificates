# B03 Phase 8 Pre-Registration — third assignment family (FROZEN 2026-09-08)

## Registered finding (documentation-level, from sealed artifacts)

The metrics_summary.csv of the sealed lung leg records
`segmented_cell_stain_frac = 0.978`: the lung vendor assignment was produced by
the XOA v3.0.0 cell-segmentation-staining (multimodal) workflow. The breast
Rep1 vendor assignment (XOA v1.0.1) is nucleus-expansion only. The existing
design therefore already spans two vendor assignment families, and breast_s6
(XOA v3.0.0) will provide a multimodal-base + Proseg-probe comparison as part
of Phase 7.

## Remaining gap this phase closes

A re-segmentation engine that is neither the vendor pipeline nor Proseg, applied
to a section we have certified. Two candidate sources, in priority order:

1. **STHELAR cell-boundary inventory (S-BIAD2146).** The STHHELAR release ships
   `cell_boundaries` polygons per slide (10x-origin, CC BY 4.0). For any slide
   whose tissue appears in our certified set (breast, lung), the shipped
   boundaries define a third assignment source. Registration requires confirming,
   before use, that the shipped boundaries differ from the corresponding 10x
   bundle boundaries (if identical, they are not a third family and the phase
   reports that null result).
2. **A public multimodal-segmentation Xenium section carrying a pair gene set
   that overlaps our tested pairs** (e.g. pancreas multimodal, XOA v2.0.0),
   analyzed end-to-end with the crop rule and the standard battery.

## Procedure (per source)

- Pre-register the specific slide/section and the exact assignment-source file
  before computing anything.
- Compute certified intervals once per (pair, region) under the third assignment
  as the base A0 (labels via the same frozen transfer protocol; regions frozen
  from that A0); report coverage of the ORIGINAL certificates' intervals and the
  fraction of certified directions preserved.
- Gates: descriptive-first. If coverage falls below the E2 gate (90%), the
  finding is reported as a calibration result (the uncertainty set's empirical
  coverage under a third family), not a failure of the framework; the
  pre-registered interpretation clause applies (mechanism named, no silent drop).
- All deviations DEV-020+.

## Outputs

One supplement section ("sensitivity to the segmentation algorithm family"),
one manuscript sentence in 7.4.4 replacing "Proseg is the only re-segmentation
engine we ran" with the measured third-family result, and one figure-panel
candidate if the result is informative.

## Effort

1-2 days machine-light work after Phase 7 lands; entirely public data.
