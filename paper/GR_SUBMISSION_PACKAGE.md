# GR Spatial Omics SI — Submission Package (working draft, 2026-09-08)

Submission portal: submit.genome.org. Keyword on submission: "Spatial Omics Special
Issue". Deadline 2026-11-01. Formatting measured 2026-09-11 against the current
manuscript: abstract 264 words (trim to <= 250 at submission); running title 49 chars;
4 figures <= 6; main text ~6.5k words (GR Research Article norms ~4,000 words plus
methods; trim Section 4.1/4.2 methods detail toward the supplement, keep the results
narrative intact).

## Cover letter (draft)

Dear Dr. Sussman and the editors of the Spatial Omics Special Issue,

Most published cell-cell communication claims from imaging-based spatial
transcriptomics rest on a quantity no reader can currently audit: which cell owns
which molecule. We show, on public Xenium sections from two tumor types measured
with panels differing 16-fold, that this inferential step can be turned into a
per-claim certificate: an interval valid over an explicitly declared
uncertainty set of molecule-to-cell assignments, reported only when a null matched
to the worst-case layer confirms the direction.

The biological yield is concrete. On a breast-cancer section, 19 of 60
ligand-receptor contrasts pre-selected as segmentation-fragile collapse to
non-identifiable, including contrasts conventional tests call strong interactions.
On a lung-cancer section, ESR1-PGR mutual exclusion is certified in all three
regions with near-point-identified intervals and replicates fully under an
independent segmentation algorithm, while the same statistic in breast is the
fragile configuration whose cross-segmentation significance fails. One autocrine
certificate (CDH1-EGFR) mirrors breast's CDH1-ERBB2 across tumors, panels, and
instrument generations. The uncertainty set, the pair lists, the null, and the
gates were frozen in pre-registration before any expression was computed; every
number traces to a hashed artifact.

The work matches the special issue's call for computational and statistical
methods for spatial data analysis applied to cell-cell communication in cancer,
and the full pipeline runs in minutes on a workstation. Code, pre-registrations,
and all result artifacts are public (github.com/amywang2010/b03-certificates);
the two datasets are public 10x Genomics Xenium sections.

This manuscript is not under consideration elsewhere and all authors approve
submission. A preprint DOI will be supplied at submission per the journal's
preprint policy.

Sincerely,
Amy Wang, on behalf of all authors

## Statements to paste into the portal

**Data availability.** All data are publicly available: Xenium FFPE Human Breast
Cancer Rep1 (10x Genomics, Xenium v1.0.1; URLs and byte sizes in the ledger,
DEV-013) and Xenium Prime 5K FFPE Human Lung Cancer (10x Genomics technote bundle,
Nov 2024; DEV-013). A third section (Xenium Prime FFPE Human Breast
Cancer, Oct 2024, "breast_s6") is under acquisition for a panel-controlled
replication leg (Phase 7 prereg frozen); if its analysis completes before
submission it is included, otherwise it is reported as ongoing work.

**Code availability.** All code, pre-registrations, deviation ledger, and SHA-256
manifests for every result artifact are available at
https://github.com/amywang2010/b03-certificates. The certificate computation
requires only Python 3.12+ (numpy, pandas, pyarrow); the exact software
environment is recorded in the repository README.

**Author contributions (CRediT).** Conceptualization: A.W. Methodology,
Software, Validation, Formal analysis, Investigation, Data curation: A.W.
Writing - original draft: A.W. Writing - review & editing: A.W. Funding
acquisition: A.W.

**Funding.** [list grants or state "No external funding was used."]

**Competing interests.** The authors declare no competing interests.

**Suggested reviewers (3-5, outside author institutions).** [To select at
submission: candidates are authors of the papers we cite, e.g. in the
segmentation-uncertainty and spatial-methods community, with no collaboration
history; final list confirmed at submission time.]

## Submission checklist

- [ ] bioRxiv DOI minted and inserted (Sep 23-24 per plan)
- [ ] Figures 1-4 uploaded as separate files (600 dpi TIFF/PDF per GR specs)
- [ ] Supplement (B03_SUPPLEMENT.md converted per GR format) with S-tables
- [ ] Phase 6 (B=10k) results integrated before submission
- [ ] Phase 7 (breast_s6) results integrated before submission
- [ ] Cover letter final, keyword set, editor names correct
- [ ] ORCID for every author registered in the portal
- [ ] All artifact manifests regenerated after final artifact freeze
