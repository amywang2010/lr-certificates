# Certified worst-case bounds for ligand–receptor statistics over segmentation uncertainty

Code, pre-registrations, and hashed result artifacts for the manuscript submitted to the
Genome Research Spatial Omics Special Issue.

## What the paper claims

Every cell–cell communication statistic from imaging-based spatial transcriptomics
inherits the ambiguity of cell segmentation: which cell owns which molecule. This work
converts that ambiguity into a per-claim certificate: an interval valid over an explicitly
declared uncertainty set of molecule-to-cell assignments (boundary-local reassignment,
mask erosion/dilation, per-gene leakage kernels), reported only when a permutation null
matched to the worst-case layer confirms the direction. On public Xenium sections from two
tumor types (breast, Rep1 313-plex panel; lung, Prime 5K panel), 19 of 60 fragile breast
contrasts collapse to non-identifiable while ESR1–PGR mutual exclusion certifies in all
lung regions and replicates under an independent segmentation algorithm.

## Repository map

| Path | Content |
|---|---|
| `paper/B03_MANUSCRIPT.md` | Manuscript (main text) |
| `paper/B03_SUPPLEMENT.md` | Supplement (datasets, gates, verification receipts, S8.1 independent re-derivation) |
| `paper/figures/` | Figures 1–4 (PNG 600 dpi + PDF), generated only from sealed artifacts |
| `B03_PREREGISTRATION.md`, `B03_PHASE*_PREREG.md` | Frozen pre-registrations (one per phase, each dated before its compute) |
| `B03_PHASE*_REPORT.md`, pre-registration addenda | Per-phase reports and dated amendment records (deposited with the journal submission) |
| `src/` | Pipeline: load/concordance (`b03_load.py`, `b03_scout.py`), certified bounds (`b03_scout.py` inline construction), matched null (`b03_robust_null_v3.py`), cross-segmentation (`b03_lung_*.py`, `b03_phase5_*`), kernel margins (`b03_phase4_*`), neighbor effect (`b03_phase9_neighbor.py`) |
| `results/` | Hashed artifacts: `scout_bounds.csv`, `null_robust_matched.csv`, `b10k_breast/` (B=10,000 null + checkpoint), `exactness_verification.csv`, `sensitivity/` (independent re-derivation), phase reports |
| `B03_THEORY.md` | Bound construction and exactness argument |

## Provenance chain

Every number in the manuscript traces to a SHA-256-hashed artifact (`results/*.sha256`
manifests per phase). Raw data are public 10x Genomics Xenium downloads (breast Rep1,
lung Prime 5K); loaders reproduce the vendor assignment from molecules + polygons and
verify against the vendor cell-by-gene matrix before any statistic is computed
(`B03_DATA_AUDIT.md`, gates G1–G3). The interval construction was independently
re-derived after the initial analysis; the re-derivation reproduces every certified
verdict bit-exactly and shows the reported intervals are conservative supersets
(`results/sensitivity/`, Supplement S8.1).

## Reproduction

Each phase script runs standalone on a workstation (16 GB RAM class); the permutation
nulls checkpoint every 50 permutations and resume exactly. Expected runtimes: certified
bounds ~2 min; B=10,000 matched null ~12 h on 4 workers (checkpointed, resumable);
cross-segmentation battery < 1 h.
