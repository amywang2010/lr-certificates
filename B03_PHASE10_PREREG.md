# B03 Phase 10 Pre-Registration — corrected-rule primary numbers (FROZEN 2026-09-11, before any compute)

## Registered motivation (documentation-level, from sealed artifacts)

DEV-031 established that the sealed interval construction overstated gain ends (loss ends
exact) on both tissues; verdicts and p-values are invariant, and the corrected per-gene-slot
rule certifies strictly more (3 additional interval-negative rows). To make every primary
number one-rule-consistent end to end, the matched nulls and all derived artifacts are
re-derived under the corrected rule. Nothing here changes any registered *definition*
(estimand, uncertainty set, permutation scheme, FDR); only the interval-construction
implementation inside the null changes, to the DEV-031 corrected rule.

## What runs (in order, sequential; single-null-per-machine rule)

1. **P10-1 breast corrected null (b03_robust_null_v4.py, B=10,000).**
   - Identical permutation machinery to the sealed v3 (blocks = region × density decile;
     seed 20260906; add-one p convention p=(1+#{...})/(B+1); BH within region over the pair
     family, q<=0.10 thresholds; parent + 2 spawn workers; pre-flight determinism perm-0
     bitwise across parent/workers; 10 named permutation-invariance checks).
   - THE ONE CHANGE: interval_obs and the per-permutation interval builder use per-gene
     own-slot aggregation (DEV-031 corrected rule) with the DEV-016 capacity formula
     (never binds on breast; asserted).
   - Observed-interval cross-check (startup assertion): the corrected builder's all-eligible
     intervals must reproduce results/sensitivity/breast_tight_bounds_sensitivity.csv
     BIT-EXACTLY (75/75), superseding the sealed scout_bounds comparison of v3.
   - Output: results/b10k_breast_v4/ with b10k_null.csv (same columns), checkpoints at
     every 50 permutations (B03_P10_CKPT override allowed; default 50), ledger JSON.
2. **P10-2 lung corrected null (same script, lung paths, B=10,000).** DEV-016 caps bind
   (197 cells) and are applied identically in observed and permuted intervals (one rule).
   Output: xenium_lung/crop/data/b10k_v4/.
3. **P10-3 verdict-stability vs sealed.** Registered reporting, not a gate: certified sets
   at B=10k under v4 vs the sealed v3 certified sets; expected by DEV-031: zero losses,
   possible gains (MS4A1->CD274 r0/r1 breast; CDH1->EGFR r1 lung). Any LOSS would be a
   registered anomaly (halt, human review — it would contradict DEV-031's monotonicity
   result).
4. **P10-4 artifact refresh.** scout_final_matched.csv (breast), scout_final_lung.csv,
   phase4_kernel_margins.csv (both tissues), paper_fig1/2/3 data CSVs, coverage tables:
   re-derived from v4 p/q and corrected intervals. Figures regenerated; every changed
   number in the manuscript gets its new value with the DEV-031 citation.

## Gates (all pre-registered, mechanical)

- G10-1: v4 observed-interval cross-check bitwise vs the DEV-031 corrected CSVs (75 + 54).
- G10-2: perm-0 bitwise determinism across parent + workers (both tissues).
- G10-3: checkpoint continuity — resumed runs must reproduce checkpoint rows exactly.
- G10-4 (registered expectation, halt-on-violation): v4 certified_neg set ⊇ sealed
  certified_neg set, and v4 certified_pos set ⊇ sealed certified_pos set (monotonicity;
  losses prohibited by DEV-031's verified direction).

## Interpretation (frozen)

- Primary numbers for the manuscript are the v4 outputs. Sealed v3 artifacts remain in the
  repository untouched as the audit trail (DEV-031 documents their role).
- The 3 newly interval-certified rows gain claim certification iff v4 q_neg <= 0.10 for
  them (mechanically decided by P10-1/P10-2; no discretion).
- No other interpretation change of any kind. BioRxiv/GR text cites DEV-031 + this prereg
  for the one-rule consistency of all primary numbers.

## Resource plan

- Sequential breast-then-lung (memory envelope identical to sealed runs: parent + 2 workers).
- Launch tonight after breast_s6 assembly completes or immediately, whichever is later in
  wall-clock; expected ~10-14 h breast, ~2-3 h lung on the registered machine.
