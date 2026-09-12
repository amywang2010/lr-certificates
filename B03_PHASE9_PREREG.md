# B03 Phase 9 Pre-Registration — Neighbor-Effect Statistic under Certified Bounds (FROZEN 2026-09-08, before any compute)

## Motivation (registered before any compute; no neighbor artifact exists yet)

Every certificate reported so far uses the co-expression statistic T(A,K) of the
locked scope (B03_LOCK): per-region sender/receiver type means of ligand and
receptor. A referee objection applies: a within-region co-expression signal can
arise without any spatial adjacency between sender and receiver cells. Phase 9
adds the adjacency-resolved statistic and shows the same exact-bound machinery
applies to it unchanged.

## Statistic (locked)

For region rho, pair (L,R), sender type S, receiver type R, and adjacency radius
d (registered value: d = 15 um, the Xenium in-cell resolution scale; sensitivity
at d = 25 um reported):

    N_L^adj(rho) = # transcripts of L in rho assigned to cells of type S
                   that have >= 1 cell of type R within distance d of their
                   polygon centroid
    n_S^adj(rho) = # cells of type S in rho with >= 1 type-R neighbor within d

    E_s^adj(g) = N_L^adj / n_S^adj  (analogous definition for R on receivers)

    T^adj(A,K) = log2(E_s^adj(L)) + log2(E_r^adj(R))   (same functional form)

Adjacency is computed ONCE from frozen vendor centroids and frozen labels
(geometry only; independent of the assignment uncertainty set A), so the
neighbor/n denominators are FIXED NUMBERS, not functions of A or K.

## Bound transfer (why the machinery applies unchanged)

Conditional on the frozen adjacency mask M (cell i is neighbor-eligible or not),
the statistic is the same ratio of per-cell-type transcript counts as T(A,K),
restricted to the eligible subset. Every operation in the certified interval
construction (max-flow over the bipartite transcript-to-cell graph, per-gene
leakage intervals, band enumeration) acts on transcript assignment variables;
intersecting with a fixed mask only deletes edges/nodes from the same graph.
Therefore:

  T^adj_lo = min over A in U of T^adj(A,K restricted to U)  is computed by the
  SAME max-flow with the graph restricted to eligible cells; the exactness proof
  transfers verbatim because eligibility is A-independent.

Registered check G9-1: recompute the co-expression bound for one breast pair
with a trivially all-eligible mask and require bitwise equality with the sealed
scout bounds (proves the mask plumbing changed nothing else).

## Testable set (locked)

The 18-pair lung testable set and the 20-pair breast set (B03_LOCK), all
regions. Same preregistered positive/negative control pairs.

## Null and FDR (locked)

The matched worst-case permutation null of the sealed pipeline, applied to the
restricted statistic: for each permutation b, recompute eligibility-restricted
counts under the permuted labels (mask recomputed from the permuted labels --
eligibility follows labels under permutation, exactly as cell counts do), and
compute the same worst-case object. Add-one p-values; BH within region;
q <= 0.10. B = 1,000 for the phase-9 screen; B = 10,000 only for rows that
certify (resolution upgrade registered in advance).

## Gates (frozen before compute)

- G9-1 mask-immutability bitwise check (above). FAIL = quarantine + audit.
- G9-2 alignment: eligible-cell counts recomputed twice (once from centroids,
  once from polygon pairwise distances at d) must agree exactly.
- G9-3 verdict reporting: certified_pos / certified_neg / interval-certified
  only, same definitions as the sealed pipeline; no reinterpretation.
- Pre-specified interpretation rule (registered): if adjacency certification is
  STRICTLY RARER than co-expression certification (e.g. <= half the rows),
  that is reported as a substantive finding (localization burden), not as a
  failure of the method. Both outcomes are publishable rows in the paper.

## Execution plan

Runs AFTER the Phase 6 breast B=10k completes (single-null-per-machine
envelope). New script b03_phase9_neighbor.py; reuses certified-interval
machinery via import from the sealed modules; estimated wall time: minutes
(graph restriction only; one max-flow per pair-region-row as in the sealed
scout, which ran in minutes). Any deviation from this prereg is logged in
B03_DEVIATIONS.md before results are interpreted.

## Amendment A1 (2026-09-08, BEFORE any compute — logged in ledger as DEV-023)

G9-2 as frozen above ("recomputed twice: once from centroids, once from polygon
pairwise distances at d, must agree exactly") is unexecutable as stated: centroid
KD-tree queries and polygon-pairwise distances are different geometries, and no
exact identity between them can hold. The registered check is amended to:
adjacency masks recomputed by per-cell brute-force radius enumeration on the SAME
centroid metric must equal the KD-tree masks exactly, on a 2,000-cell random
sample per (sender-type, receiver-type) key, plus all cells in the
boundary-adjacent band sample. No other section of this prereg is amended.

## Addendum (2026-09-11, FROZEN before corrected compute) — DEV-031 convention correction

The DEV-031 audit (B03_DEVIATIONS.md) established that the sealed breast bounds were
produced by the scout's inline vectorized construction whose reduceat segments absorb
interleaved other-gene donor slots (gains overstated; losses exact); corrected per-gene-slot
bounds reproduce all counts bit-exactly, flip zero verdicts, tighten widths, and newly
interval-certify 3 rows (whose null transfer is registered for Phase 10). This addendum
amends Phase 9 BEFORE its observed-side compute, which had not yet successfully run:

1. **Interval construction.** Phase 9 uses the corrected per-gene-slot construction with
   the DEV-016 capacity formula (the sealed null's own formula), NOT
   `b03_certify.interval_bounds_region` (dead code; id-keyed; cannot execute).
2. **G9-1 (amended).** With an all-eligible mask, the corrected construction must
   reproduce the DEV-031 corrected bounds (results/sensitivity/*_tight_bounds_sensitivity.csv)
   BIT-EXACTLY on every (dL_min, dL_max, dR_min, dR_max), AND the sealed scout_bounds.csv
   intervals must contain the corrected intervals row-wise (superset property — the
   soundness direction check).
3. **G9-0 and G9-2 unchanged** (label-artifact integrity; adjacency determinism with
   explicit self-exclusion, DEV-023).
4. **Outputs.** phase9_neighbor.csv (radius 15 um primary; 25 um registered sensitivity
   pass), columns carry the radius. p-values remain gated on the separate matched null
   (b03_phase9_null.py) under the SAME corrected construction, registered at execution
   time per the single-null-per-machine rule.
5. **No other change.** Eligibility from frozen centroids+labels; sentinel-masked counted
   sets; denominators nl_adj/nr_adj fixed by geometry; exactness-preservation argument
   unchanged (masked restriction of the same bipartite graph).
