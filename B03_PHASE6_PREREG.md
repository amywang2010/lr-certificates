# B03 Phase 6 Pre-Registration — B=10,000 null resolution (FROZEN 2026-09-08, before any compute)

## Motivation (registered, not post hoc)

Manuscript limitation 7.4.5: at B=1000 the permutation p-value floors at
9.99e-4 (add-one convention, Phipson & Smyth 2010); several headline certified
rows sit at the floor. This phase sharpens resolution 10-fold for every row in
both tissues. No verdict-generating rule changes.

## Frozen parameters

- B = 10,000 permutations per tissue; region- and density-preserving label
  permutations, identical scheme to Phase 1/5 (same code path, `b03_robust_null_v3.py`
  and its crop-parameterized lung driver `b03_lung_null.py`).
- Seed stream UNCHANGED: perm p uses `default_rng(20260905 + 100003 + p)`. Perms
  0..999 are therefore bitwise identical to the B=1000 runs. This yields a
  registered determinism check: the per-permutation T statistics for perms 0..999
  must match the Phase 1/5 checkpoint rows exactly (max |delta| = 0). Failure of
  this check kills the phase and triggers a code audit before any output is used.
- Workers: 4 (breast), 4 (lung); checkpointing every 50 perms to `*.tmp.npz`;
  fresh output directories (`results/b10k_breast/`, `xenium_lung/crop/data/b10k/`)
  so the sealed B=1000 artifacts remain untouched.
- Uncertainty set, capacity caps (DEV-016 formula), pair lists, regions, labels,
  epsilon = 0.5: unchanged from the frozen per-tissue preregistrations.
- p = (1 + #{T_b <=/< T_obs}) / (B + 1); BH within region at q <= 0.10; certified
  direction rules unchanged. The certification verdicts are expected to be
  IDENTICAL to the B=1000 runs (all current verdicts have q <= 0.10 with large
  margin); a verdict flip would be a registered anomaly to be investigated, not
  silently absorbed.

## Outputs

`b10k_null.csv` per tissue with per-row p/q at B=10,000; determinism check JSON;
one report section appended to the supplement (S-table) and one sentence each in
Abstract/Results where the p-floor was disclosed. Manuscript limitation 7.4.5
replaced by the sharpened result.

## Runtime estimate (resource audit)

Lung: B=1000 took 417 s on the crop -> ~70 min at B=10,000. Breast: B=1000 took
74 min -> ~12 h. Sequential chain (lung first), 4 workers each, 8-core machine:
total ~13 h, completing overnight. Disk 176 GB free; memory envelope within the
proven Phase 1/5 footprint (no new allocations; B enters loop count only).
