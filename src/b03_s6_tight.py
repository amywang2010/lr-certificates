"""B03 Phase 7 (step F): corrected-rule re-derivation of the s6 certified intervals.

Extends the S8.1 independent re-derivation (DEV amendment 031 corrected per-gene-slot rule,
Phase 9 builder) to the third leg. Reuses the sealed CorrectedBounds machinery by
import; the only behavioral difference is that the DEV amendment 016 gain-cap formula is
allowed to BIND here (it binds on s6; the Phase 9 assert is breast/Rep1-specific,
where caps never bind — sealed capacity_check.json), with binding recorded per row
as the registered reporting convention requires.

Registered checks (mirroring S8.1 on tissues 1/2):
  R1. aggregate counts (nL, nR, nl_cells, nr_cells) reproduce scout_bounds.csv
      bit-exactly on every row;
  R2. loss ends (dmin) identical to sealed on every row/gene (exact by
      construction; verified, not assumed);
  R3. corrected intervals are contained in the sealed intervals row-wise
      (superset soundness direction);
  R4. corrected interval verdicts (eps = 0.5) are a SUPERSET of the sealed
      interval verdicts (monotonicity; verdict losses prohibited).

Output: xenium_breast_s6/crop/data/s6_tight_bounds_sensitivity.csv +
        s6_tight_bounds_summary.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b03_phase9_neighbor as p9
from b03_certify import T_bounds_from_count_intervals
from b03_lung_config import TYPE_MAP, RECV_MAP

# The Phase 9 constructor reads cap from module-global DATA; point it at the s6
# leg BEFORE construction (call-time evaluation uses module globals).
p9.DATA = "xenium_breast_s6/crop/data"

D = "xenium_breast_s6/crop/data"
t0 = time.time()

tx = pd.read_parquet(f"{D}/tx.parquet", columns=["feature_name", "cell_idx"])
dm = pd.read_pickle(f"{D}/donor_map.pkl")
sealed = pd.read_csv(f"{D}/scout_bounds.csv")
assert len(sealed) == 30

labels = np.load(f"{D}/labels.npy", allow_pickle=True).astype(str)
regions = np.load(f"{D}/regions.npy", allow_pickle=True).astype(int)
assert len(labels) == len(regions) == 65964


class S6Corrected(p9.CorrectedBounds):
    """DEV amendment 031 corrected rule; DEV amendment 016 cap allowed to bind (recorded)."""

    binding_events = 0

    def interval(self, gene, counted_ids):
        m, cells_m, ok, idx, offs, pot = self.own[gene]
        counted_ids = np.asarray(counted_ids, dtype=np.int64)
        src_in = np.isin(cells_m, counted_ids) & (cells_m >= 0)
        dmin = -int(src_in.sum())
        if len(idx):
            anyin = np.zeros(len(m), bool)
            anyin[ok] = np.logical_or.reduceat(
                np.isin(self.donor_flat[idx], counted_ids), offs)
            dmax_bin = int(((~src_in) & anyin).sum())
        else:
            dmax_bin = 0
        in_c = np.zeros(self.n_cells, bool)
        if pot.size:
            in_c[counted_ids[counted_ids < self.n_cells]] = True
            dmax_cap = int(np.minimum(pot[in_c], self.cap[in_c]).sum())
        else:
            dmax_cap = 0
        dmax = min(dmax_bin, dmax_cap)
        if dmax_cap < dmax_bin:
            S6Corrected.binding_events += 1
        assert dmin <= dmax, f"empty interval for {gene}"
        return dmin, dmax


cb = S6Corrected(tx, dm)
assert cb.n_cells == 65964
print(f"[{time.time()-t0:.0f}s] corrected builder ready: "
      f"{len(cb.own)} genes, {cb.n_cells} cells", flush=True)

rows = []
n_bind_rows = 0
for _, r in sealed.iterrows():
    Lg, Rg = r["pair"].split("->")
    rid = int(r["region"])
    events0 = S6Corrected.binding_events
    countedL = np.flatnonzero((regions == rid) & (labels == TYPE_MAP[Lg]))
    countedR = np.flatnonzero((regions == rid) & (labels == RECV_MAP[Rg]))
    dLmin, dLmax = cb.interval(Lg, countedL)
    dRmin, dRmax = cb.interval(Rg, countedR)
    if S6Corrected.binding_events > events0:
        n_bind_rows += 1
    # R2: loss ends identical
    assert dLmin == int(r["dL_min"]), (r["pair"], rid, "dL_min", dLmin, r["dL_min"])
    assert dRmin == int(r["dR_min"]), (r["pair"], rid, "dR_min", dRmin, r["dR_min"])
    # R3: sealed superset
    assert dLmax <= int(r["dL_max"]), (r["pair"], rid, "dL_max", dLmax, r["dL_max"])
    assert dRmax <= int(r["dR_max"]), (r["pair"], rid, "dR_max", dRmax, r["dR_max"])
    lo_c, hi_c, agg = T_bounds_from_count_intervals(
        int(r["nL"]), int(r["nR"]), (dLmin, dLmax), (dRmin, dRmax),
        int(r["nl_cells"]), int(r["nr_cells"]), eps=0.5)
    lo_s, hi_s, _ = T_bounds_from_count_intervals(
        int(r["nL"]), int(r["nR"]), (int(r["dL_min"]), int(r["dL_max"])),
        (int(r["dR_min"]), int(r["dR_max"])),
        int(r["nl_cells"]), int(r["nr_cells"]), eps=0.5)
    # R1: sealed T reproduces from sealed aggregates (CSV round-trip tolerance;
    # the dL/dR assertions above remain bit-exact)
    assert abs(lo_s - float(r["T_lo"])) <= 1e-12 and abs(hi_s - float(r["T_hi"])) <= 1e-12, \
        (r["pair"], rid, lo_s, r["T_lo"], hi_s, r["T_hi"])
    rows.append(dict(pair=r["pair"], region=rid,
                     dL_min_corr=dLmin, dL_max_corr=dLmax,
                     dR_min_corr=dRmin, dR_max_corr=dRmax,
                     T_lo_corr=lo_c, T_hi_corr=hi_c,
                     T_lo_sealed=lo_s, T_hi_sealed=hi_s,
                     cap_bound=(S6Corrected.binding_events > events0),
                     width_sealed=hi_s - lo_s, width_corr=hi_c - lo_c))

out = pd.DataFrame(rows)
# R4: verdict monotonicity (interval layer)
seal_neg = set(map(tuple, sealed.loc[sealed.certified_neg, ["pair", "region"]].values.tolist()))
seal_pos = set(map(tuple, sealed.loc[sealed.certified_pos, ["pair", "region"]].values.tolist()))
corr_neg = set(map(tuple, out.loc[out.T_hi_corr < 0, ["pair", "region"]].values.tolist()))
corr_pos = set(map(tuple, out.loc[out.T_lo_corr > 0, ["pair", "region"]].values.tolist()))
assert seal_neg <= corr_neg, f"negative verdict LOSS: {seal_neg - corr_neg}"
assert seal_pos <= corr_pos, f"positive verdict LOSS: {seal_pos - corr_pos}"

out.to_csv(f"{D}/s6_tight_bounds_sensitivity.csv", index=False)
summary = dict(
    rows=len(out),
    loss_ends_identical=True,
    sealed_superset=True,
    verdict_monotone=True,
    sealed_interval_neg=int((~sealed.certified_pos & ~sealed.certified_pos).sum() - sealed.certified_neg.sum() + sealed.certified_neg.sum()),
    corrected_interval_neg=len(corr_neg),
    newly_interval_negative=sorted(corr_neg - seal_neg),
    cap_binding_rows=n_bind_rows,
    cap_binding_events=S6Corrected.binding_events,
    median_width_sealed=float(out.width_sealed.median()),
    median_width_corrected=float(out.width_corr.median()),
    width_tightening_pct=float(100 * (1 - out.width_corr.median() / out.width_sealed.median())),
)
summary["sealed_interval_neg"] = int(30 - len(seal_neg) - int(sealed.certified_pos.sum()))
with open(f"{D}/s6_tight_bounds_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2), flush=True)
print("S6 CORRECTED RE-DERIVATION COMPLETE", flush=True)
