"""B03 Phase 4: leakage-kernel sensitivity — certification margins per
B03_PHASE4_PREREG.md (FROZEN before compute; amendment record).

Model: scalar kernel multiplier c >= 1 scales the scout's geometric worst-case
flows: effective count interval per gene
    N_G in [max(0, nG + c*dG_min), nG + c*dG_max]
c = 1 reproduces the scout kernel exactly (Assertion A validates this against
scout_bounds.csv on all 75 rows).

Estimand: per certified row, the critical multiplier c* at which the certificate
loses its sign (bisection on a monotone function; exact to 1e-10). c* = inf for
kernel-immune rows (certificate survives every c up to the search bound, or the
relevant flow is zero).

Machinery: T_bounds_from_count_intervals IMPORTED from b03_certify.py (no
copy-modify). Inputs: scout_bounds.csv (frozen deltas), scout_final_matched.csv
(certified flags). Zero heavy compute; smoke mode included via B03_P4_SMOKE.
"""

import json
import math
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, "B03_project/src")
from b03_certify import T_bounds_from_count_intervals

RES = os.environ.get("B03_P4_RES", "B03_project/results")
N_ROWS = int(os.environ.get("B03_P4_ROWS", "75"))
N_CERT = int(os.environ.get("B03_P4_CERT", "46"))
SMOKE = os.environ.get("B03_P4_SMOKE") == "1"
OUT = f"{RES}/_phase4_dryrun" if SMOKE else RES
C_MAX = 100.0
TOL = 1e-13  # bisection bracket width; residual |f(c*)| <= slope*TOL stays << 1e-10

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:.1f}s] {msg}", flush=True)


def T_at(c, row, direction):
    """Scout machinery evaluated under kernel multiplier c. direction: 'lo'/'hi'."""
    dL = (c * row.dL_min, c * row.dL_max)
    dR = (c * row.dR_min, c * row.dR_max)
    T_lo, T_hi, _ = T_bounds_from_count_intervals(
        row.nL, row.nR, dL, dR, row.nl_cells, row.nr_cells)
    return T_lo if direction == "lo" else T_hi


def bisect_root(f, lo, hi, tol):
    """Root of monotone-decreasing f in [lo, hi] with f(lo) > 0 >= f(hi)."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def main():
    sb = pd.read_csv(f"{RES}/scout_bounds.csv")
    assert len(sb) == N_ROWS
    fin = pd.read_csv(f"{RES}/scout_final_matched.csv")
    m = sb.merge(fin[["pair", "region", "certified_pos_matched",
                      "certified_neg_matched"]], on=["pair", "region"],
                 validate="1:1")
    assert len(m) == N_ROWS
    assert {c for c in ("certified_pos", "certified_neg")} <= set(m.columns)

    # ---------- Assertion A + B: self-validation at c = 1 ----------
    T_lo_1 = np.array([T_at(1.0, r, "lo") for r in m.itertuples()])
    T_hi_1 = np.array([T_at(1.0, r, "hi") for r in m.itertuples()])
    da = max(np.abs(T_lo_1 - m.T_lo.to_numpy()).max(),
             np.abs(T_hi_1 - m.T_hi.to_numpy()).max())
    assert da <= 1e-9, f"Assertion A FAILED: max |dT| at c=1 is {da}"
    # B1: interval signs match the scout's INTERVAL flags (certified_pos/neg in
    # scout_bounds.csv = T_lo>0 / T_hi<0).
    sign_ok = (((T_lo_1 > 0) == m.certified_pos.to_numpy()).all() and
               ((T_hi_1 < 0) == m.certified_neg.to_numpy()).all())
    assert sign_ok, "Assertion B1 FAILED: interval-sign mismatch at c=1"
    # B2: matched-null certifications are a SUBSET of interval certifications
    # (amendment record history: 37/9 matched vs 39/12 interval; the 5 dropped rows are
    # matched-null-insignificant, not sign-flipped).
    subset_ok = ((~m.certified_pos_matched | m.certified_pos).all() and
                 (~m.certified_neg_matched | m.certified_neg).all())
    assert subset_ok, "Assertion B2 FAILED: matched flags not a subset of interval flags"
    log(f"Assertions A+B PASSED: c=1 reproduces scout intervals (max |dT| = {da:.2e}); "
        f"interval signs match; matched set is a subset")

    # ---------- certified rows: critical multiplier c* ----------
    rows = []
    for r in m.itertuples():
        if r.certified_pos_matched:
            direction = "lo"
        elif r.certified_neg_matched:
            direction = "hi"
        else:
            continue
        f = lambda c: T_at(c, r, direction)  # monotone: lo decreasing, hi increasing
        # verify monotonicity numerically on the bracket (Assertion C)
        grid = np.linspace(1.0, C_MAX, 11)
        vals = np.array([f(c) for c in grid])
        if direction == "lo":
            assert (np.diff(vals) <= 1e-12).all(), f"Assertion C FAILED (lo) {r.pair} r{r.region}"
        else:
            assert (np.diff(vals) >= -1e-12).all(), f"Assertion C FAILED (hi) {r.pair} r{r.region}"
        v_hi = f(C_MAX)
        if direction == "lo" and v_hi > 0:
            c_star = math.inf
        elif direction == "hi" and v_hi < 0:
            c_star = math.inf
        else:
            # bisect on a DECREASING function g: for 'lo', g = f; for 'hi', g = -f
            # (T_hi is increasing; -T_hi is decreasing with -f(1) > 0 >= -f(C_MAX))
            if direction == "lo":
                g = f
            else:
                g = lambda c, _f=f: -_f(c)
            c_star = bisect_root(g, 1.0, C_MAX, TOL)
            v_root = f(c_star)
            assert abs(v_root) <= 1e-10, \
                f"Assertion D FAILED: |T(c*)| = {abs(v_root)} at {r.pair} r{r.region}"
        rows.append(dict(pair=r.pair, region=int(r.region), direction=direction,
                         T_lo_scout=r.T_lo, T_hi_scout=r.T_hi,
                         c_star=c_star,
                         control=bool(r.control)))
    cert = pd.DataFrame(rows)
    assert len(cert) == N_CERT, f"expected {N_CERT} certified rows, got {len(cert)}"
    log(f"c* computed for {len(cert)} certified rows "
        f"({int(np.isinf(cert.c_star).sum())} kernel-immune)")

    # ---------- secondary: calibration upside at c = 0.5 ----------
    gained = []
    for r in m.itertuples():
        if r.certified_pos_matched or r.certified_neg_matched:
            continue
        t_lo_half = T_at(0.5, r, "lo")
        t_hi_half = T_at(0.5, r, "hi")
        if t_lo_half > 0 or t_hi_half < 0:
            gained.append(r.pair)
    log(f"secondary: {len(gained)}/29 non-certified rows would certify at c=0.5")

    # ---------- summary ----------
    cs = cert.c_star.to_numpy()
    fin_vals = cs[np.isfinite(cs)]
    summary = dict(
        prereg="B03_PHASE4_PREREG.md v1.0 (frozen before compute)",
        smoke=SMOKE,
        n_certified=46,
        kernel_immune=int(np.isinf(cs).sum()),
        c_star_median=float(np.median(fin_vals)),
        c_star_min=float(fin_vals.min()),
        frac_ge=dict((str(t), float((cs >= t).mean()))
                     for t in (1.25, 1.5, 2.0, 3.0)),
        control_median_c_star=float(cert[cert.control].c_star.replace(
            [np.inf], np.nan).dropna().median()),
        disputed_median_c_star=float(cert[~cert.control].c_star.replace(
            [np.inf], np.nan).dropna().median()),
        secondary_certify_at_half=len(gained),
        secondary_gained_pairs=sorted(set(gained)),
        assertions=dict(A_selfvalidation=f"PASSED (max |dT| = {da:.2e})",
                        B_signs="PASSED",
                        C_monotonicity="PASSED (all rows)",
                        D_bisection="PASSED (|T(c*)| <= 1e-10)"),
    )
    os.makedirs(OUT, exist_ok=True)
    cert.to_csv(f"{OUT}/phase4_kernel_margins.csv", index=False)
    with open(f"{OUT}/phase4_summary.json", "w") as fj:
        json.dump(summary, fj, indent=2)
    log(json.dumps({k: summary[k] for k in
                    ("kernel_immune", "c_star_median", "c_star_min", "frac_ge",
                     "control_median_c_star", "disputed_median_c_star",
                     "secondary_certify_at_half")}))
    print("PHASE4 KERNEL SENSITIVITY COMPLETE", flush=True)


if __name__ == "__main__":
    main()
