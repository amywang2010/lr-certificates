"""B03: the claim layer, from stored matched-null p-values to claim flags.

The direction of a claim is read from the interval sign, so the one-sided
p-value that decides it is a selected p-value: under the null either tail could
have produced the selection, and the union bound over the two directions makes
$2p$ the valid p-value (Pr(p <= a) <= 2a).  A tempting shortcut at this point is
to argue that the Benjamini-Yekutieli harmonic factor H_m >= 2 absorbs that
factor.  It does not.  BY's guarantee is stated for valid p-values, and its
thresholds i q / (m H_m) are the BH thresholds for H_m p; the selection factor
sits outside that structure.  Charging it explicitly gives

    p_adj = min(1, 2 p_selected),      q_adj = BY(p_adj) within region,

i.e. thresholds p_selected <= i q / (2 m H_m).  A row carries a claim when its
interval excludes zero in a direction and q_adj <= 0.10 in that same direction.
BH on the undoubled selected p-values is reported alongside for comparability
only; it carries no guarantee, because the sign filter is a function of the same
data as the p-value.

Inputs are the shipped matched-null tables; p-values are not recomputed here, so
this script cannot drift from the null layer.  Two of the shipped tables store p
to three decimals, which rounds the add-one floor 1/(B+1) to 0.000; a stored zero
is therefore replaced by the smallest positive p in the same column, which is
that leg's floor.  The claim set is unchanged by the substitution.

Output per leg:

    <leg>/claim_layer.csv   pair, region, T_lo, T_hi, p_pos, p_neg, p_sel, p_adj,
                            q_adj_pos, q_adj_neg, q_bh_unadj_pos,
                            q_bh_unadj_neg, claim_pos, claim_neg

Usage:  python src/b03_claim_layer.py
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QLVL = 0.10

LEGS = [
    ("breast", os.path.join(ROOT, "results", "scout_final_matched.csv"),
     "p_pos_raw", "p_neg_raw"),
    ("lung", os.path.join(ROOT, "lung", "scout_final_matched.csv"),
     "p_pos", "p_neg"),
    ("breast_s6", os.path.join(ROOT, "breast_s6", "scout_final_matched.csv"),
     "p_pos_raw", "p_neg_raw"),
]


def bh(p):
    """BH adjusted p-values within one family."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    return out


def by(p):
    """BY adjusted p-values within one family (BH times the harmonic factor)."""
    p = np.asarray(p, float)
    m = len(p)
    return np.clip(bh(p) * np.sum(1.0 / np.arange(1, m + 1)), 0, 1)


def repair_floor(p, replaced):
    """Replace stored zeros (three-decimal rounding of the add-one floor)."""
    p = np.asarray(p, float).copy()
    zero = p == 0.0
    pos = p[~zero]
    if zero.any() and len(pos):
        p[zero] = pos.min()
        replaced.append(int(zero.sum()))
    return p


def claim_layer(df, ppos, pneg):
    out = df[["pair", "region", "T_lo", "T_hi"]].copy()
    out["sign_pos"] = out.T_lo > 0
    out["sign_neg"] = out.T_hi < 0
    zeros = []
    p_pos = repair_floor(df[ppos], zeros)
    p_neg = repair_floor(df[pneg], zeros)
    out["p_pos"] = p_pos
    out["p_neg"] = p_neg
    out["p_sel"] = np.where(out.sign_pos.to_numpy(), p_pos, p_neg)
    out["p_adj"] = np.minimum(1.0, 2.0 * out.p_sel)
    for col in ("q_adj_pos", "q_adj_neg", "q_bh_unadj_pos", "q_bh_unadj_neg"):
        out[col] = np.nan
    for _, g in out.groupby("region"):
        # one family per (region, direction): each direction gets its own
        # harmonic factor, and the doubling is charged once per row
        out.loc[g.index, "q_adj_pos"] = by(np.minimum(1.0, 2.0 * g.p_pos.to_numpy()))
        out.loc[g.index, "q_adj_neg"] = by(np.minimum(1.0, 2.0 * g.p_neg.to_numpy()))
        out.loc[g.index, "q_bh_unadj_pos"] = bh(g.p_pos.to_numpy())
        out.loc[g.index, "q_bh_unadj_neg"] = bh(g.p_neg.to_numpy())
    out["claim_pos"] = out.sign_pos & (out.q_adj_pos <= QLVL)
    out["claim_neg"] = out.sign_neg & (out.q_adj_neg <= QLVL)
    return out, zeros


def attach(leg_dir, table, cl):
    """Write the claim columns into the leg's shipped matched-null table."""
    add = cl[["pair", "region", "p_adj", "q_adj_pos", "q_adj_neg",
              "claim_pos", "claim_neg"]].copy()
    keep = [c for c in table.columns
            if c not in ("p_adj", "q_adj_pos", "q_adj_neg", "claim_pos", "claim_neg")]
    out = table[keep].merge(add, on=["pair", "region"], how="left", validate="1:1")
    assert len(out) == len(table)
    out.to_csv(os.path.join(leg_dir, "scout_final_matched.csv"), index=False)
    return out


def downstream(summary):
    """Every headline number that depends on which rows carry a claim."""
    rep = {}
    b = pd.read_csv(os.path.join(ROOT, "results", "scout_final_matched.csv"))
    cb = b.claim_pos.astype(bool) | b.claim_neg.astype(bool)
    rep["breast_rows_without_claim"] = int((~cb).sum())
    rep["breast_no_claim_straddling"] = int((~cb & ~(b.T_lo > 0) & ~(b.T_hi < 0)).sum())
    rep["breast_no_claim_with_sign"] = int(
        (~cb & ((b.T_lo > 0) | (b.T_hi < 0))).sum())
    ctrl = b[b.control.astype(bool)]
    rep["breast_control_claim_pos"] = int(ctrl.claim_pos.sum())
    rep["breast_control_claim_neg"] = int(ctrl.claim_neg.sum())
    rep["breast_control_no_claim"] = int((~ctrl.claim_pos & ~ctrl.claim_neg).sum())
    rep["breast_claim_rows_same_class"] = int(
        cb[b.pair.isin(["PECAM1->KDR", "ESR1->PGR", "PGR->ESR1", "ERBB2->EGFR"])].sum())

    conv = pd.read_csv(os.path.join(ROOT, "results", "baseline_comparison",
                                    "conventional_null_breast.csv"))
    conv["call"] = np.where(conv.q_pos <= 0.10, "positive",
                            np.where(conv.q_neg <= 0.10, "negative", "none"))
    m = b.merge(conv[["pair", "region", "call"]], on=["pair", "region"])
    m["claim"] = np.where(m.claim_pos, "positive",
                          np.where(m.claim_neg, "negative", "none"))
    rep["cross_tab"] = {
        f"{a}->{c}": int(n) for (a, c), n in
        m.groupby(["call", "claim"]).size().items()}
    rep["breast_enriched_no_claim"] = int(
        ((m.call == "positive") & (m.claim == "none")).sum())
    rep["breast_enriched_claim_negative"] = int(
        ((m.call == "positive") & (m.claim == "negative")).sum())

    cov = pd.read_csv(os.path.join(ROOT, "results", "phase2_proseg_coverage.csv"))
    cov = cov.merge(b[["pair", "region", "claim_pos", "claim_neg"]],
                    on=["pair", "region"], how="left")
    cbc = cov.claim_pos.astype(bool) | cov.claim_neg.astype(bool)
    rep["breast_coverage_claim_rows"] = f"{int(cov.inside[cbc].sum())}/{int(cbc.sum())}"
    rep["breast_coverage_all_rows"] = f"{int(cov.inside.sum())}/{len(cov)}"
    for _c in ("T", "T_lo", "T_hi", "inside"):
        cov[_c] = pd.to_numeric(cov[_c], errors="coerce")
    cov["inside"] = cov["inside"].astype(bool)
    sgn = np.sign(cov["T"].to_numpy(float)) * np.sign(
        cov.T_lo.to_numpy(float) + cov.T_hi.to_numpy(float))
    rep["breast_sign_agreement_claim_rows"] = f"{int((sgn[cbc] > 0).sum())}/{int(cbc.sum())}"

    p3 = pd.read_csv(os.path.join(ROOT, "results", "phase3_proseg_null.csv"))
    p3["q_pos"] = pd.to_numeric(p3["q_pos"], errors="coerce")
    p3["q_neg"] = pd.to_numeric(p3["q_neg"], errors="coerce")
    p3 = p3.merge(b[["pair", "region", "claim_pos", "claim_neg"]],
                  on=["pair", "region"], how="left")
    p3c = p3.claim_pos.astype(bool) | p3.claim_neg.astype(bool)
    ok = np.where(p3.claim_pos, p3.q_pos <= 0.10, p3.q_neg <= 0.10)
    rep["breast_proseg_null_pass_claim_rows"] = f"{int(ok[p3c].sum())}/{int(p3c.sum())}"

    kb = pd.read_csv(os.path.join(ROOT, "results", "phase4_kernel_margins.csv"))
    kb = kb.merge(b[["pair", "region", "claim_pos", "claim_neg"]],
                  on=["pair", "region"], how="left")
    kbc = kb[(kb.claim_pos.astype(bool)) | (kb.claim_neg.astype(bool))]
    rep["breast_margin_rows"] = int(len(kbc))
    rep["breast_margin_median"] = round(float(kbc.c_star.median()), 3)
    rep["breast_margin_min"] = round(float(kbc.c_star.min()), 3)
    rep["breast_margin_ge125"] = round(float((kbc.c_star >= 1.25).mean()), 4)
    rep["breast_margin_ge2"] = round(float((kbc.c_star >= 2.0).mean()), 4)

    nc = pd.read_csv(os.path.join(ROOT, "results", "coverage_check.csv"))
    nc = nc[nc.nucleus_outside_U.astype(bool)].merge(
        b[["pair", "region", "claim_pos", "claim_neg"]], on=["pair", "region"])
    rep["nucleus_only_rows"] = int(len(nc))
    rep["nucleus_only_claim_rows"] = int((nc.claim_pos | nc.claim_neg).sum())

    lu = pd.read_csv(os.path.join(ROOT, "lung", "scout_final_matched.csv"))
    lc = lu.claim_pos.astype(bool) | lu.claim_neg.astype(bool)
    rep["lung_claim_rows_same_class"] = int(
        lc[lu.pair.isin(["ESR1->PGR", "PGR->ESR1", "CDH1->EGFR"])].sum())
    ep = lu[lu.pair.isin(["ESR1->PGR", "PGR->ESR1"])]
    rep["lung_ESR1_PGR_max_q_adj"] = round(float(ep.q_adj_neg.max()), 4)
    lcov = pd.read_csv(os.path.join(ROOT, "lung", "phase5_e2_coverage.csv"))
    lcov = lcov.merge(lu[["pair", "region", "claim_pos", "claim_neg"]],
                      on=["pair", "region"], how="left")
    lcbc = lcov.claim_pos.astype(bool) | lcov.claim_neg.astype(bool)
    rep["lung_coverage_claim_rows"] = f"{int(lcov.inside[lcbc].sum())}/{int(lcbc.sum())}"
    rep["lung_coverage_all_rows"] = f"{int(lcov.inside.sum())}/{len(lcov)}"
    kl = pd.read_csv(os.path.join(ROOT, "lung", "phase4_kernel_margins.csv"))
    kl = kl.merge(lu[["pair", "region", "claim_pos", "claim_neg"]],
                  on=["pair", "region"], how="left")
    klc = kl[(kl.claim_pos.astype(bool)) | (kl.claim_neg.astype(bool))]
    rep["lung_margin_rows"] = int(len(klc))
    rep["lung_margin_median"] = round(float(klc.c_star.median()), 3)
    rep["lung_margin_min"] = round(float(klc.c_star.min()), 3)

    s = pd.read_csv(os.path.join(ROOT, "breast_s6", "scout_final_matched.csv"))
    rep["s6_claim_rows"] = int((s.claim_pos.astype(bool) | s.claim_neg.astype(bool)).sum())
    return rep


def main():
    summary = {}
    for name, path, ppos, pneg in LEGS:
        df = pd.read_csv(path)
        cl, zeros = claim_layer(df, ppos, pneg)
        leg_dir = os.path.dirname(path)
        cl.to_csv(os.path.join(leg_dir, "claim_layer.csv"), index=False)
        attach(leg_dir, df, cl)
        summary[name] = {
            "rows": int(len(cl)),
            "regions": int(cl.region.nunique()),
            "rows_per_family": int(len(cl) / cl.region.nunique()),
            "interval_pos": int(cl.sign_pos.sum()),
            "interval_neg": int(cl.sign_neg.sum()),
            "claim_pos": int(cl.claim_pos.sum()),
            "claim_neg": int(cl.claim_neg.sum()),
            "claim_pos_bh_unadj": int((cl.sign_pos & (cl.q_bh_unadj_pos <= QLVL)).sum()),
            "claim_neg_bh_unadj": int((cl.sign_neg & (cl.q_bh_unadj_neg <= QLVL)).sum()),
            "max_q_adj_claimed": float(
                np.max(np.r_[cl.loc[cl.claim_pos, "q_adj_pos"],
                             cl.loc[cl.claim_neg, "q_adj_neg"]], initial=0.0)),
            "stored_zeros_repaired": int(sum(zeros)),
        }
        print(f"{name}: {summary[name]}")
    summary["downstream"] = downstream(summary)
    with open(os.path.join(ROOT, "results", "claim_layer_summary.json"), "w",
              newline="\n") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    print(json.dumps(summary["downstream"], indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
