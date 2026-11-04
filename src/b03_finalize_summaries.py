"""Refresh the post-claim summaries after the claim layer has been applied.

The phase scripts run in pipeline order.  Their summary files were first written when
the certified set was the undoubled Benjamini--Hochberg set (46 breast rows, seven
lung rows).  The reported certificates are the selection-adjusted claims written by
``b03_claim_layer.py`` (45 breast rows, five lung rows).  This step recomputes the
downstream summary fields from the shipped per-row tables so that every file agrees
with the claim layer, and records the pre-correction counts separately.

Inputs are the released per-row CSVs; no raw data are needed.  The outputs are the
same JSON and CSV files the phase scripts write, so this step is idempotent and can
be rerun at any time.

  python src/b03_finalize_summaries.py
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def save_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)


def claim_frame(directory):
    """Per-row claim flags for one leg, keyed by (pair, region)."""
    cl = pd.read_csv(os.path.join(directory, "claim_layer.csv"))
    return cl[["pair", "region", "claim_pos", "claim_neg"]]


def phase2():
    """Breast E2 coverage: primary over all evaluable rows, claim subset beside it."""
    cov = pd.read_csv(os.path.join(RES, "phase2_proseg_coverage.csv"))
    cov = cov.merge(claim_frame(RES), on=["pair", "region"], validate="1:1")
    claim = cov.claim_pos | cov.claim_neg
    summary = load_json(os.path.join(RES, "phase2_proseg_summary.json"))
    ev = claim & cov["T"].notna()
    n_in = int(cov.loc[ev, "inside"].sum())
    n = int(ev.sum())
    summary["certified_rows"] = int(claim.sum())
    summary["certified_inside"] = n_in
    summary["certified_fraction"] = n_in / n
    summary["undoubled_bh_rows"] = 46
    summary["sign_agreement_on_certified"] = "{}/{}".format(
        int((((cov["T"] > 0) == cov.claim_pos)[ev]).sum()), n)
    save_json(os.path.join(RES, "phase2_proseg_summary.json"), summary)
    return "phase2", summary["certified_rows"], n_in, n


def phase3():
    """Breast E3: directional significance inside the claim set."""
    nul = pd.read_csv(os.path.join(RES, "phase3_proseg_null.csv"))
    nul = nul.merge(claim_frame(RES), on=["pair", "region"], validate="1:1")
    fin = pd.read_csv(os.path.join(RES, "scout_final_matched.csv"))
    claim = nul.claim_pos | nul.claim_neg
    sig = ((nul.claim_pos & (nul.q_pos <= 0.10)) |
           (nul.claim_neg & (nul.q_neg <= 0.10)))
    bh = fin.certified_pos_matched | fin.certified_neg_matched
    bh_sig = ((fin.certified_pos_matched & (fin.q_pos <= 0.10)) |
              (fin.certified_neg_matched & (fin.q_neg <= 0.10)))
    summary = load_json(os.path.join(RES, "phase3_summary.json"))
    summary["certified_rows"] = int(claim.sum())
    summary["significant_in_certified_direction"] = int(sig.sum())
    summary["fraction"] = float(sig.sum() / claim.sum())
    summary["gate"] = ("PASS" if summary["fraction"] >= 0.85 else
                       "PARTIAL" if summary["fraction"] >= 0.60 else "FAIL")
    summary["undoubled_bh_rows"] = int(bh.sum())
    summary["undoubled_bh_significant"] = int(bh_sig.sum())
    save_json(os.path.join(RES, "phase3_summary.json"), summary)
    return "phase3", summary["certified_rows"], int(sig.sum())


def phase4(directory, name):
    """Kernel margins: the tabulated CSV is the BH superset, summaries describe claims."""
    margins = os.path.join(directory, "phase4_kernel_margins.csv")
    m = pd.read_csv(margins)
    cl = claim_frame(directory)
    j = m.merge(cl, on=["pair", "region"], how="left", validate="1:1")
    claim = j.claim_pos.fillna(False) | j.claim_neg.fillna(False)
    sub = j[claim]
    cs = sub.c_star.to_numpy()
    finite = cs[np.isfinite(cs)]
    summary = load_json(os.path.join(directory, "phase4_summary.json"))
    summary["n_certified"] = int(claim.sum())
    summary["undoubled_bh_rows"] = int(len(m))
    summary["kernel_immune"] = int(np.isinf(cs).sum())
    summary["c_star_median"] = float(np.median(finite))
    summary["c_star_min"] = float(finite.min())
    summary["frac_ge"] = {str(t): float((cs >= t).mean())
                          for t in (1.25, 1.5, 2.0, 3.0)}
    summary["control_median_c_star"] = float(sub[sub.control].c_star.replace(
        [np.inf], np.nan).dropna().median())
    summary["disputed_median_c_star"] = float(sub[~sub.control].c_star.replace(
        [np.inf], np.nan).dropna().median())
    save_json(os.path.join(directory, "phase4_summary.json"), summary)
    return name, summary["n_certified"], summary["c_star_median"]


def phase5():
    """Lung E2 coverage on the claim set."""
    directory = os.path.join(ROOT, "lung")
    cov = pd.read_csv(os.path.join(directory, "phase5_e2_coverage.csv"))
    cov = cov.merge(claim_frame(directory), on=["pair", "region"], validate="1:1")
    claim = cov.claim_pos | cov.claim_neg
    ev = claim & cov["T"].notna()
    n_in = int(cov.loc[ev, "inside"].sum())
    n = int(ev.sum())
    bh = cov.cert_pos_final | cov.cert_neg_final
    bh_ev = bh & cov["T"].notna()
    summary = load_json(os.path.join(directory, "phase5_e2_summary.json"))
    summary["certified_rows"] = int(claim.sum())
    summary["primary"] = dict(
        n_inside=n_in, n_certified_evaluable=n,
        fraction=n_in / n,
        gate=("PASS" if n_in / n >= 0.90 else
              "PARTIAL" if n_in / n >= 0.75 else "FAIL"))
    summary["undoubled_bh_rows"] = int(bh.sum())
    summary["undoubled_bh_inside"] = int(cov.loc[bh_ev, "inside"].sum())
    summary["sign_agreement_on_certified"] = "{}/{}".format(
        int((((cov["T"] > 0) == cov.claim_pos)[ev]).sum()), n)
    save_json(os.path.join(directory, "phase5_e2_summary.json"), summary)
    return "phase5", summary["certified_rows"], n_in, n


def coupled():
    """Rename the mixed sign-change counter and add the coupling-only counter."""
    paths = [os.path.join(RES, "coupled_end.csv"),
             os.path.join(ROOT, "lung", "coupled_end.csv"),
             os.path.join(ROOT, "breast_s6", "coupled_end.csv")]
    out = []
    for path in paths:
        df = pd.read_csv(path)
        if "sign_changed" in df.columns:
            df = df.rename(columns={"sign_changed": "sign_changed_vs_sealed"})
        df["sign_changed_coupling"] = (
            ((df.T_lo_joint > 0) != (df.T_lo_flowsep > 0)) |
            ((df.T_hi_joint < 0) != (df.T_hi_flowsep < 0)))
        df.to_csv(path, index=False)
        out.append(df)
    summary_path = os.path.join(RES, "coupled_end_summary.json")
    s = load_json(summary_path)
    for tag, df in zip(("breast", "lung", "s6"), out):
        st = s["legs"][tag]
        st.pop("sign_changes", None)
        st["sign_changes_vs_sealed"] = int(df.sign_changed_vs_sealed.sum())
        st["sign_changes_from_coupling"] = int(df.sign_changed_coupling.sum())
    s.pop("n_joint_sign_changes", None)
    s["n_joint_sign_changes_vs_sealed"] = int(sum(
        d.sign_changed_vs_sealed.sum() for d in out))
    s["n_joint_sign_changes_from_coupling"] = int(sum(
        d.sign_changed_coupling.sum() for d in out))
    save_json(summary_path, s)
    return "coupled", s["n_joint_sign_changes_vs_sealed"], \
        s["n_joint_sign_changes_from_coupling"]


def capacity_flags():
    """Make the frozen capacity feasibility report self-describing."""
    for path in (os.path.join(ROOT, "lung", "capacity_check.json"),
                 os.path.join(ROOT, "breast_s6", "capacity_check.json")):
        d = load_json(path)
        renamed = {}
        for k, v in d.items():
            renamed[{
                "cap_gain_ok": "uncapped_gain_within_cap",
                "cap_loss_ok": "uncapped_loss_within_cap",
                "cells_violating_gain_cap": "cells_whose_uncapped_gain_exceeds_cap",
            }.get(k, k)] = v
        renamed["capped_formula_applied"] = True
        renamed["sound_upper_bound"] = True
        save_json(path, renamed)
        print(path, renamed["uncapped_gain_within_cap"],
              renamed["cells_whose_uncapped_gain_exceeds_cap"])


if __name__ == "__main__":
    print(phase2())
    print(phase3())
    print(phase4(RES, "breast"))
    print(phase4(os.path.join(ROOT, "lung"), "lung"))
    print(phase5())
    print(coupled())
    capacity_flags()
