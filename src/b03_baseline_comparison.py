"""B03 baseline comparison: descriptive cross-tabs from sealed artifacts only.

Implements B03_BASELINE_COMPARISON.md (frozen 2026-09-13) exactly:
  Comparison 1 (breast): conventional point-estimate permutation calls
    (quarantined null_perm_INVALID_layer_mismatch_DEV003.csv, B=1,000) vs the
    two-layer certificates at the same screen resolution (scout_final_matched.csv).
  Comparison 2 (all sections): certified vs not x Proseg-T inside vs outside the
    certified interval, per section and pooled; Fisher exact p, odds ratio with
    95% profile likelihood CI (closed-form Woolf logit CI on the OR, reported as
    such; zero cells handled by the Haldane-Anscombe 0.5 correction, recorded).

Descriptive only: no verdict changes, no certified-set modification. Every output
row names its source artifact. Writes:
  results/baseline_comparison/comparison1_breast_conventional_vs_certified.csv
  results/baseline_comparison/comparison2_instability_by_section.csv
  results/baseline_comparison/summary.json
"""
from __future__ import annotations

import json
import math
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESEARCH = os.path.dirname(ROOT)
OUT = os.path.join(RESEARCH, "B03_project", "results", "baseline_comparison")

SOURCES = {
    "breast_final": "B03_project/results/scout_final_matched.csv",
    "breast_conventional": "B03_project/results/baseline_comparison/conventional_null_breast.csv",
    "breast_e2": "B03_project/results/phase2_proseg_coverage.csv",
    "lung_final": "xenium_lung/crop/data/scout_final_lung.csv",
    "lung_e2": "xenium_lung/crop/data/phase5_e2_coverage.csv",
    "s6_final": "xenium_breast_s6/crop/data/scout_final_s6.csv",
    "s6_e2": "xenium_breast_s6/crop/data/phase7_e2_coverage.csv",
}


def p(name: str) -> str:
    return os.path.join(RESEARCH, name)


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Exact two-sided Fisher p via the probability-of-table-in-tail sum."""
    def logchoose(n: int, k: int) -> float:
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def loghyper(x: int) -> float:
        return (logchoose(a + b, x) + logchoose(c + d, a + c - x)
                - logchoose(a + b + c + d, a + c))

    obs = loghyper(a)
    lo = max(0, a + c - (c + d))
    hi = min(a + b, a + c)
    weights = {x: math.exp(loghyper(x) - obs) for x in range(lo, hi + 1)}
    total = sum(weights.values())
    tail = sum(w for w in weights.values() if w <= 1.0 + 1e-9)
    return min(1.0, tail / total)


def woolf_or_ci(a: int, b: int, c: int, d: int):
    """Odds ratio with Woolf logit 95% CI; Haldane-Anscombe 0.5 correction if any zero cell.
    Returns (None, None, None) when a 2x2 margin is zero (OR undefined)."""
    if min(a + b, c + d, a + c, b + d) == 0:
        return None, None, None
    corrected = min(a, b, c, d) == 0
    a2, b2, c2, d2 = (a + 0.5, b + 0.5, c + 0.5, d + 0.5) if corrected else (a, b, c, d)
    orr = (a2 * d2) / (b2 * c2)
    se = math.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    lo, hi = orr * math.exp(-1.959963984540054 * se), orr * math.exp(1.959963984540054 * se)
    return orr, lo, hi


def main() -> None:
    os.makedirs(OUT, exist_ok=True)

    # ---------- Comparison 1 (breast) ----------
    conv = pd.read_csv(p(SOURCES["breast_conventional"]))[["pair", "region", "q_pos", "q_neg", "T0"]]
    conv.columns = ["pair", "region", "qpos_conv", "qneg_conv", "T0_conv"]
    fin = pd.read_csv(p(SOURCES["breast_final"]))
    assert len(conv) == 75 and len(fin) == 75, (len(conv), len(fin))
    m = fin.merge(conv, on=["pair", "region"])
    assert len(m) == 75, "merge lost rows: key mismatch between final and conventional CSVs"
    assert (m["T0"] - m["T0_conv"]).abs().max() < 1e-9, "T0 mismatch vs conventional run"

    def conv_call(r) -> str:
        if r["qpos_conv"] <= 0.10:
            return "positive"
        if r["qneg_conv"] <= 0.10:
            return "negative"
        return "none"

    m["conv_call"] = m.apply(conv_call, axis=1)
    # two-layer certificates from the final CSV's verdict flags
    m["cert_call"] = "none"
    m.loc[m["certified_pos_matched"].astype(bool), "cert_call"] = "positive"
    m.loc[m["certified_neg_matched"].astype(bool), "cert_call"] = "negative"
    assert not (m["certified_pos_matched"].astype(bool)
                & m["certified_neg_matched"].astype(bool)).any()
    ct1 = pd.crosstab(m["conv_call"], m["cert_call"])
    conv_pos_any = int((m["conv_call"] != "none").sum())
    conv_pos_not_cert = int(((m["conv_call"] != "none") & (m["cert_call"] == "none")).sum())
    conv_pos_cert_dir = int(
        (
            ((m["conv_call"] == "positive") & (m["cert_call"] == "positive"))
            | ((m["conv_call"] == "negative") & (m["cert_call"] == "negative"))
        ).sum()
    )
    conv_pos_flip = int(
        (
            ((m["conv_call"] == "positive") & (m["cert_call"] == "negative"))
            | ((m["conv_call"] == "negative") & (m["cert_call"] == "positive"))
        ).sum()
    )
    flip_rows = m[
        ((m["conv_call"] == "positive") & (m["cert_call"] == "negative"))
        | ((m["conv_call"] == "negative") & (m["cert_call"] == "positive"))
    ][["pair", "region", "T0"]]
    c1 = ct1.reset_index().melt(id_vars="conv_call", var_name="cert_call", value_name="n")
    c1["source"] = SOURCES["breast_conventional"] + " + " + SOURCES["breast_final"]
    c1.to_csv(os.path.join(OUT, "comparison1_breast_conventional_vs_certified.csv"), index=False)
    flip_rows.to_csv(os.path.join(OUT, "comparison1_direction_flips.csv"), index=False)

    # ---------- Comparison 2 (per section + pooled) ----------
    rows = []
    pooled = {"cert": 0, "inst": 0, "cert_inst": 0, "n": 0}
    for section, sf, ef in [
        ("breast", SOURCES["breast_final"], SOURCES["breast_e2"]),
        ("lung", SOURCES["lung_final"], SOURCES["lung_e2"]),
        ("breast_s6", SOURCES["s6_final"], SOURCES["s6_e2"]),
    ]:
        e = pd.read_csv(p(ef))
        inst_col = "inside"
        e["instable"] = ~e[inst_col].astype(bool)
        key = ["pair", "region"]
        if section == "breast_s6":
            # s6's coverage file has no verdict columns; merge from the final CSV and
            # assert the certified set is empty (A1-g, frozen).
            f = pd.read_csv(p(sf))
            merged = e.merge(f[key + ["certified_pos_matched", "certified_neg_matched"]],
                             on=key, how="left")
            merged["certified"] = (merged["certified_pos_matched"].astype(bool)
                                   | merged["certified_neg_matched"].astype(bool))
            assert merged["certified"].sum() == 0, "A1-g violation: s6 certified rows exist"
        elif section == "breast":
            merged = e.copy()
            merged["certified"] = (merged["certified_pos_matched"].astype(bool)
                                   | merged["certified_neg_matched"].astype(bool))
        else:
            merged = e.copy()
            merged["certified"] = (merged["cert_pos_final"].astype(bool)
                                   | merged["cert_neg_final"].astype(bool))
        a = int((merged["certified"] & merged["instable"]).sum())
        b = int((merged["certified"] & ~merged["instable"]).sum())
        c = int((~merged["certified"] & merged["instable"]).sum())
        d = int((~merged["certified"] & ~merged["instable"]).sum())
        n = a + b + c + d
        assert n == len(merged) == (75 if section == "breast" else 54 if section == "lung" else 30)
        fpr = fisher_two_sided(a, b, c, d)
        orr, orlo, orhi = woolf_or_ci(a, b, c, d)
        r4 = (lambda v: None if v is None else round(v, 4))
        rows.append({
            "section": section, "n": n, "certified": a + b, "not_certified": c + d,
            "certified_unstable": a, "certified_stable": b,
            "notcert_unstable": c, "notcert_stable": d,
            "fisher_p": round(fpr, 6), "odds_ratio": r4(orr),
            "or_lo95": r4(orlo), "or_hi95": r4(orhi),
            "zero_cell_corrected": min(a, b, c, d) == 0,
            "source": ef,
        })
        pooled["cert"] += a + b
        pooled["inst"] += c
        pooled["cert_inst"] += a
        pooled["n"] += n
    a = pooled["cert_inst"]
    b_ = pooled["cert"] - pooled["cert_inst"]
    c_ = pooled["inst"]
    d_ = pooled["n"] - pooled["cert"] - c_
    assert a + b_ + c_ + d_ == pooled["n"] and min(a, b_, c_, d_) >= 0
    rows.append({
        "section": "pooled", "n": pooled["n"], "certified": pooled["cert"],
        "not_certified": pooled["n"] - pooled["cert"],
        "certified_unstable": a, "certified_stable": b_,
        "notcert_unstable": c_, "notcert_stable": d_,
        "fisher_p": round(fisher_two_sided(a, b_, c_, d_), 6),
        "odds_ratio": r4(woolf_or_ci(a, b_, c_, d_)[0]),
        "or_lo95": r4(woolf_or_ci(a, b_, c_, d_)[1]),
        "or_hi95": r4(woolf_or_ci(a, b_, c_, d_)[2]),
        "zero_cell_corrected": min(a, b_, c_, d_) == 0,
        "source": "phase2_proseg_coverage.csv + phase5_e2_coverage.csv + phase7_e2_coverage.csv",
    })
    c2 = pd.DataFrame(rows)
    c2.to_csv(os.path.join(OUT, "comparison2_instability_by_section.csv"), index=False)

    summary = {
        "protocol": "B03_BASELINE_COMPARISON.md v1.0 (frozen 2026-09-13)",
        "comparison1": {
            "n": 75,
            "conventional_calls_any": conv_pos_any,
            "conventional_pos_not_certified": conv_pos_not_cert,
            "conventional_pos_same_direction_certified": conv_pos_cert_dir,
            "conventional_pos_direction_flip": conv_pos_flip,
            "flip_rows": flip_rows.to_dict("records"),
            "crosstab": ct1.to_dict(),
        },
        "comparison2": rows,
    }
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
