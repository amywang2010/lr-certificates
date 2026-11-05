"""Baseline comparison: descriptive cross-tabs built from the released artifacts.

The two comparisons:
  Comparison 1 (breast): conventional point-estimate permutation calls
    (a point-statistic null at B=1,000, computed by
    src/b03_baseline_conventional_null.py) vs the
    selection-charged claims at the same screen resolution (scout_final_matched.csv).
  Comparison 2 (all sections): certified vs not x Proseg-T inside vs outside the
    certified interval, per section and pooled. Fisher exact p, odds ratio with
    95% profile likelihood CI (closed-form Woolf logit CI on the OR, with zero
    cells handled by the Haldane-Anscombe 0.5 correction, recorded).

Descriptive only: no verdict changes, no certified-set modification. Every output
row names its source artifact. Writes:
  baseline_comparison/comparison1_breast_conventional_vs_certified.csv
  baseline_comparison/comparison2_instability_by_section.csv
  baseline_comparison/summary.json
"""
from __future__ import annotations

import json
import math
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "baseline_comparison")

SOURCES = {
    "breast_final": "results/breast/scout_final_matched.csv",
    "breast_conventional": "baseline_comparison/conventional_null_breast.csv",
    "breast_e2": "coverage/phase2_proseg_coverage.csv",
    "lung_final": "results/lung/scout_final_lung.csv",
    "lung_e2": "coverage/phase5_e2_coverage.csv",
    "s6_final": "results/breast_s6/scout_final_s6.csv",
    "s6_e2": "coverage/breast_s6_phase7_e2_coverage.csv",
    "breast_claims": "results/breast/claim_layer.csv",
    "lung_claims": "results/lung/claim_layer.csv",
    "s6_claims": "results/breast_s6/claim_layer.csv",
}

# The certificate set of record is the selection-charged claim layer, not the
# pre-doubling flags stored in the matched-null tables.
CLAIM_FOR = {"breast": "breast_claims", "lung": "lung_claims", "breast_s6": "s6_claims"}


def p(name: str) -> str:
    return os.path.join(ROOT, name)


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
    # claims from the selection-charged layer, which the final
    # table already carries (written by b03_claim_layer.py); the pre-doubling flags are
    # still present in that file and must not be used here.
    assert len(m) == 75 and m["claim_pos"].notna().all(), "claim layer does not cover the 75 rows"
    assert int(m["claim_pos"].astype(bool).sum()) == 36, "claim layer drifted from 36 positives"
    m["cert_call"] = "none"
    m.loc[m["claim_pos"].astype(bool), "cert_call"] = "positive"
    m.loc[m["claim_neg"].astype(bool), "cert_call"] = "negative"
    assert not (m["claim_pos"].astype(bool) & m["claim_neg"].astype(bool)).any()
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
        # the certificate set of record is the claim layer of that leg
        merged = e.merge(
            pd.read_csv(p(SOURCES[CLAIM_FOR[section]]))[key + ["claim_pos", "claim_neg"]],
            on=key, how="left", validate="1:1")
        assert merged["claim_pos"].notna().all(), f"{section}: claim layer does not cover the leg"
        merged["certified"] = (merged["claim_pos"].astype(bool)
                               | merged["claim_neg"].astype(bool))
        if section == "breast_s6":
            assert merged["certified"].sum() == 0, "frozen rule violated: s6 certified rows exist"
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
        "source": SOURCES["breast_e2"] + " + " + SOURCES["lung_e2"] + " + " + SOURCES["s6_e2"],
    })
    c2 = pd.DataFrame(rows)
    c2.to_csv(os.path.join(OUT, "comparison2_instability_by_section.csv"), index=False)

    summary = {
        "protocol": "conventional co-expression permutation baseline against the selection-charged claims, descriptive cross-tabs (paper Table 1)",
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
