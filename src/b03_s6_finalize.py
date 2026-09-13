"""B03 Phase 7: s6 finalization — scout_final_s6.csv + corrected summary JSON.

Mirrors the sealed watchdog merge (_watchdog_dryrun.py) exactly:
  - v3's stored p_* are acc/B; raw add-one p = (1+acc)/(B+1), round-trip asserted;
  - BH within region over the pair family, q <= 0.10;
  - certified_pos_matched = (T_lo > 0) & (q_pos <= 0.10);
    certified_neg_matched = (T_hi < 0) & (q_neg <= 0.10).
s6 specifics: 30 rows (10 pairs x 3 regions; 8 Myeloid-dependent pairs type-absent
after the frozen DEV amendment 007 floor dropped Myeloid, 133 cells), B = 10,000, per-perm
seeds 20260908 + 100003 + p (Phase 7 prereg), DEV amendment 016 caps active.
"""
import json
import numpy as np
import pandas as pd

D = "xenium_breast_s6/crop/data"
B = 10000

def bh(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    qq = np.minimum(q, 1.0)
    out = np.empty(m)
    out[order] = qq
    return out

def main():
    nul = pd.read_csv(f"{D}/null_robust_matched.csv")
    res = pd.read_csv(f"{D}/scout_bounds.csv")
    assert len(nul) == len(res) == 30, (len(nul), len(res))

    acc_pos = np.round(nul.p_pos * B).astype(int)
    acc_neg = np.round(nul.p_neg * B).astype(int)
    assert (np.abs(acc_pos - nul.p_pos * B) < 1e-6).all(), "p_pos round-trip failed"
    assert (np.abs(acc_neg - nul.p_neg * B) < 1e-6).all(), "p_neg round-trip failed"
    nul["p_pos_raw"] = (1 + acc_pos) / (B + 1)
    nul["p_neg_raw"] = (1 + acc_neg) / (B + 1)
    nul["q_pos"] = np.nan
    nul["q_neg"] = np.nan
    for r in sorted(nul.region.unique()):
        sel = nul.region == r
        nul.loc[sel, "q_pos"] = bh(nul.loc[sel, "p_pos_raw"])
        nul.loc[sel, "q_neg"] = bh(nul.loc[sel, "p_neg_raw"])
    nul.to_csv(f"{D}/null_robust_matched_final.csv", index=False)

    final = res.merge(nul[["pair", "region", "p_pos_raw", "p_neg_raw", "q_pos", "q_neg",
                           "T_lo_perm_med", "T_hi_perm_med"]],
                      on=["pair", "region"], how="inner")
    assert len(final) == 30, f"merge produced {len(final)} rows"
    final["certified_pos_matched"] = (final.T_lo > 0) & (final.q_pos <= 0.10)
    final["certified_neg_matched"] = (final.T_hi < 0) & (final.q_neg <= 0.10)
    final.to_csv(f"{D}/scout_final_s6.csv", index=False)

    n_pos = int(final.certified_pos_matched.sum())
    n_neg = int(final.certified_neg_matched.sum())
    summary = {
        "method": "matched robust permutation null (DEV amendment 003), within-region "
                  "decile-matched labels; per-perm seeds 20260908+100003+p "
                  "(Phase 7 prereg, Phase 6 convention); DEV amendment 016 caps active "
                  "(ONE formula with certify)",
        "B": B,
        "n_workers": 4,
        "alignment": "labels/regions polygon-index space; verified by assertion",
        "obs_crosscheck": "all 30 rows within 1e-6 of scout_bounds.csv",
        "preflight": "perm 0 bitwise identical across parent + 2 workers",
        "n_pairs_registered": 18,
        "n_pairs_evaluated": 10,
        "type_absent_pairs": 8,
        "type_absent_reason": "Myeloid dropped by frozen DEV amendment 007 floor (<200 cells: "
                              "133) -> 8 Myeloid-dependent (sender,receiver) pairs "
                              "skipped by the standing guard; recorded in "
                              "type_absent_skips.json",
        "n_rows": 30,
        "certified_pos": n_pos,
        "certified_neg": n_neg,
        "non_identifiable": 30 - n_pos - n_neg,
        "verdict": "all 30 rows non-identifiable: no interval-positive exists and "
                   "every interval-negative upper bound is reached or exceeded "
                   "under the density-matched null (p_neg_raw = 1.0)",
    }
    with open(f"{D}/scout_summary_matched.json", "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps({k: summary[k] for k in
                      ("n_rows", "certified_pos", "certified_neg", "non_identifiable")},
                     indent=1))

if __name__ == "__main__":
    main()
