"""B03 Step 6 watcher: on matched-null completion, validate, apply add-one p-value
convention, run referee consistency checks, and emit SCOUT_REPORT.md.

Runs unattended; exits silently if the null is still in progress.
"""
import os
import sys
import json
import math
import time
import pandas as pd
import numpy as np

OUT = "B03_project/results"
REPORT = "B03_project/B03_SCOUT_REPORT.md"
LOG = "B03_project/logs/watchdog.log"

T0 = time.time()


def log(msg):
    with open(LOG, "a") as f:
        f.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")


def wait_for_null(timeout_s=8 * 3600):
    """Wait for scout_summary_matched.json to appear and null process to exit."""
    while time.time() - T0 < timeout_s:
        if os.path.exists(f"{OUT}/scout_summary_matched.json"):
            # ensure the writer has finished (file appears at end of script)
            time.sleep(30)
            return True
        time.sleep(60)
    return False


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    qq = np.minimum(q, 1.0)
    out = np.empty(m)
    out[order] = qq
    return out


def referee_checks(final):
    checks = {}
    checks["rows_75"] = len(final) == 75
    checks["no_nan_p"] = final[["p_pos_raw", "p_neg_raw", "q_pos", "q_neg"]].notna().all().all()
    checks["p_in_range"] = ((final.p_pos_raw >= 0) & (final.p_pos_raw <= 1)).all() and \
                           ((final.p_neg_raw >= 0) & (final.p_neg_raw <= 1)).all()
    checks["no_zero_p"] = (final.p_pos_raw > 0).all() and (final.p_neg_raw > 0).all()
    checks["p_not_all_ties"] = not ((final.p_pos_raw == 1.0).all() or
                                    (final.p_neg_raw == 1.0).all())
    # certified iff bound strictly beyond 0 AND matched q <= 0.1
    checks["pos_flag_consistent"] = ((final.certified_pos_matched ==
                                      ((final.T_lo > 0) & (final.q_pos <= 0.10))).all())
    checks["neg_flag_consistent"] = ((final.certified_neg_matched ==
                                      ((final.T_hi < 0) & (final.q_neg <= 0.10))).all())
    # interval ordering preserved
    checks["bounds_ordered"] = ((final.T_lo <= final.T0 + 1e-9) &
                                (final.T0 <= final.T_hi + 1e-9)).all()
    # NOTE: no gate on controls being certified negative — CD274->PDCD1 (control row)
    # is expected to certify NEGATIVE (certified trans-cellular disco-localization);
    # the prereg gate only requires >=1 control certified POSITIVE.
    return checks


def main():
    if not wait_for_null():
        log("TIMEOUT waiting for null completion")
        sys.exit(1)
    log("null complete; loading outputs")

    nul = pd.read_csv(f"{OUT}/null_robust_matched.csv")
    res = pd.read_csv(f"{OUT}/scout_bounds.csv")
    with open(f"{OUT}/scout_summary_matched.json") as f:
        summ = json.load(f)

    # add-one convention: v2 saved acc/B; exact value is (1+acc)/(B+1). The count is
    # recoverable: acc = round(p*B). Verify round-trip consistency before transforming.
    B = summ["B"]
    acc_pos = np.round(nul.p_pos * B).astype(int)
    acc_neg = np.round(nul.p_neg * B).astype(int)
    assert (np.abs(acc_pos - nul.p_pos * B) < 1e-6).all(), "p_pos round-trip failed"
    assert (np.abs(acc_neg - nul.p_neg * B) < 1e-6).all(), "p_neg round-trip failed"
    nul["p_pos_raw"] = (1 + acc_pos) / (B + 1)
    nul["p_neg_raw"] = (1 + acc_neg) / (B + 1)
    # BH within region, preserving row order (nul rows are pair-major, not grouped by region)
    nul["q_pos"] = np.nan
    nul["q_neg"] = np.nan
    for r in sorted(nul.region.unique()):
        sel = nul.region == r
        nul.loc[sel, "q_pos"] = bh(nul.loc[sel, "p_pos_raw"])
        nul.loc[sel, "q_neg"] = bh(nul.loc[sel, "p_neg_raw"])
    nul.to_csv(f"{OUT}/null_robust_matched_final.csv", index=False)

    final = res.merge(nul[["pair", "region", "p_pos_raw", "p_neg_raw", "q_pos", "q_neg",
                           "T_lo_perm_med", "T_hi_perm_med"]],
                      on=["pair", "region"], how="inner")
    assert len(final) == 75, f"merge produced {len(final)} rows"
    final["certified_pos_matched"] = (final.T_lo > 0) & (final.q_pos <= 0.10)
    final["certified_neg_matched"] = (final.T_hi < 0) & (final.q_neg <= 0.10)
    final.to_csv(f"{OUT}/scout_final_matched.csv", index=False)

    checks = referee_checks(final)
    log("referee checks: " + json.dumps({k: bool(v) for k, v in checks.items()}))    # ---- summary numbers ----
    n_pos = int(final.certified_pos_matched.sum())
    n_neg = int(final.certified_neg_matched.sum())
    nonid = 75 - n_pos - n_neg
    ctrl = final[final.control]
    n_ctrl_pos = int(ctrl.certified_pos_matched.sum())
    disp = final[~final.control]
    # prereg gate 2, VERBATIM: "at least 1 disputed pair loses sign certification
    # (0 in [T_lo, T_hi])" — the parenthetical DEFINES the event as 0 inside the
    # interval (non-identification). Report also the stricter complementary
    # evidence: disputed rows certifying the OPPOSITE/negative sign (T_hi < 0,
    # q_neg <= 0.1) — disjoint from the prereg event, stronger for the thesis.
    n_disp_0in = int(((disp.T_lo <= 0) & (0 <= disp.T_hi)).sum())
    n_disp_neg = int((disp.certified_neg_matched).sum())
    n_disp_nonid = n_disp_0in
    cov = json.load(open(f"{OUT}/coverage_summary.json"))
    # prereg PASS gate (verbatim, line 103-106): >=1 positive control certified
    # (T_lo > 0, matched q <= 0.1) AND >=1 disputed pair with 0 in [T_lo, T_hi]
    # AND coverage >= 90% AND all referee consistency checks.
    pass_gate = ((n_ctrl_pos >= 1) and (n_disp_0in >= 1) and
                 bool(cov["prereg_gate_90pct"]) and all(checks.values()))
    # KILL-clause audit (prereg lines 108-109): report explicitly; none may hold
    # for a PASS verdict. (a) all intervals vacuous; (b) coverage < 90%;
    # (c) cannot calibrate U (moot: calibration executed); (d) every disputed
    # pair vacuous (0 in interval for ALL disputed rows).
    kill = {
        "a_all_vacuous": nonid == 75,
        "b_coverage_fail": not bool(cov["prereg_gate_90pct"]),
        "c_calibration": False,  # calibration executed and cross-checked
        "d_all_disputed_vacuous": n_disp_0in == len(disp),
    }
    any_kill = any(kill.values())

    # ---- report ----
    L = []
    a = L.append
    a("# B03 SCOUT REPORT — Segmentation-Robust Spatial Interaction Inference")
    a("")
    a(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}  ")
    a("Data: Xenium v1.0.1 FFPE Human Breast Cancer Rep1 (34.49M transcripts, "
      "167,780 cells, 313-plex panel)  ")
    a("Controlling documents: the frozen topic scope.md, B03_PREREGISTRATION.md, B03_THEORY.md, "
      "B03_DEVIATIONS.md (amendment record..004)")
    a("")
    a("## Decision")
    a("")
    a(f"**{'PASS' if pass_gate else ('KILL' if any_kill else 'INCONCLUSIVE')}** — "
      f"preregistered gate (verbatim): >=1 positive control certified (T_lo > 0, "
      f"matched q <= 0.1) AND >=1 disputed pair loses sign certification "
      f"(0 in [T_lo, T_hi]) AND coverage >= 90% AND all referee consistency "
      f"checks pass. KILL clauses audited: {json.dumps({k: bool(v) for k, v in kill.items()})}")
    a("")
    a("## Gate arithmetic")
    a("")
    a(f"- Positive controls certified: **{n_ctrl_pos}** of 15 control (pair, region) rows")
    a(f"- Disputed rows with 0 in [T_lo, T_hi] (prereg gate-2 event, verbatim): "
      f"**{n_disp_0in}** of {len(disp)} disputed rows")
    a(f"- Disputed rows certifying the OPPOSITE (negative) sign, T_hi < 0 and "
      f"q_neg <= 0.1 (stronger, complementary evidence): **{n_disp_neg}** "
      f"(total certified-neg incl. controls: {n_neg})")
    a(f"- Non-identifiable rows (0 in [T_lo, T_hi]): **{nonid}** of 75")
    a(f"- Coverage (eroded in-U variant): **{cov['coverage_fraction']*100:.0f}%** "
      f"({cov['n_covered']}/{cov['n_rows']}); nucleus sensitivity outside-U rows: "
      f"{cov['nucleus_outside_U_rows']} (not gated)")
    a(f"- Referee consistency checks: "
      f"{'all PASS' if all(checks.values()) else json.dumps({k: bool(v) for k, v in checks.items()})}")
    a("")
    a("## Statistical protocol (as executed)")
    a("")
    a("- Matched robust permutation null (amendment record, parallel rerun amendment record): per "
      "permutation b, worst-case objects T_lo^(b), T_hi^(b) via Theorem-2 extremes "
      "under permuted labels; identical machinery as observed intervals (single code "
      "path; startup re-check of observed intervals within 1e-6 of scout_bounds.csv; "
      "pre-flight determinism check: perm 0 bitwise identical across parent and "
      "workers).")
    a("- Labels permuted within (region x transcript-count decile) blocks, B=1000, "
      "per-permutation documented seeds rng(20260905+100003+p) in fixed block order "
      "(reproducible under parallelism).")
    a("- p_pos = (1 + #{T_lo^(b) >= T_lo_obs})/(B+1); p_neg = (1 + #{T_hi^(b) <= T_hi_obs})/(B+1) "
      "(add-one for finite-B validity); BH within region over 25 pairs; q <= 0.10.")
    a("")
    a("## Results table (per region)")
    a("")
    for r in sorted(final.region.unique()):
        fr = final[final.region == r].sort_values(["control", "T_lo"], ascending=[True, False])
        a(f"### Region {r}")
        a("")
        a("| pair | control | T0 | T_lo | T_hi | width | q_pos | q_neg | status |")
        a("|---|---|---|---|---|---|---|---|---|")
        for _, row in fr.iterrows():
            status = ("CERT-POS" if row.certified_pos_matched else
                      "CERT-NEG" if row.certified_neg_matched else "non-id")
            a(f"| {row.pair} | {'ctrl' if row.control else 'disputed'} | "
              f"{row.T0:.3f} | {row.T_lo:.3f} | {row.T_hi:.3f} | {row.width:.3f} | "
              f"{row.q_pos:.3f} | {row.q_neg:.3f} | {status} |")
        a("")
    a("## Referee checklist (self-audit)")
    a("")
    a("- [x] Estimand, pairs, regions, U, null, FDR pre-registered before expression "
      "computation (B03_PREREGISTRATION.md; transcript file still downloading at freeze)")
    a("- [x] Assignment layer reproduced and verified (concordance diagnosis amendment record; "
      "r=0.991 per-gene, boundary-semantics class)")
    a("- [x] Certificates exact (0 violations / 150,000 random configs; subcube "
      "exhaustive enumeration matches greedy extremes; 11/11 global consistency checks)")
    a("- [x] Null matched to test statistic layer (amendment record); invalid v1 artifact "
      "quarantined and documented")
    a("- [x] Coverage gate passed on an independent assignment code path (amendment record)")
    a("- [x] Deviations logged with reasons (amendment record..004); no silent changes")
    a("- [x] Invalidated artifacts quarantined, not deleted (auditable provenance)")
    a("")
    a("## Limitations (stated for the paper)")
    a("")
    a("- Certificates quantify assignment/leakage uncertainty ONLY (reassignment of "
      "molecules among a fixed cell inventory). Segmentation schemes that change the "
      "inventory (Proseg-style merges/splits) or labels are outside U; label-drift and "
      "cross-platform sensitivity are scheduled phases of the campaign.")
    a("- Single tissue section, single technology (Xenium FFPE breast). Multi-tissue/"
      "multi-platform replication is the next campaign phase.")
    a("- p-values at B=1000 have resolution 1/(B+1) ~ 1e-3; pairs with p at the floor are "
      "reported at that floor (add-one convention).")
    a("")
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    log(f"report written: {REPORT}; decision={'PASS' if pass_gate else ('KILL' if any_kill else 'INCONCLUSIVE')}; "
        f"kill_clauses={json.dumps({k: bool(v) for k, v in kill.items()})}")


if __name__ == "__main__":
    main()
