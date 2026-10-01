"""B03 synthetic-study paper-grade analysis (post-fleet, zero new compute).

Reads only shipped per-run artifacts (run_manifest.json, scout_final_matched.csv)
and the shipped real-section table (results/scout_final_matched.csv).

Truth-model note (amends spec S4 as executed). Arm N permutes ONLY the 21
tested feature names; class labels, geometry, donor sets, regions, and the
null machinery are unchanged. A row (L -> R) of a pseudo-section decomposes
by the class map:

  sigma = 1 (class-preserving): TYPE_MAP[L] == TYPE_MAP[pre(L)] and
    RECV_MAP[R] == RECV_MAP[pre(R)], where pre() is the pre-image of the name
    permutation. The row's content assignment equals the real section's, so
    its statistic is BITWISE the real-section statistic of the pre-image pair
    (same counted sets, same donor map); if the pre-image pair is one of the
    25 tested pairs, its matched-null p/q are bitwise the real-section values
    too (same null seeds and blocks).
  sigma = 0 (mixed-class): the row mixes content across classes; no
    real-section analogue exists.

Consequences:
  - Every certificate (either direction) on every arm is a sound superset
    statement about its own pseudo-section (Propositions 1-2 logic); the
    certificates never certify anything false about the section they score.
  - Content-invariance receipt: whenever a row's sender (resp. receiver)
    slot is class-preserving, the row's sender (resp. receiver) count is
    BITWISE the real section's count of the pre-image gene in that region.
    Verified over all class-preserving slots of all 6 Arm-N runs.
  - "False certification" is therefore truth-model-relative:
    * construction-truth model: 0 by soundness, witnessed by the slot-level
      content-invariance receipt;
    * naive truth model ("every row null by construction"): the certified
      rows are the false certifications, and they are exclusively sigma = 0
      mixed-class rows (inherited real spatial structure, not machinery
      errors). sigma = 1 certified rows would be inherited real-section
      certifications.

Outputs (written to results/synthstudy/):
  v1_slot_identity.csv  per-run slot-level content-invariance tallies (receipt)
  v1_sigma_audit.csv   per Arm-N run x row: sigma, pre-image, verdicts, matches
  v1_fcr_table.csv     per-run + pooled FCR under both truth models
  v2_dose_curve.csv    Arm E: measured injected mass vs injected-pair T_lo
  v3_power_curve.csv   Arm M: achieved mass share vs injected-pair outcome
  synthstudy_analysis.sha256  receipt over all written files
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent
ROOT = SRC.parent
OUT = ROOT / "results"
STUDY = OUT / "synthstudy"

sys.path.insert(0, str(SRC))
import b03_synthstudy as SS  # module-level constants only; no data load

REAL_FINAL = OUT / "scout_final_matched.csv"
N_RUNS = [f"N{i:02d}" for i in range(1, 7)]
E_RUNS = ["E001F10", "E001F40", "E002F10", "E002F40", "E003F10", "E003F40"]
M_RUNS = ["M0801", "M0802", "M0401", "M0402", "M0101", "M0102",
          "M0051", "M0052"]


def load_real():
    real = pd.read_csv(REAL_FINAL)
    assert len(real) == 75, f"real table rows {len(real)}"
    assert real.certified_pos_matched.sum() == 37
    assert real.certified_neg_matched.sum() == 9
    return real


def pi_preimage(pi_seed):
    """Pre-image map of the Arm-N name permutation.

    build_pseudo_section does:
        rng = np.random.default_rng(pi_seed)
        perm = rng.permutation(21)
        remap[tested_codes] = tested_codes[perm]
    i.e. content of original gene TESTED[i] now carries name TESTED[perm[i]].
    pre[j] = i (the gene whose content the name TESTED[j] now holds).
    """
    rng = np.random.default_rng(pi_seed)
    perm = rng.permutation(len(SS.TESTED))
    pre = np.empty(len(SS.TESTED), dtype=np.int64)
    pre[perm] = np.arange(len(SS.TESTED))
    return {SS.TESTED[j]: SS.TESTED[pre[j]] for j in range(len(SS.TESTED))}


def class_of(g):
    """Class of a tested gene in the union class map. The 5 dual-role genes
    (CD274, CD3D, ERBB2, ESR1, PGR) carry the same class in both maps."""
    c = SS.TYPE_MAP.get(g) or SS.RECV_MAP[g]
    if g in SS.TYPE_MAP and g in SS.RECV_MAP:
        assert SS.TYPE_MAP[g] == SS.RECV_MAP[g], f"inconsistent class for {g}"
    return c


def sigma_row(L, R, pre):
    """1 if the row's declared classes are preserved by the permutation.

    Declared classes come from the role maps (TYPE_MAP for the sender slot,
    RECV_MAP for the receiver slot); content classes come from the union map
    applied to the pre-image names."""
    return int(class_of(pre[L]) == SS.TYPE_MAP[L]
               and class_of(pre[R]) == SS.RECV_MAP[R])


def slot_identity_check(real_bounds):
    """Slot-level content-invariance receipt: for any Arm-N row whose sender
    slot is class-preserving (class_of(pre(L)) == TYPE_MAP[L]), the row's nL
    must bitwise equal the real section's sender count of pre(L) in that
    region; likewise for class-preserving receiver slots. The real count maps
    are pair-invariant (verified here), so this is an exact comparison.
    Returns (per-run tallies, list of mismatches)."""
    sender_counts, recv_counts = {}, {}
    for _, r in real_bounds.iterrows():
        sender_counts.setdefault(r.pair.split("->")[0], {})[int(r.region)] = int(r.nL)
        recv_counts.setdefault(r.pair.split("->")[1], {})[int(r.region)] = int(r.nR)
    for a, d in sender_counts.items():
        assert len(d) == 3, ("sender count not pair-invariant", a, d)
    for b, d in recv_counts.items():
        assert len(d) == 3, ("receiver count not pair-invariant", b, d)
    tallies, bad = [], []
    for rid in N_RUNS:
        m = json.loads((STUDY / rid / "run_manifest.json").read_text())
        pre = pi_preimage(m["spec"]["seed"])
        nb = pd.read_csv(STUDY / rid / "scout_bounds.csv")
        t = dict(run_id=rid, s_ok=0, s_tot=0, s_skip=0, r_ok=0, r_tot=0, r_skip=0)
        for _, row in nb.iterrows():
            L, R = row.pair.split("->")
            reg = int(row.region)
            pL, pR = pre[L], pre[R]
            if class_of(pL) == SS.TYPE_MAP[L]:
                if pL in sender_counts:
                    t["s_tot"] += 1
                    if sender_counts[pL][reg] == int(row.nL):
                        t["s_ok"] += 1
                    else:
                        bad.append((rid, "S", row.pair, reg, pL, int(row.nL),
                                    sender_counts[pL][reg]))
                else:
                    t["s_skip"] += 1
            if class_of(pR) == SS.RECV_MAP[R]:
                if pR in recv_counts:
                    t["r_tot"] += 1
                    if recv_counts[pR][reg] == int(row.nR):
                        t["r_ok"] += 1
                    else:
                        bad.append((rid, "R", row.pair, reg, pR, int(row.nR),
                                    recv_counts[pR][reg]))
                else:
                    t["r_skip"] += 1
        tallies.append(t)
    return pd.DataFrame(tallies), bad


def audit_arm_n(real):
    rows = []
    for rid in N_RUNS:
        rd = STUDY / rid
        m = json.loads((rd / "run_manifest.json").read_text())
        assert m["spec"]["arm"] == "N"
        pi_seed = m["spec"]["seed"]
        fin = pd.read_csv(rd / "scout_final_matched.csv")
        assert len(fin) == 75
        pre = pi_preimage(pi_seed)
        pair_set = {f"{a}->{b}" for a, b in SS.ALL_PAIRS}
        for _, r in fin.iterrows():
            L, R = r.pair.split("->")
            sg = sigma_row(L, R, pre)
            preL, preR = pre[L], pre[R]
            pre_pair = f"{preL}->{preR}"
            comparable = pre_pair in pair_set
            tlo_match = thi_match = verd_match = None
            tlo_diff = thi_diff = np.nan
            if comparable:
                rr = real[(real.pair == pre_pair) & (real.region == r.region)]
                assert len(rr) == 1, f"pre-image {pre_pair} r{r.region} lookup"
                rr = rr.iloc[0]
                tlo_diff = abs(float(r.T_lo) - float(rr.T_lo))
                thi_diff = abs(float(r.T_hi) - float(rr.T_hi))
                tlo_match = bool(float(r.T_lo) == float(rr.T_lo))
                thi_match = bool(float(r.T_hi) == float(rr.T_hi))
                verd_match = bool(
                    bool(r.certified_pos_matched) == bool(rr.certified_pos_matched)
                    and bool(r.certified_neg_matched) == bool(rr.certified_neg_matched))
            rows.append(dict(
                run_id=rid, seed=pi_seed, pair=r.pair, region=int(r.region),
                sigma=sg,
                slotL_preserving=int(class_of(preL) == SS.TYPE_MAP[L]),
                slotR_preserving=int(class_of(preR) == SS.RECV_MAP[R]),
                preimage_pair=pre_pair, preimage_in_panel=comparable,
                T_lo=float(r.T_lo), T_hi=float(r.T_hi),
                cert_pos=bool(r.certified_pos_matched),
                cert_neg=bool(r.certified_neg_matched),
                pre_cert_pos=bool(rr.certified_pos_matched) if comparable else None,
                pre_cert_neg=bool(rr.certified_neg_matched) if comparable else None,
                T_lo_exact=tlo_match, T_hi_exact=thi_match,
                T_lo_diff=tlo_diff, T_hi_diff=thi_diff,
                verdict_match=verd_match))
    return pd.DataFrame(rows)


def fcr_table(aud):
    out = []
    for rid in N_RUNS:
        a = aud[aud.run_id == rid]
        n = len(a)
        tot = int(a.cert_pos.sum() + a.cert_neg.sum())
        s1 = a[a.sigma == 1]
        s0 = a[a.sigma == 0]
        comp = a[a.preimage_in_panel]
        out.append(dict(
            run_id=rid, n_rows=n,
            n_cert_total=tot,
            n_cert_pos=int(a.cert_pos.sum()),
            n_cert_neg=int(a.cert_neg.sum()),
            fcr_naive=tot / n,
            n_cert_sigma1=int(s1.cert_pos.sum() + s1.cert_neg.sum()),
            n_cert_sigma0=int(s0.cert_pos.sum() + s0.cert_neg.sum()),
            fcr_sigma0=(int(s0.cert_pos.sum() + s0.cert_neg.sum())) / max(len(s0), 1),
            n_comparable=len(comp),
            n_stat_exact=int(comp.T_lo_exact.sum() + comp.T_hi_exact.sum()),
            n_stat_compared=2 * len(comp),
            n_verdict_match=int(comp.verdict_match.sum()),
        ))
    df = pd.DataFrame(out)
    n_sigma1_total = int(sum((aud[aud.run_id == r].sigma == 1).sum()
                             for r in N_RUNS))
    n_sigma0_total = int(df.n_rows.sum()) - n_sigma1_total
    pooled = dict(run_id="POOLED", n_rows=int(df.n_rows.sum()),
                  n_cert_total=int(df.n_cert_total.sum()),
                  n_cert_pos=int(df.n_cert_pos.sum()),
                  n_cert_neg=int(df.n_cert_neg.sum()),
                  fcr_naive=df.n_cert_total.sum() / df.n_rows.sum(),
                  n_cert_sigma1=int(df.n_cert_sigma1.sum()),
                  n_cert_sigma0=int(df.n_cert_sigma0.sum()),
                  # sigma-0 denominator: total rows minus sigma-1 rows
                  fcr_sigma0=int(df.n_cert_sigma0.sum()) / n_sigma0_total,
                  n_sigma1_rows=n_sigma1_total,
                  n_comparable=int(df.n_comparable.sum()),
                  n_stat_exact=int(df.n_stat_exact.sum()),
                  n_stat_compared=int(df.n_stat_compared.sum()),
                  n_verdict_match=int(df.n_verdict_match.sum()))
    df = pd.concat([df, pd.DataFrame([pooled])], ignore_index=True)
    df["n_sigma1_rows"] = df["n_sigma1_rows"].fillna(0).astype(int)
    return df


def dose_and_power():
    dose = []
    for rid in E_RUNS:
        rd = STUDY / rid
        m = json.loads((rd / "run_manifest.json").read_text())
        assert m["spec"]["arm"] == "E"
        inj = m["info"]["injection"]
        mass = int(inj["swapped_L"]) + int(inj["swapped_R"])
        fin = pd.read_csv(rd / "scout_final_matched.csv")
        ij = fin[fin.pair == f"{SS.INJ_PAIR[0]}->{SS.INJ_PAIR[1]}"]
        assert len(ij) == 3, f"{rid} injected-pair rows {len(ij)}"
        for _, r in ij.iterrows():
            dose.append(dict(
                run_id=rid, arm="E", seed=m["spec"]["seed"], f=inj["f"],
                injected_mass=mass, swapped_L=int(inj["swapped_L"]),
                swapped_R=int(inj["swapped_R"]), region=int(r.region),
                T_lo=float(r.T_lo), T_hi=float(r.T_hi),
                q_pos=float(r.q_pos),
                cert_pos=bool(r.certified_pos_matched),
                cert_neg=bool(r.certified_neg_matched)))
    dose = pd.DataFrame(dose).sort_values(["injected_mass", "region"])
    power = []
    for rid in M_RUNS:
        rd = STUDY / rid
        m = json.loads((rd / "run_manifest.json").read_text())
        assert m["spec"]["arm"] == "M"
        inj = m["info"]["injection"]
        mass = int(inj["swapped_L"]) + int(inj["swapped_R"])
        fin = pd.read_csv(rd / "scout_final_matched.csv")
        ij = fin[fin.pair == f"{SS.INJ_PAIR[0]}->{SS.INJ_PAIR[1]}"]
        assert len(ij) == 3, f"{rid} injected-pair rows {len(ij)}"
        power.append(dict(
            run_id=rid, arm="M", target_pct=m["spec"]["target_pct"],
            seed=m["spec"]["seed"], achieved_share=m["diag"]["tested_share"],
            injected_mass=mass, swapped_L=int(inj["swapped_L"]),
            swapped_R=int(inj["swapped_R"]),
            T_lo_med=float(ij.T_lo.median()),
            T_lo_min=float(ij.T_lo.min()), T_lo_max=float(ij.T_lo.max()),
            q_pos_med=float(ij.q_pos.median()),
            n_cert_pos_pair=int(ij.certified_pos_matched.sum()),
            n_cert_neg_pair=int(ij.certified_neg_matched.sum()),
            run_cert_pos=m["n_cert_pos"], run_cert_neg=m["n_cert_neg"]))
    power = pd.DataFrame(power).sort_values(["target_pct", "seed"])
    return dose, power


def sha256_of(p):
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def main():
    real = load_real()
    print("real section: 75 rows,",
          int(real.certified_pos_matched.sum()), "cert-pos /",
          int(real.certified_neg_matched.sum()), "cert-neg")

    tallies, bad = slot_identity_check(pd.read_csv(OUT / "scout_bounds.csv"))
    tallies.to_csv(STUDY / "v1_slot_identity.csv", index=False)
    assert not bad, f"content-invariance violations: {bad[:5]}"
    print("\nV1 slot-level content-invariance receipt (class-preserving slots,"
          " bitwise vs real section):")
    print(tallies.to_string(index=False))
    print(f"  TOTAL exact: {int(tallies.s_ok.sum() + tallies.r_ok.sum())}/"
          f"{int(tallies.s_tot.sum() + tallies.r_tot.sum())}; "
          "mismatches: 0 (asserted)")

    aud = audit_arm_n(real)
    aud.to_csv(STUDY / "v1_sigma_audit.csv", index=False)

    comp = aud[aud.preimage_in_panel]
    s1_comp = comp[comp.sigma == 1]
    n_exact = int(s1_comp.T_lo_exact.sum() + s1_comp.T_hi_exact.sum())
    print(f"\nV1 sigma audit: {len(aud)} rows over {len(N_RUNS)} replicates")
    print(f"  sigma=1 rows: {(aud.sigma == 1).sum()} "
          f"(pre-image in 25-pair panel: {len(s1_comp)}; "
          f"bitwise statistic matches: {n_exact}/"
          f"{2 * len(s1_comp)} among those)")
    print(f"  sigma=0 rows: {(aud.sigma == 0).sum()} "
          f"(mixed-class content; no real-section analogue - the "
          f"{len(comp) - len(s1_comp)} panel-matching rows among them are NOT "
          f"expected to match real statistics)")
    s1c = aud[(aud.sigma == 1) & (aud.cert_pos | aud.cert_neg)]
    print(f"  sigma=1 certified rows: {len(s1c)} "
          f"(pos {int(s1c.cert_pos.sum())} / neg {int(s1c.cert_neg.sum())})")
    s0c = aud[(aud.sigma == 0) & (aud.cert_pos | aud.cert_neg)]
    print(f"  sigma=0 certified rows: {len(s0c)} "
          f"(pos {int(s0c.cert_pos.sum())} / neg {int(s0c.cert_neg.sum())})")
    mech = s0c.preimage_pair.str.contains("ERBB2").sum()
    print(f"  of the certified sigma=0 rows, {mech}/{len(s0c)} carry ERBB2 "
          f"content in the pre-image (the dominant tested-gene mass)")
    print("  per-run sigma=0 certified counts:")
    for rid in N_RUNS:
        s0r = s0c[s0c.run_id == rid]
        print(f"    {rid}: {len(s0r)} (pos {int(s0r.cert_pos.sum())} / "
              f"neg {int(s0r.cert_neg.sum())})")

    fcr = fcr_table(aud)
    fcr.to_csv(STUDY / "v1_fcr_table.csv", index=False)
    print("\nV1 FCR (per run and pooled):")
    print(fcr.to_string(index=False))

    dose, power = dose_and_power()
    dose.to_csv(STUDY / "v2_dose_curve.csv", index=False)
    power.to_csv(STUDY / "v3_power_curve.csv", index=False)

    print("\nV2 dose curve (Arm E, injected pair KRT8->ERBB2, 3 regions/run, "
          "ambient share 13.84%):")
    g = dose.groupby(["run_id", "injected_mass"], as_index=False).agg(
        T_lo_min=("T_lo", "min"), T_lo_med=("T_lo", "median"),
        T_lo_max=("T_lo", "max"), q_pos_med=("q_pos", "median"),
        n_cert_pos=("cert_pos", "sum"))
    print(g.to_string(index=False))
    neg_runs = g[g.T_lo_max <= 0]
    pos_runs = g[g.T_lo_min > 0]
    print(f"  all-negative up to mass {neg_runs.injected_mass.max() if len(neg_runs) else None}; "
          f"all-positive from mass {pos_runs.injected_mass.min() if len(pos_runs) else None}")

    print("\nV3 power curve (Arm M):")
    print(power.to_string(index=False))

    files = ["v1_sigma_audit.csv", "v1_slot_identity.csv", "v1_fcr_table.csv",
             "v2_dose_curve.csv", "v3_power_curve.csv"]
    receipt = "\n".join(f"{sha256_of(STUDY / f)}  {f}" for f in files) + "\n"
    (STUDY / "synthstudy_analysis.sha256").write_text(receipt)
    print("\nreceipt written: synthstudy_analysis.sha256")


if __name__ == "__main__":
    main()
