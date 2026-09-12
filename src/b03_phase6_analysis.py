"""B03 Phase 6 analysis: b10k checkpoints -> p/q CSVs + verdict-stability check.

Registered conventions (B03_PHASE6_PREREG.md; identical to sealed runs):
  p_pos = (1 + #{T_lo^(b) >= T_lo_obs}) / (B+1);  p_neg = (1 + #{T_hi^(b) <= T_hi_obs}) / (B+1)
  BH within region over the pair family, q <= 0.10.
  Verdict: certified_pos = T_lo_obs > 0 AND q_pos <= 0.10;
           certified_neg = T_hi_obs < 0 AND q_neg <= 0.10.

Ground-truth columns per tissue (pre-flight-validated 2026-09-08):
  lung:   scout_final_lung.csv, p_pos_ref/p_neg_ref/q_pos_ref/q_neg_ref,
          verdicts cert_pos_final/cert_neg_final (null_robust_matched.csv is an
          intermediate and is NOT the sealed source; pre-flight caught this).
  breast: scout_final_matched.csv, p_pos_raw/p_neg_raw/q_pos/q_neg,
          verdicts certified_pos_matched/certified_neg_matched.

Stability check (registered): the certified row set at B=10,000 must match the
sealed B=1,000 certified sets EXACTLY (both tissues). Any flip is a registered
anomaly: print loudly, quarantine outputs, exit 2.
"""
import json
import os

import numpy as np
import pandas as pd

RES = "B03_project/results"
LUN = "xenium_lung/crop/data"
B = 10000

TISSUES = {
    "lung": dict(
        out_dir=f"{LUN}/b10k",
        bounds_csv=f"{LUN}/scout_bounds.csv",
        sealed_csv=f"{LUN}/scout_final_lung.csv",
        p_cols=("p_pos_ref", "p_neg_ref"),
        q_cols=("q_pos_ref", "q_neg_ref"),
        cert_cols=("cert_pos_final", "cert_neg_final"),
    ),
    "breast": dict(
        out_dir=f"{RES}/b10k_breast",
        bounds_csv=f"{RES}/scout_bounds.csv",
        sealed_csv=f"{RES}/scout_final_matched.csv",
        p_cols=("p_pos_raw", "p_neg_raw"),
        q_cols=("q_pos", "q_neg"),
        cert_cols=("certified_pos_matched", "certified_neg_matched"),
    ),
}


def bh_qvalues(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def tissue(tag, cfg):
    ck = np.load(f"{cfg['out_dir']}/null_v3_checkpoint.npz")
    Tlo, Thi, done = ck["Tlo"], ck["Thi"], ck["done"]
    n_rows = Tlo.shape[1]
    order = np.argsort(done)
    assert len(order) == B, f"{tag}: expected {B} perms, got {len(order)}"
    Tlo, Thi = Tlo[order], Thi[order]  # rows now in perm-id order 0..B-1

    bounds = pd.read_csv(cfg["bounds_csv"])
    final = pd.read_csv(cfg["sealed_csv"])
    key = ["pair", "region"]
    assert len(bounds) == n_rows and len(final) == n_rows, f"{tag}: row-count mismatch"

    # v3 iterates scout_bounds.csv rows in file order; that order is
    # authoritative for the checkpoint's columns and must not be re-sorted.
    b = bounds[key + ["T_lo", "T_hi"]].rename(columns={"T_lo": "T_lo_b", "T_hi": "T_hi_b"})
    f = final[key + list(cfg["cert_cols"])].rename(
        columns={cfg["cert_cols"][0]: "cert_pos_sealed", cfg["cert_cols"][1]: "cert_neg_sealed"})
    merged = b.merge(f, on=key, how="left", validate="one_to_one")
    assert len(merged) == n_rows
    assert merged["cert_pos_sealed"].notna().all(), f"{tag}: sealed final missing rows"

    p_pos = (1 + (Tlo >= merged["T_lo_b"].to_numpy()[None, :]).sum(axis=0)) / (B + 1)
    p_neg = (1 + (Thi <= merged["T_hi_b"].to_numpy()[None, :]).sum(axis=0)) / (B + 1)
    q_pos, q_neg = bh_qvalues(p_pos), bh_qvalues(p_neg)

    out = merged[key].copy()
    out["T_lo"] = merged["T_lo_b"]
    out["T_hi"] = merged["T_hi_b"]
    out["p_pos_b10k"] = p_pos
    out["p_neg_b10k"] = p_neg
    out["q_pos_b10k"] = q_pos
    out["q_neg_b10k"] = q_neg
    out["cert_pos_b10k"] = (out["T_lo"] > 0) & (out["q_pos_b10k"] <= 0.10)
    out["cert_neg_b10k"] = (out["T_hi"] < 0) & (out["q_neg_b10k"] <= 0.10)

    # ---- verdict-stability check vs sealed B=1k ----
    sealed_pos = merged["cert_pos_sealed"].astype(bool).to_numpy()
    sealed_neg = merged["cert_neg_sealed"].astype(bool).to_numpy()
    new_pos = out["cert_pos_b10k"].to_numpy(bool)
    new_neg = out["cert_neg_b10k"].to_numpy(bool)
    flip_pos = [(merged["pair"].iloc[i], int(merged["region"].iloc[i]))
                for i in np.flatnonzero(new_pos != sealed_pos)]
    flip_neg = [(merged["pair"].iloc[i], int(merged["region"].iloc[i]))
                for i in np.flatnonzero(new_neg != sealed_neg)]

    summary = {
        "tissue": tag,
        "B": B,
        "n_rows": n_rows,
        "min_p_pos": float(p_pos.min()),
        "min_p_neg": float(p_neg.min()),
        "n_cert_pos_sealed": int(sealed_pos.sum()),
        "n_cert_neg_sealed": int(sealed_neg.sum()),
        "n_cert_pos_b10k": int(new_pos.sum()),
        "n_cert_neg_b10k": int(new_neg.sum()),
        "flips_pos": flip_pos,
        "flips_neg": flip_neg,
        "verdict_flips_vs_b1000": len(flip_pos) + len(flip_neg),
    }
    os.makedirs(f"{cfg['out_dir']}/analysis", exist_ok=True)
    out.to_csv(f"{cfg['out_dir']}/analysis/b10k_null.csv", index=False)
    return summary


def main():
    results = {}
    flips_total = 0
    for tag, cfg in TISSUES.items():
        if not os.path.exists(f"{cfg['out_dir']}/null_v3_checkpoint.npz"):
            print(f"[skip] {tag}: checkpoint not present yet")
            continue
        s = tissue(tag, cfg)
        results[tag] = s
        flips_total += s["verdict_flips_vs_b1000"]
        print(f"{tag}: pos {s['n_cert_pos_sealed']}->{s['n_cert_pos_b10k']}, "
              f"neg {s['n_cert_neg_sealed']}->{s['n_cert_neg_b10k']}, "
              f"flips={s['verdict_flips_vs_b1000']} {s['flips_pos'] + s['flips_neg']}")

    os.makedirs(f"{RES}/phase6", exist_ok=True)
    with open(f"{RES}/phase6/phase6_summary.json", "w") as fh:
        json.dump(results, fh, indent=1)

    if flips_total > 0:
        print(f"REGISTERED ANOMALY: {flips_total} verdict flips vs B=1000. "
              f"Outputs quarantined; code audit required before any use.", flush=True)
        raise SystemExit(2)
    print("PHASE 6 ANALYSIS COMPLETE: verdicts stable, p-floor resolved.", flush=True)


if __name__ == "__main__":
    main()
