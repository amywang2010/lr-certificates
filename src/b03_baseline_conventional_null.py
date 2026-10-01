"""B03 conventional point-null (breast): the field's standard first pass, computed
fresh (registered amendment of the conventional protocol).

Why recomputed: the quarantined first-generation null artifact is degenerate by
construction (q_pos = 1.0 on all rows; worst-case bounds trivially dominate point
statistics), which is exactly why it was killed for OUR test. As a description of
the conventional pipeline it must be a real run of the conventional test.

Registered machinery reused verbatim from the sealed matched null
(src/b03_robust_null_v3.py): within-(region x transcript-decile) label permutation
blocks (same construction, same registration order), per-permutation seeds
rng(SEED_BASE + p) with SEED_BASE = 20260905 + 100003, B = 1,000, frozen cell-count
denominators (nl_cells, nr_cells from scout_bounds.csv), epsilon = 0.5, BH within
region. The single difference is the tested object: the point statistic
    T(A) = log2(N_L/nLc + eps) + log2(N_R/nRc + eps)
under permuted labels, tested against T(A0) — the conventional co-expression
permutation test (CellPhoneDB/CellChat-family null), no worst-case layer.

Cross-implementation check: the recomputed T(A0) must match the sealed scout T0
column (max |delta| < 1e-9) before any p-value is produced.

Output: results/baseline_comparison/conventional_null_breast.csv
        (pair, region, T0, p_pos, p_neg, q_pos, q_neg — add-one p-values).
"""
from __future__ import annotations

import json
import math
import os
import pickle
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_certify import ALL_PAIRS  # registered pair list, sealed
from b03_robust_null_v3 import TYPE_MAP, RECV_MAP  # registered sender/receiver maps

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")
OUTDIR = os.path.join(RES, "baseline_comparison")
EPS = 0.5
B_PERM = 1000
SEED_BASE = 20260905 + 100003

tested = sorted({g for p in ALL_PAIRS for g in p})
gidx = {g: i for i, g in enumerate(tested)}


def bh_within_region(df: pd.DataFrame, pcol: str) -> pd.Series:
    out = []
    for region_id, gdf in df.groupby("region"):
        m = len(gdf)
        order = np.argsort(gdf[pcol].to_numpy())
        p_sorted = gdf[pcol].to_numpy()[order]
        q = np.minimum.accumulate((p_sorted * m / np.arange(1, m + 1))[::-1])[::-1]
        qq = np.minimum(q, 1.0)
        res_q = np.empty(m)
        res_q[order] = qq
        out.append(pd.Series(res_q, index=gdf.index))
    return pd.concat(out).sort_index()


def main() -> None:
    t0 = time.time()
    tx = pd.read_parquet(os.path.join(DATA, "tx.parquet"))
    cells_meta = pd.read_parquet(os.path.join(DATA, "cells_meta.parquet"))
    labels = np.load(os.path.join(DATA, "labels.npy"), allow_pickle=True)
    region_of_cell = np.load(os.path.join(DATA, "regions.npy"))
    res = pd.read_csv(os.path.join(RES, "scout_bounds.csv"))

    n_cells = len(cells_meta)
    assert n_cells == len(labels) == len(region_of_cell), "length misalignment"
    vc = pd.Series(labels).value_counts()
    assert vc.get("Breast cancer") == 60459 and vc.get("T cells") == 20006, \
        f"label counts differ from scout run: {vc.head(3).to_dict()}"

    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()
    Gmat = np.zeros((n_cells, len(tested)), dtype=np.int64)
    for g in tested:
        sel_g = (feat == g) & (cell >= 0)
        Gmat[:, gidx[g]] = np.bincount(cell[sel_g], minlength=n_cells)

    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    n_types = len(types_needed)
    type_id = {t: i for i, t in enumerate(types_needed)}
    lab0 = np.array([type_id.get(l, -1) for l in labels], dtype=np.int64)

    # registered blocks: (region, decile of total transcripts) — verbatim from v3
    base_tc = cells_meta.transcript_counts.to_numpy(np.int64)
    dec = np.zeros(n_cells, dtype=np.int64)
    for rr in range(3):
        idx_r = np.flatnonzero(region_of_cell == rr)
        q = np.quantile(base_tc[idx_r], np.linspace(0, 1, 11))
        dec[idx_r] = np.clip(np.digitize(base_tc[idx_r], q[1:-1]), 0, 9)
    blocks = []
    for rr in range(3):
        for d in range(10):
            idx = np.flatnonzero((region_of_cell == rr) & (dec == d))
            if len(idx) > 1:
                blocks.append(idx)
    print(f"blocks: {len(blocks)}; sizes {min(map(len, blocks))}-{max(map(len, blocks))}",
          flush=True)

    # keys: (pair, region) -> target type slots; denominators frozen from scout_bounds
    keys = []
    for _, r in res.iterrows():
        Lg, Rg = r.pair.split("->")
        region_id = int(r.region)
        keys.append((r.pair, region_id,
                     type_id[TYPE_MAP[Lg]], type_id[RECV_MAP[Rg]],
                     int(r.nl_cells), int(r.nr_cells)))
    assert len(keys) == 75

    # per (region, type) cell index lists under A0, and observed T0
    def N_of(lab, gi, target, region_id):
        return int(Gmat[(region_of_cell == region_id) & (lab == target), gi].sum())

    T0 = np.empty(len(keys))
    for k, (pair, region_id, tS, tR, nLc, nRc) in enumerate(keys):
        Lg, Rg = pair.split("->")
        NL = N_of(lab0, gidx[Lg], tS, region_id)
        NR = N_of(lab0, gidx[Rg], tR, region_id)
        T0[k] = math.log2(NL / nLc + EPS) + math.log2(NR / nRc + EPS)
    scout = pd.read_csv(os.path.join(RES, "scout_final_matched.csv"))
    mine = pd.DataFrame({"pair": [k[0] for k in keys],
                         "region": [k[1] for k in keys], "T0_mine": T0})
    mm = mine.merge(scout[["pair", "region", "T0"]], on=["pair", "region"], how="inner")
    assert len(mm) == 75, f"key merge lost rows: {len(mm)}"
    dmax = float(np.max(np.abs(mm["T0_mine"] - mm["T0"])))
    assert dmax < 1e-9, f"T0 cross-implementation check failed (max |d| = {dmax:.3e})"
    print(f"[{time.time()-t0:.0f}s] T0 verified vs sealed scout (max |d| < 1e-9)", flush=True)

    # precompute region masks once; per permutation, group cells by (region, perm label)
    region_masks = [np.flatnonzero(region_of_cell == rr) for rr in range(3)]

    def permuted_counts(lab_perm):
        """Return dict (region, type) -> Gmat row-sum vector (len n_genes)."""
        out = {}
        for rr in range(3):
            idx = region_masks[rr]
            labs_r = lab_perm[idx]
            order = np.argsort(labs_r, kind="stable")
            srt = Gmat[idx[order]]
            lab_sorted = labs_r[order]
            boundaries = np.flatnonzero(np.diff(lab_sorted)) + 1
            starts = np.r_[0, boundaries]
            sums = np.add.reduceat(srt, starts, axis=0)
            uniq = lab_sorted[starts]
            for row, t in zip(sums, uniq):
                if t >= 0:
                    out[(rr, int(t))] = row
        return out

    rng0 = np.random.default_rng(SEED_BASE + 0)
    lab_perm = lab0.copy()
    for idx in blocks:
        pm = rng0.permutation(len(idx))
        lab_perm[idx] = lab_perm[idx[pm]]
    counts0 = permuted_counts(lab_perm)

    T = np.empty((B_PERM, len(keys)))
    for p in range(B_PERM):
        rng = np.random.default_rng(SEED_BASE + p)
        lab_perm = lab0.copy()
        for idx in blocks:
            pm = rng.permutation(len(idx))
            lab_perm[idx] = lab_perm[idx[pm]]
        cnt = permuted_counts(lab_perm)
        for k, (pair, region_id, tS, tR, nLc, nRc) in enumerate(keys):
            Lg, Rg = pair.split("->")
            NL = cnt.get((region_id, tS), np.zeros(len(tested), dtype=np.int64))[gidx[Lg]]
            NR = cnt.get((region_id, tR), np.zeros(len(tested), dtype=np.int64))[gidx[Rg]]
            T[p, k] = math.log2(NL / nLc + EPS) + math.log2(NR / nRc + EPS)
        if (p + 1) % 200 == 0:
            print(f"[{time.time()-t0:.0f}s] perm {p+1}/{B_PERM}", flush=True)

    p_pos = (1 + (T >= T0[None, :]).sum(axis=0)) / (B_PERM + 1)
    p_neg = (1 + (T <= T0[None, :]).sum(axis=0)) / (B_PERM + 1)

    df = pd.DataFrame(keys, columns=["pair", "region", "tS", "tR", "nLc", "nRc"])
    df["T0"] = T0
    df["p_pos"] = p_pos
    df["p_neg"] = p_neg
    df["q_pos"] = bh_within_region(df, "p_pos")
    df["q_neg"] = bh_within_region(df, "p_neg")
    df = df.sort_values(["pair", "region"])[["pair", "region", "T0", "p_pos", "p_neg",
                                             "q_pos", "q_neg"]]
    os.makedirs(OUTDIR, exist_ok=True)
    df.to_csv(os.path.join(OUTDIR, "conventional_null_breast.csv"), index=False)
    print(f"done in {time.time()-t0:.0f}s -> conventional_null_breast.csv", flush=True)


if __name__ == "__main__":
    main()
