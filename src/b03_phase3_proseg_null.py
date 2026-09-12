"""B03 Phase 3: Proseg-side matched permutation null — per B03_PHASE3_PREREG.md
(FROZEN 2026-09-07 09:40, before any Phase-3 compute; amendment record).

Question closed: are the 46 matched-certified directional signs ALSO significant
under the independent segmentation's own null (within-region decile-matched label
permutation of Proseg labels)? Gate: >= 85% of certified rows significant at
within-region BH q <= 0.10 in the certified direction.

Machinery reuse (prereg):
  - EPS, t_log2, ALL_PAIRS, TYPE_MAP, RECV_MAP: IMPORTED from b03_scout.py.
  - bh_within_region: IMPORTED from b03_robust_null_v3.py (identical BH code the
    scout used; module import is side-effect-free, main is guarded).
  - Label transfer + region mapping: identical algorithm to b03_phase2_analysis.py
    (verbatim logic; correctness enforced by Assertion A below, which pins the
    reproduction to the phase-2 outputs on all 75 rows).

Hard assertions (abort on failure; prereg section "Hard assertions"):
  A. Label/region reproduction: re-derived per-(pair,region) sender/receiver cell
     counts nL/nR must equal phase2_proseg_coverage.csv EXACTLY for all 75 rows.
  B. Observed-statistic crosscheck: T_obs recomputed here equals the coverage CSV
     `T` column within 1e-9 for all 75 rows.
  C. Permutation invariance: per-region per-type cell counts preserved under
     permutation, checked on 10 pre-specified permutations.
  D. Pre-flight determinism: perm 0 computed twice (parent + worker) — bitwise
     equality required before the sweep starts.

Null scheme (prereg, identical structure to scout amendment record/amendment record):
  blocks = (region) x (decile of per-cell TOTAL foreground Proseg counts over ALL
  named genes; per-entry rint then sum — the pipeline's count convention);
  np.quantile linspace(0,1,11), digitize on q[1:-1], clip 0..9; blocks with >=2
  cells; labels permuted within blocks in fixed (region, decile) order;
  B = 1000; per-perm seed default_rng(20260907 + 100003 + p).

p-values: raw counts reported; inference on add-one p = (count+1)/(B+1); BH
within region applied to add-one p (documented here; scout's watchdog used the
same effective convention).

Smoke mode (B03_P3_SMOKE=1): B=20, outputs to results/_phase3_dryrun/ — validates
the full path end-to-end BEFORE the real run, without touching production paths.
"""

import gzip
import json
import math
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.spatial import cKDTree
import multiprocessing as mp

sys.path.insert(0, "B03_project/src")
from b03_scout import ALL_PAIRS, TYPE_MAP, RECV_MAP, EPS, t_log2
from b03_robust_null_v3 import bh_within_region

DATA = "B03_project/data"
RES = "B03_project/results"
# Tissue-2 reuse (amendment record): paths and expected row count env-overridable.
DATA = os.environ.get("B03_P3_DATA", DATA)
RES = os.environ.get("B03_P3_RES", RES)
N_ROWS = int(os.environ.get("B03_P3_ROWS", "75"))
SMOKE = os.environ.get("B03_P3_SMOKE") == "1"
SEED = 20260907
SEED_BASE = SEED + 100003
B_PERM = 20 if SMOKE else 1000
N_WORKERS = 4
CKPT_EVERY = 50
OUT = f"{RES}/_phase3_dryrun" if SMOKE else RES
# pre-specified permutations for invariance Assertion C
CHECK_PERMS = {0, 50, 123, 250, 375, 500, 625, 750, 875, 999} & set(range(B_PERM))

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:.0f}s] {msg}", flush=True)


# ---------------------------------------------------------------------------
# shared worker state
# ---------------------------------------------------------------------------
_W = {}


def worker_init(Gm, cell_slot, blocks, obs, n_slots):
    _W["Gm"] = Gm
    _W["cell_slot"] = cell_slot
    _W["blocks"] = blocks
    _W["obs"] = obs  # list of (giL, giR, tS, tR, nL, nR, T_obs)
    _W["n_slots"] = n_slots


def permute_and_score(p):
    cell_slot = _W["cell_slot"]
    Gm = _W["Gm"]
    obs = _W["obs"]
    rng = np.random.default_rng(SEED_BASE + p)
    lab = cell_slot.copy()
    for idx in _W["blocks"]:
        pm = rng.permutation(len(idx))
        lab[idx] = lab[idx[pm]]
    if p in CHECK_PERMS:  # Assertion C: invariance of per-region type counts
        n_slots = _W["n_slots"]
        for r in range(3):
            m_r = (cell_slot // n_slots) == r
            b_before = np.bincount(cell_slot[m_r], minlength=n_slots * 3)
            b_after = np.bincount(lab[m_r], minlength=n_slots * 3)
            assert (b_before == b_after).all(), f"invariance violated at perm {p} region {r}"
    T = np.empty(len(obs))
    masks = {}
    for k, (giL, giR, tS, tR, nL, nR, _To) in enumerate(obs):
        mS = masks.get(tS)
        if mS is None:
            mS = lab == tS
            masks[tS] = mS
        mR = masks.get(tR)
        if mR is None:
            mR = lab == tR
            masks[tR] = mR
        NL = int(Gm[mS, giL].sum())
        NR = int(Gm[mR, giR].sum())
        T[k] = math.log2(NL / nL + EPS) + math.log2(NR / nR + EPS)
    return p, T


def main():
    # ---------- frozen vendor-side artifacts (identical to phase 2) ----------
    vm = np.load(f"{DATA}/vendor_markers.npz", allow_pickle=True)
    centroids = vm["centroids"].astype(np.float64)
    type_names = [str(t) for t in vm["type_names"]]
    marker_genes = [str(g) for g in vm["marker_genes"]]
    scale = float(vm["scale"])
    vendor_xy = vm["vendor_centroids_xy"]
    regions_vendor = vm["regions_row"].astype(int)

    cov = pd.read_csv(f"{RES}/phase2_proseg_coverage.csv")
    assert len(cov) == N_ROWS

    tested = sorted({g for p in ALL_PAIRS for g in p})
    gidx = {g: i for i, g in enumerate(tested)}

    # ---------- proseg outputs (schema verified from source in phase 2) ----------
    with gzip.open(f"{DATA}/proseg_out/counts", "rb") as f:
        Coo = mmread(f).tocoo()
    cells_p = pd.read_parquet(f"{DATA}/proseg_out/cells")
    genes_p = pd.read_parquet(f"{DATA}/proseg_out/genes")
    log(f"mtx {Coo.shape} nnz={Coo.nnz:,}; cells {cells_p.shape}; genes {genes_p.shape}")
    assert Coo.shape[0] == len(cells_p)
    assert Coo.shape[1] in (len(genes_p), len(genes_p) - 1), "unexpected gene layout"

    is_bg = cells_p["original_cell_id"].astype(str) == "0"
    n_bg = int(is_bg.sum())
    if n_bg:
        keep = (~is_bg).to_numpy().nonzero()[0]
        Coo = Coo.tocsr()[keep].tocoo()
        cells_p = cells_p[~is_bg].reset_index(drop=True)
    n_pc = len(cells_p)
    log(f"background rows dropped: {n_bg} (phase-2 observed 0); universe {n_pc}")

    gene_names = genes_p["gene"].astype(str).to_numpy()
    if Coo.shape[1] == len(gene_names) - 1:
        gene_names_eff = gene_names[:-1]
    else:
        gene_names_eff = gene_names
    need = marker_genes + [g for g in tested if g not in set(marker_genes)]
    gene_pos = {g: i for i, g in enumerate(gene_names_eff)}
    missing = [g for g in need if g not in gene_pos]
    assert not missing, f"genes missing: {missing[:5]}"
    cols = np.array([gene_pos[g] for g in need])
    Sub = Coo.tocsr()[:, cols].tocoo()
    n_need = len(need)
    Mall = np.zeros((n_pc, n_need), dtype=np.int32)
    np.add.at(Mall, (Sub.row.astype(np.int64), Sub.col.astype(np.int64)),
              np.rint(Sub.data).astype(np.int64))
    need_pos = {g: i for i, g in enumerate(need)}
    M = Mall[:, [need_pos[g] for g in marker_genes]].astype(np.int32)
    Gm = Mall[:, [need_pos[g] for g in tested]].astype(np.int32)
    del Mall, Sub, Coo

    # decile covariate: per-cell total over ALL named genes (prereg), per-entry rint
    Csr_full = None
    with gzip.open(f"{DATA}/proseg_out/counts", "rb") as f:
        Full = mmread(f).tocoo()
    if n_bg:
        Full = Full.tocsr()[(~is_bg).to_numpy().nonzero()[0]].tocoo()
    tot = np.zeros(n_pc, dtype=np.int64)
    np.add.at(tot, Full.row.astype(np.int64), np.rint(Full.data).astype(np.int64))
    del Full
    log(f"decile covariate computed (per-cell total, median {np.median(tot):.0f})")

    # ---------- label transfer (verbatim phase-2 algorithm) ----------
    Ctot = M.sum(axis=1).astype(np.int64)
    has_marker = Ctot > 0
    Xp = np.log1p((M / np.maximum(Ctot, 1)[:, None]) * scale).astype(np.float64)
    Xn = Xp / np.maximum(np.linalg.norm(Xp, axis=1, keepdims=True), 1e-12)
    Cn = centroids / np.maximum(np.linalg.norm(centroids, axis=1, keepdims=True), 1e-12)
    labels = np.full(n_pc, "Unknown", dtype=object)
    sim = Xn[has_marker] @ Cn.T
    labels[has_marker] = np.array(type_names, dtype=object)[sim.argmax(axis=1)]
    log(f"labels reproduced: {pd.Series(labels).value_counts().head(4).to_dict()}")

    # ---------- region mapping (verbatim phase-2 rule 1) ----------
    pc_xy = cells_p[["centroid_x", "centroid_y"]].to_numpy(np.float64)
    tree = cKDTree(vendor_xy)
    _dist, nn = tree.query(pc_xy, k=1)
    regions = regions_vendor[nn]

    # ---------- Assertion A: reproduction pinned to phase-2 outputs ----------
    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    n_types = len(types_needed)
    SHIFT = 1
    UNKNOWN_SLOT = 0
    type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
    n_slots = n_types + SHIFT
    cell_slot = regions * n_slots + np.array(
        [type_shift.get(l, UNKNOWN_SLOT) for l in labels], dtype=np.int64)

    obs = []
    for _, r in cov.iterrows():
        Lg, Rg = r.pair.split("->")
        region_id = int(r.region)
        tS = region_id * n_slots + type_shift[TYPE_MAP[Lg]]
        tR = region_id * n_slots + type_shift[RECV_MAP[Rg]]
        mS = cell_slot == tS
        mR = cell_slot == tR
        nL, nR = int(mS.sum()), int(mR.sum())
        assert nL == int(r.nL) and nR == int(r.nR), \
            f"Assertion A FAILED at {r.pair} r{region_id}: ({nL},{nR}) vs ({r.nL},{r.nR})"
        N_L = int(Gm[mS, gidx[Lg]].sum())
        N_R = int(Gm[mR, gidx[Rg]].sum())
        T_obs = t_log2(N_L, N_R, nL, nR)
        # NOTE: r["T"] not r.T — .T on a pandas row is the transpose property
        assert abs(T_obs - float(r["T"])) <= 1e-9, \
            f"Assertion B FAILED at {r.pair} r{region_id}: {T_obs} vs {r['T']}"
        obs.append((gidx[Lg], gidx[Rg], tS, tR, nL, nR, T_obs))
    log("Assertions A+B PASSED: labels/regions/T_obs reproduce phase-2 outputs "
        "exactly on all 75 rows")

    # ---------- blocks (scout convention) ----------
    dec = np.zeros(n_pc, dtype=np.int64)
    for rr in range(3):
        idx_r = np.flatnonzero(regions == rr)
        q = np.quantile(tot[idx_r], np.linspace(0, 1, 11))
        dec[idx_r] = np.clip(np.digitize(tot[idx_r], q[1:-1]), 0, 9)
    blocks = []
    for rr in range(3):
        for d in range(10):
            idx = np.flatnonzero((regions == rr) & (dec == d))
            if len(idx) > 1:
                blocks.append(idx)
    log(f"blocks: {len(blocks)}; sizes {min(map(len, blocks))}-{max(map(len, blocks))}")

    # ---------- pre-flight determinism (Assertion D) ----------
    ctx = mp.get_context("spawn")
    init_args = (Gm, cell_slot, blocks, obs, n_slots)
    pool = ctx.Pool(N_WORKERS, initializer=worker_init, initargs=init_args)
    worker_init(*init_args)
    ref = permute_and_score(0)
    w1 = pool.apply(permute_and_score, (0,))
    assert np.array_equal(ref[1], w1[1]), "pre-flight: parent vs worker mismatch"
    w2 = pool.apply(permute_and_score, (0,))
    assert np.array_equal(w1[1], w2[1]), "pre-flight: worker-worker mismatch"
    log("Assertion D PASSED: perm 0 bitwise identical across parent + 2 workers")

    # ---------- parallel sweep with checkpointing ----------
    T_all = np.full((B_PERM, len(obs)), np.nan)
    T_all[0] = ref[1]
    done = {0}
    ckpt_path = f"{OUT}/phase3_null_checkpoint.npz"
    if not SMOKE and os.path.exists(ckpt_path):
        d0 = np.load(ckpt_path)
        prev = set(d0["done"].tolist())
        for p in prev:
            T_all[p] = d0["T"][p]
        assert np.array_equal(T_all[0], ref[1]), "resume: perm 0 mismatch"
        done |= prev
        log(f"resumed from checkpoint: {len(done)}/{B_PERM}")
    todo = [p for p in range(1, B_PERM) if p not in done]
    t_last = time.time()
    for p, T in pool.imap_unordered(permute_and_score, todo, chunksize=4):
        T_all[p] = T
        done.add(p)
        if len(done) % CKPT_EVERY == 0 and not SMOKE:
            tmp = ckpt_path + ".tmp.npz"
            with open(tmp, "wb") as fh:
                np.savez(fh, T=T_all, done=np.array(sorted(done), dtype=np.int64))
            os.replace(tmp, ckpt_path)
            rate = (time.time() - t_last) / CKPT_EVERY
            t_last = time.time()
            log(f"{len(done)}/{B_PERM} perms | {rate:.2f}s/perm eff | "
                f"ETA {rate * (B_PERM - len(done)) / 60:.1f} min")
    pool.close()
    pool.join()
    if not SMOKE:
        with open(ckpt_path, "wb") as fh:
            np.savez(fh, T=T_all, done=np.array(sorted(done), dtype=np.int64))
    assert not np.isnan(T_all).any(), "missing permutations"
    log(f"all {B_PERM} permutations complete")

    # ---------- p-values, BH, gate ----------
    os.makedirs(OUT, exist_ok=True)
    T_obs_arr = np.array([o[6] for o in obs])
    cnt_pos = (T_all >= T_obs_arr[None, :] - 1e-9).sum(axis=0)
    cnt_neg = (T_all <= T_obs_arr[None, :] + 1e-9).sum(axis=0)
    rows = []
    for k, o in enumerate(obs):
        p_pos_add = (cnt_pos[k] + 1) / (B_PERM + 1)
        p_neg_add = (cnt_neg[k] + 1) / (B_PERM + 1)
        rows.append(dict(pair=cov.pair.iloc[k], region=int(cov.region.iloc[k]),
                         T_obs=T_obs_arr[k],
                         p_pos_raw=cnt_pos[k] / B_PERM, p_neg_raw=cnt_neg[k] / B_PERM,
                         p_pos=p_pos_add, p_neg=p_neg_add,
                         T_perm_med=float(np.median(T_all[:, k]))))
    nul = pd.DataFrame(rows)
    nul["q_pos"] = bh_within_region(nul, "p_pos")
    nul["q_neg"] = bh_within_region(nul, "p_neg")

    cert_pos = cov.certified_pos_matched.to_numpy(bool)
    cert_neg = cov.certified_neg_matched.to_numpy(bool)
    sig_pos = cert_pos & (nul.q_pos.to_numpy() <= 0.10)
    sig_neg = cert_neg & (nul.q_neg.to_numpy() <= 0.10)
    n_cert = int(cert_pos.sum() + cert_neg.sum())
    n_sig = int(sig_pos.sum() + sig_neg.sum())
    frac = n_sig / n_cert
    gate = "PASS" if frac >= 0.85 else ("PARTIAL" if frac >= 0.60 else "FAIL")

    nul.to_csv(f"{OUT}/phase3_proseg_null.csv", index=False)
    summary = dict(
        prereg="B03_PHASE3_PREREG.md v1.0 (frozen 2026-09-07 09:40)",
        smoke=SMOKE, B=B_PERM, n_workers=N_WORKERS,
        seed_base=SEED_BASE,
        assertions=dict(
            A_label_region_reproduction="PASSED (75/75 nL/nR exact)",
            B_obs_crosscheck="PASSED (75/75 within 1e-9)",
            C_invariance=f"PASSED on perms {sorted(CHECK_PERMS)}",
            D_preflight_determinism="PASSED (bitwise, parent+2 workers)"),
        certified_rows=n_cert,
        significant_in_certified_direction=n_sig,
        fraction=frac, gate=gate,
        p_convention="add-one (count+1)/(B+1); BH within region on add-one p",
    )
    with open(f"{OUT}/phase3_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log(json.dumps({k: summary[k] for k in
                    ("certified_rows", "significant_in_certified_direction",
                     "fraction", "gate")}))
    print("PHASE3 NULL COMPLETE", flush=True)


if __name__ == "__main__":
    main()
