"""B03 Step 4 v3: matched robust permutation null — HARDWARE-AWARE PARALLEL RERUN.

Scientific protocol is IDENTICAL to v2 (see DEV-003); only the sampling loop changed:

  R1 PARALLELISM: permutations distributed over 4 worker processes (8 physical cores,
     5.7 GB free RAM measured before launch). Same B=1000, same within-region
     decile-matched scheme, same tie convention (Tlo >= T_lo_o - 1e-9), same interval
     machinery (cross-checked vs scout_bounds at 8.9e-16 in v2 and re-checked here).
  R2 REPRODUCIBLE PERMUTATION STREAM: permutation p is drawn from
     np.random.default_rng(20260905 + 100003 + p), blocks in the same fixed
     registration order as v2. (v2 used one sequential stream; a mid-run parallel
     continuation of that stream is not reproducible, so the rerun uses per-perm
     documented seeds. Documented as DEV-005. Estimator validity is unaffected:
     every permutation is a valid uniform within-block draw under H0.)
  R3 CHECKPOINTING: results saved every 50 completed permutations (atomic tmp+rename),
     so interruption never loses more than one chunk. v2's partial block_sums.npy
     (perms 0..239 of the old stream) is quarantined, not reused.
  R4 PRE-FLIGHT DETERMINISM CHECK: perm 0 computed twice (parent + a worker); bitwise
     equality required before the full run starts.

Outputs (identical schema to v2 so the watchdog consumes them unchanged):
  results/null_robust_matched.csv      (p_pos/p_neg as raw count/B — watchdog applies
                                        the add-one convention)
  results/scout_summary_matched.json
"""
import numpy as np
import pandas as pd
import pickle, json, time, math, os, sys
import multiprocessing as mp

# Tissue-2 reuse (DEV-017): seed base, dirs, and worker count overridable by env so
# spawned workers (which re-import this module fresh) inherit the registered
# tissue-2 configuration. Defaults unchanged -> tissue-1 reruns bitwise identical.
SEED = int(os.environ.get("B03_NULL_SEED", "20260905"))
SEED_BASE = SEED + 100003

EPS = 0.5
DATA = os.environ.get("B03_NULL_DATA", "B03_project/data")
OUT = os.environ.get("B03_NULL_OUT", "B03_project/results")
B_PERM = int(os.environ.get("B03_NULL_B", "1000"))
N_WORKERS = int(os.environ.get("B03_NULL_WORKERS", "4"))
CKPT_EVERY = int(os.environ.get("B03_NULL_CKPT", "50"))

POSITIVE_CONTROLS = [("CD274", "PDCD1"), ("CXCL12", "CXCR4"), ("CCL5", "CCR7"),
                     ("PECAM1", "KDR"), ("CD163", "CD3D")]
DISPUTED = [("ERBB2", "EGFR"), ("ERBB2", "PDCD1"), ("ESR1", "PGR"), ("PGR", "ESR1"),
            ("CD274", "CTLA4"), ("CD274", "CD8A"), ("CXCL12", "CCR7"), ("CCL5", "CXCR4"),
            ("PECAM1", "PDCD1"), ("PECAM1", "CTLA4"), ("CD163", "CCR7"), ("CD163", "CXCR4"),
            ("ACTA2", "EGFR"), ("ACTA2", "ERBB2"), ("KRT8", "ERBB2"), ("CDH1", "EGFR"),
            ("CDH1", "ERBB2"), ("CD3D", "ERBB2"), ("MS4A1", "CD274"), ("CD68", "CD274")]
ALL_PAIRS = POSITIVE_CONTROLS + DISPUTED
TYPE_MAP = {
    "CD274": "Breast cancer", "CXCL12": "Fibroblasts", "CCL5": "Macrophages",
    "PECAM1": "Endothelial cells", "CD163": "Macrophages", "ERBB2": "Breast cancer",
    "ESR1": "Breast cancer", "PGR": "Breast cancer", "ACTA2": "Smooth muscle cells",
    "KRT8": "Breast glandular cells", "CDH1": "Breast glandular cells",
    "CD3D": "T cells", "MS4A1": "B cells", "CD68": "Macrophages",
}
RECV_MAP = {
    "PDCD1": "T cells", "CXCR4": "T cells", "CCR7": "T cells",
    "KDR": "Endothelial cells", "CD3D": "T cells", "EGFR": "Breast cancer",
    "PGR": "Breast cancer", "ESR1": "Breast cancer", "CTLA4": "T cells",
    "CD8A": "T cells", "ERBB2": "Breast cancer", "CD274": "Breast cancer",
}

# ---------------------------------------------------------------------------
# shared state (set once in parent, transferred to workers via initializer)
# ---------------------------------------------------------------------------
_W = {}
# DEV-016 cap support (tissue 2): when B03_NULL_CAP=1, gain capacities are enforced
# with the SAME formula as the certify step: dmax = min(dmax_unc, sum_{c in target}
# min(pot_c(g), cap_c)). pot_by_gene/cap_per_cell arrive via worker_init. Default off
# -> tissue-1 behavior identical.
CAPS = os.environ.get("B03_NULL_CAP", "0") == "1"


def worker_init(Gmat, mov_cell, mov_gid, donor_flat, gene_slices, blocks, obs_arr,
                cid, caps=None):
    _W["Gmat"] = Gmat
    _W["mov_cell"] = mov_cell
    _W["donor_flat"] = donor_flat
    _W["gene_slices"] = gene_slices
    _W["blocks"] = blocks
    _W["obs"] = obs_arr  # list of (giL, giR, tS, tR, nLc, nRc, T_lo_o, T_hi_o)
    _W["cid"] = cid      # region-encoded frozen labels (permutation source array)
    _W["caps"] = caps    # None or (pot_by_gene ndarray (n_genes, n_cells), cap_per_cell)


def interval_for(gi, target, lab):
    Gmat = _W["Gmat"]
    mov_cell = _W["mov_cell"]
    donor_flat = _W["donor_flat"]
    mv, ps, pe = _W["gene_slices"][gi]
    cells_m = mov_cell[mv]
    src_in = np.where(cells_m >= 0, lab[np.maximum(cells_m, 0)], -999) == target
    d_in = lab[donor_flat] == target
    offs = ps - ps[0]
    seg = d_in[ps[0]:pe[-1]]
    if len(offs) and len(seg) and len(seg) >= len(offs):
        anyin = np.asarray(np.logical_or.reduceat(seg, offs))[:len(mv)]
        if len(anyin) < len(mv):
            anyin = np.pad(anyin, (0, len(mv) - len(anyin)))
    else:
        anyin = np.zeros(len(mv), bool)
    dmax = int(((~src_in) & anyin).sum())
    if _W.get("caps") is not None:
        pot, capc = _W["caps"]
        in_t = lab == target
        dmax_cap = int(np.minimum(pot[gi][in_t], capc[in_t]).sum())
        dmax = min(dmax, dmax_cap)
    dmin = int(src_in.sum())
    N = int(Gmat[np.flatnonzero(lab == target), gi].sum())
    return N, dmin, dmax


def permute_and_score(p):
    """One permutation: draw from seed SEED_BASE+p, score all 75 (pair,region) keys."""
    obs_arr = _W["obs"]
    cid = _W["cid"]
    lab_perm = cid.copy()
    rng = np.random.default_rng(SEED_BASE + p)
    for idx in _W["blocks"]:
        pm = rng.permutation(len(idx))
        lab_perm[idx] = lab_perm[idx[pm]]
    n_keys = len(obs_arr)
    Tlo = np.empty(n_keys)
    Thi = np.empty(n_keys)
    for k in range(n_keys):
        giL, giR, tS, tR, nLc, nRc, T_lo_o, T_hi_o = obs_arr[k]
        NL, dLmin, dLmax = interval_for(giL, tS, lab_perm)
        NR, dRmin, dRmax = interval_for(giR, tR, lab_perm)
        Tlo[k] = (math.log2(max(0, NL - dLmin) / nLc + EPS) +
                  math.log2(max(0, NR - dRmin) / nRc + EPS))
        Thi[k] = (math.log2((NL + dLmax) / nLc + EPS) +
                  math.log2((NR + dRmax) / nRc + EPS))
    return p, Tlo, Thi


def bh_within_region(df, pcol):
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


def main():
    t0 = time.time()
    tx = pd.read_parquet(f"{DATA}/tx.parquet")
    cells_meta = pd.read_parquet(f"{DATA}/cells_meta.parquet")
    with open(f"{DATA}/donor_map.pkl", "rb") as f:
        dm = pickle.load(f)
    labels = np.load(f"{DATA}/labels.npy", allow_pickle=True)
    region_of_cell = np.load(f"{DATA}/regions.npy")
    res = pd.read_csv(f"{OUT}/scout_bounds.csv")
    tested = sorted({g for p in ALL_PAIRS for g in p})
    gidx = {g: i for i, g in enumerate(tested)}

    # ---- alignment assertions (identical to v2) ----
    n_cells = len(cells_meta)
    assert n_cells == len(labels) == len(region_of_cell), "length misalignment"
    # Frozen label-count verification: against label_counts.json when present
    # (tissue-2 reuse, written by the tissue-2 loader), else the tissue-1 scout
    # counts (behavior identical to the original hard assertion).
    lc_path = os.path.join(DATA, "label_counts.json")
    if os.path.exists(lc_path):
        expected = json.load(open(lc_path))
        vc_d = pd.Series(labels).value_counts().to_dict()
        for k, v in expected.items():
            assert vc_d.get(k) == v, f"label count mismatch for {k}: {vc_d.get(k)} != {v}"
        print(f"label counts verified against {lc_path}", flush=True)
    else:
        vc = pd.Series(labels).value_counts()
        assert vc.get("Breast cancer") == 60459 and vc.get("T cells") == 20006, \
            f"label counts differ from scout run: {vc.head(3).to_dict()}"
    print(f"[{time.time()-t0:.0f}s] alignment verified", flush=True)

    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()
    Gmat = np.zeros((n_cells, len(tested)), dtype=np.int64)
    for g in tested:
        sel_g = (feat == g) & (cell >= 0)
        Gmat[:, gidx[g]] = np.bincount(cell[sel_g], minlength=n_cells)

    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    n_types = len(types_needed)
    SHIFT = 1
    type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
    n_slots = n_types + SHIFT
    cid = region_of_cell * n_slots + np.array(
        [type_shift.get(l, 0) for l in labels], dtype=np.int64)

    band_tested = dm["band_tested_indices"]
    donor_sets = dm["donor_sets"]
    mov_gene = feat[band_tested]
    mov_cell = cell[band_tested]
    lens = np.array([len(d) for d in donor_sets])
    donor_flat = np.concatenate([np.asarray(d, dtype=np.int64) for d in donor_sets
                                 if len(d)])
    donor_ptr = np.r_[0, np.cumsum(lens)]
    mov_gid = np.array([gidx[g] for g in mov_gene], dtype=np.int64)
    mov_by_gene = {gi: np.flatnonzero(mov_gid == gi) for gi in range(len(tested))}
    print(f"[{time.time()-t0:.0f}s] movable {len(band_tested):,}", flush=True)

    # ---- registered blocks: (region, decile of total transcripts) ----
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
    print(f"blocks: {len(blocks)}; sizes {min(map(len, blocks))}-"
          f"{max(map(len, blocks))}", flush=True)

    # ---- observed intervals, single code path (identical to v2) ----
    gene_slices = {gi: (mv, donor_ptr[mv], donor_ptr[mv + 1])
                   for gi, mv in mov_by_gene.items()}

    def interval_obs(gi, target, lab):
        mv, ps, pe = gene_slices[gi]
        cells_m = mov_cell[mv]
        src_in = np.where(cells_m >= 0, lab[np.maximum(cells_m, 0)], -999) == target
        d_in = lab[donor_flat] == target
        offs = ps - ps[0]
        seg = d_in[ps[0]:pe[-1]]
        if len(offs) and len(seg) and len(seg) >= len(offs):
            anyin = np.asarray(np.logical_or.reduceat(seg, offs))[:len(mv)]
            if len(anyin) < len(mv):
                anyin = np.pad(anyin, (0, len(mv) - len(anyin)))
        else:
            anyin = np.zeros(len(mv), bool)
        dmax = int(((~src_in) & anyin).sum())
        if CAPS and pot_by_gene is not None:
            in_t = lab == target
            dmax_cap = int(np.minimum(pot_by_gene[gi][in_t], cap_per_cell[in_t]).sum())
            dmax = min(dmax, dmax_cap)
        dmin = int(src_in.sum())
        N = int(Gmat[np.flatnonzero(lab == target), gi].sum())
        return N, dmin, dmax

    keys = []
    obs_arr = []
    pot_by_gene = None
    cap_per_cell = None
    if CAPS:
        # DEV-016: capacity inputs (same formula as the certify step)
        cap_per_cell = np.floor(0.5 * cells_meta.transcript_counts.to_numpy(np.int64)).astype(np.int64)
        pot_by_gene = np.zeros((len(tested), n_cells), dtype=np.int64)
        band_tested_p = dm["band_tested_indices"]
        donor_sets_p = dm["donor_sets"]
        mov_gid_p = np.array([gidx[g] for g in feat[band_tested_p]], dtype=np.int64)
        lens_p = np.array([len(d) for d in donor_sets_p])
        donor_flat_p = np.concatenate([np.asarray(d, dtype=np.int64) for d in donor_sets_p if len(d)])
        dptr_p = np.r_[0, np.cumsum(lens_p)]
        for gi in range(len(tested)):
            m_g = np.flatnonzero(mov_gid_p == gi)
            segs = [donor_flat_p[dptr_p[i]:dptr_p[i + 1]] for i in m_g
                    if dptr_p[i + 1] > dptr_p[i]]
            slots = np.concatenate(segs) if segs else np.array([], np.int64)
            pot_by_gene[gi] = np.bincount(slots, minlength=n_cells).astype(np.int64)
        print(f"[{time.time()-t0:.0f}s] DEV-016 caps: computed pot_by_gene "
              f"{pot_by_gene.shape}", flush=True)
    for _, r in res.iterrows():
        Lg, Rg = r.pair.split("->")
        region_id = int(r.region)
        tS = region_id * n_slots + type_shift[TYPE_MAP[Lg]]
        tR = region_id * n_slots + type_shift[RECV_MAP[Rg]]
        NL, dLmin, dLmax = interval_obs(gidx[Lg], tS, cid)
        NR, dRmin, dRmax = interval_obs(gidx[Rg], tR, cid)
        nLc, nRc = int(r.nl_cells), int(r.nr_cells)
        T_lo = (math.log2(max(0, NL - dLmin) / nLc + EPS) +
                math.log2(max(0, NR - dRmin) / nRc + EPS))
        T_hi = (math.log2((NL + dLmax) / nLc + EPS) +
                math.log2((NR + dRmax) / nRc + EPS))
        worst = max(abs(T_lo - float(r.T_lo)), abs(T_hi - float(r.T_hi)))
        assert worst < 1e-6, f"observed interval mismatch for {r.pair} r{region_id}: {worst}"
        keys.append((r.pair, region_id))
        obs_arr.append((gidx[Lg], gidx[Rg], tS, tR, nLc, nRc,
                        float(r.T_lo), float(r.T_hi)))
    T_lo_obs = np.array([o[6] for o in obs_arr])
    T_hi_obs = np.array([o[7] for o in obs_arr])
    print(f"[{time.time()-t0:.0f}s] observed intervals re-verified vs scout_bounds "
          f"(all 75 within 1e-6)", flush=True)

    # ---- pre-flight determinism check (R4) ----
    ctx = mp.get_context("spawn")
    init_args = (Gmat, mov_cell, mov_gid, donor_flat, gene_slices, blocks, obs_arr,
                 cid, (pot_by_gene, cap_per_cell) if CAPS else None)
    pool = ctx.Pool(N_WORKERS, initializer=worker_init, initargs=init_args)
    # parent-side state for its own perm-0 computation (pre-flight reference)
    worker_init(*init_args)
    ref = permute_and_score(0)  # parent computes perm 0 itself
    w1 = pool.apply(permute_and_score, (0,))
    w2 = pool.apply(permute_and_score, (0,))
    assert np.array_equal(ref[1], w1[1]) and np.array_equal(ref[2], w1[2]), \
        "pre-flight: parent vs worker mismatch"
    assert np.array_equal(w1[1], w2[1]) and np.array_equal(w1[2], w2[2]), \
        "pre-flight: worker-worker mismatch"
    print(f"[{time.time()-t0:.0f}s] pre-flight determinism check PASSED "
          f"(perm 0 bitwise identical across parent + 2 workers)", flush=True)

    # ---- parallel permutation sweep with checkpointing (R1, R3) ----
    Tlo_all = np.full((B_PERM, len(keys)), np.nan)
    Thi_all = np.full((B_PERM, len(keys)), np.nan)
    Tlo_all[0] = ref[1]
    Thi_all[0] = ref[2]
    done = {0}
    ckpt_path = f"{OUT}/null_v3_checkpoint.npz"
    # resume: np.savez auto-appends '.npz' to a requested filename, so the crash
    # left the first checkpoint at '*.tmp.npz'; its rows carry per-perm seeds and
    # are bitwise reproducible, hence valid to reuse. (DEV-005 note added.)
    legacy = ckpt_path + ".tmp.npz"
    if os.path.exists(ckpt_path):
        d0 = np.load(ckpt_path)
        prev_done = set(d0["done"].tolist())
        for p in prev_done:
            Tlo_all[p] = d0["Tlo"][p]
            Thi_all[p] = d0["Thi"][p]
        assert np.array_equal(Tlo_all[0], ref[1]) and np.array_equal(Thi_all[0], ref[2]), \
            "resume: checkpointed perm 0 disagrees with fresh parent computation"
        done |= prev_done
        print(f"[{time.time()-t0:.0f}s] resumed from checkpoint: "
              f"{len(done)}/{B_PERM} permutations already valid", flush=True)
    elif os.path.exists(legacy):
        d0 = np.load(legacy)
        fin = np.isfinite(d0["Tlo"]).all(axis=1) & np.isfinite(d0["Thi"]).all(axis=1)
        prev_done = set(np.flatnonzero(fin).tolist())
        for p in prev_done:
            Tlo_all[p] = d0["Tlo"][p]
            Thi_all[p] = d0["Thi"][p]
        assert np.array_equal(Tlo_all[0], ref[1]) and np.array_equal(Thi_all[0], ref[2]), \
            "resume: checkpointed perm 0 disagrees with fresh parent computation"
        done |= prev_done
        print(f"[{time.time()-t0:.0f}s] resumed from crash legacy checkpoint: "
              f"{len(done)}/{B_PERM} permutations already valid", flush=True)
    todo = [p for p in range(1, B_PERM) if p not in done]
    t_last = time.time()
    for p, Tlo, Thi in pool.imap_unordered(permute_and_score, todo, chunksize=4):
        Tlo_all[p] = Tlo
        Thi_all[p] = Thi
        done.add(p)
        if len(done) % CKPT_EVERY == 0:
            tmp = ckpt_path + ".tmp.npz"  # savez appends '.npz' to requested names
            with open(tmp, "wb") as fh:
                np.savez(fh, Tlo=Tlo_all, Thi=Thi_all,
                         done=np.array(sorted(done), dtype=np.int64))
            os.replace(tmp, ckpt_path)
            rate = (time.time() - t_last) / CKPT_EVERY
            t_last = time.time()
            eta_min = rate * (B_PERM - len(done)) / 60
            print(f"  {len(done)}/{B_PERM} perms | {rate:.1f}s/perm eff | "
                  f"ETA {eta_min:.0f} min", flush=True)
    pool.close()
    pool.join()
    with open(ckpt_path, "wb") as fh:
        np.savez(fh, Tlo=Tlo_all, Thi=Thi_all,
                 done=np.array(sorted(done), dtype=np.int64))
    assert not np.isnan(Tlo_all).any(), "missing permutations"
    print(f"[{time.time()-t0:.0f}s] all {B_PERM} permutations complete", flush=True)

    # ---- p-values: identical convention to v2 (watchdog applies add-one) ----
    p_pos_acc = (Tlo_all >= T_lo_obs[None, :] - 1e-9).sum(axis=0)
    p_neg_acc = (Thi_all <= T_hi_obs[None, :] + 1e-9).sum(axis=0)
    rows = []
    for k, (pair, region_id) in enumerate(keys):
        rows.append(dict(pair=pair, region=region_id,
                         p_pos=p_pos_acc[k] / B_PERM,
                         p_neg=p_neg_acc[k] / B_PERM,
                         T_lo_perm_med=float(np.median(Tlo_all[:, k])),
                         T_hi_perm_med=float(np.median(Thi_all[:, k]))))
    nul = pd.DataFrame(rows)
    nul["q_pos"] = bh_within_region(nul, "p_pos")
    nul["q_neg"] = bh_within_region(nul, "p_neg")
    nul.to_csv(f"{OUT}/null_robust_matched.csv", index=False)
    print(f"[{time.time()-t0:.0f}s] null_robust_matched.csv written", flush=True)

    summary = {
        "method": "matched robust permutation null (DEV-003), within-region "
                  "decile-matched labels; PARALLEL RERUN per DEV-005 "
                  "(per-perm seeds 20260905+100003+p, pre-flight determinism check)",
        "B": B_PERM,
        "n_workers": N_WORKERS,
        "alignment": "labels/regions polygon-index space; verified by assertion",
        "obs_crosscheck": "all 75 rows within 1e-6 of scout_bounds.csv",
        "preflight": "perm 0 bitwise identical across parent + 2 workers",
        "n_pairs": len(ALL_PAIRS), "n_regions": 3,
    }
    with open(f"{OUT}/scout_summary_matched.json", "w") as fj:
        json.dump(summary, fj, indent=2)
    print("DONE — watchdog will now produce the scout report", flush=True)


if __name__ == "__main__":
    main()
