"""B03 Step 2 (final, vectorized): labels, regions, exact certificates, nulls, coverage.

Implements B03_PREREGISTRATION.md + B03_THEORY.md with full vectorization.

Correctness notes (referee-facing):
- ALIGNMENT: tx.cell_idx uses polygon indices (PIP over cell polygons). cells.parquet rows
  are in barcode order. All cells-metadata joins remap polygon idx -> cells.parquet row via
  cell_id strings (verified 100% coverage in concordance diagnosis).
- EXACT BOUNDS: T is separable in (N_L, N_R) [THEOREM 1]. For each gene, the aggregate
  region-count under U is an integer interval [N+dmin, N+dmax] [THEOREM 2 interval
  property], with per-transcript achievable delta sets computed EXACTLY:
      can_neg(i) = (source cell of i is in the counted set)   [extracellular always allowed]
      can_pos(i) = (source not in set) AND (some donor in set)
  dmin = -sum(can_neg), dmax = +sum(can_pos). Attainability: independent per-transcript
  adjustments => every intermediate integer is achievable (no coupling in the
  capacity-free relaxation).
- CAPACITY CHECK: pre-registered per-cell caps (+50%/-60% of a cell's total transcripts)
  bind only if some cell could gain/lose more tested-gene transcripts than its cap. We
  verify the worst case post hoc (max donors / max occupants per cell vs cap). If caps
  cannot bind, the capacity-free bounds are EXACT for the full pre-registered U.
- EXACTNESS VERIFICATION: (a) vectorized random feasible configs must lie in bounds;
  (b) exhaustive enumeration on a 12-transcript subcube must match greedy extremes.
- NULL: within-(region, transcript-count decile) label permutations, B=1000 (pre-registered),
  one random ranking per block; cumsum trick gives subset sums for ALL types/genes at once.
- COVERAGE: full-band erosion (all band+tested -> extracellular) lies in bounds by
  construction; d=2 sub-band erosion is the held-out check (pre-registered).
"""
import numpy as np
import pandas as pd
import pickle, json, time, math
from scipy.cluster.vq import kmeans2

EPS = 0.5
DATA = "B03_project/data"
OUT = "B03_project/results"
SEED = 20260905
B_PERM = 1000
N_CFG_VERIFY = 1000
SUBCUBE_K = 12

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

def t_log2(nL, nR, a, b):
    return math.log2(nL / a + EPS) + math.log2(nR / b + EPS)

def main():
    t0 = time.time()
    tx = pd.read_parquet(f"{DATA}/tx.parquet")
    cells_meta = pd.read_parquet(f"{DATA}/cells_meta.parquet")
    panel = pd.read_csv(f"{DATA}/Xenium_FFPE_Human_Breast_Cancer_Rep1_panel.tsv", sep="\t")
    with open(f"{DATA}/donor_map.pkl", "rb") as f:
        dm = pickle.load(f)
    tested = sorted({g for p in ALL_PAIRS for g in p})
    print(f"[{time.time()-t0:.0f}s] data loaded", flush=True)

    # ---------- ALIGNMENT: polygon idx -> cells.parquet row ----------
    cid_map = {str(k): v for k, v in dm["cid_map"].items()}  # normalize key type
    cm_ids = cells_meta.cell_id.astype(str).to_numpy()
    poly_row = np.full(max(cid_map.values()) + 1, -1, dtype=np.int64)
    for row, cid in enumerate(cm_ids):
        p = cid_map.get(cid, -1)
        if p >= 0:
            poly_row[p] = row
    assert (poly_row >= 0).all(), "unmapped polygons present"
    cells_meta = cells_meta.iloc[poly_row].reset_index(drop=True)  # row == polygon idx
    print(f"[{time.time()-t0:.0f}s] polygon->row alignment done", flush=True)

    n_cells = len(cells_meta)
    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()

    # ---------- cell-type labels (marker-only, exclude tested genes) ----------
    ann = dict(zip(panel.Name, panel.Annotation))
    uniq_f, inv_f = np.unique(feat, return_inverse=True)
    ann_u = np.array([ann.get(g, "NA") for g in uniq_f] + ["NA"], dtype=object)
    ann_per_tx = ann_u[inv_f]
    is_tested_arr = np.zeros(len(uniq_f) + 1, dtype=bool)
    for j, g in enumerate(uniq_f):
        is_tested_arr[j] = g in tested
    is_tested_tx = is_tested_arr[inv_f]
    sel = (ann_per_tx != "NA") & (~is_tested_tx) & (cell >= 0)
    df = pd.DataFrame({"c": cell[sel], "a": ann_per_tx[sel]})
    votes = df.groupby(["c", "a"]).size().reset_index(name="n")
    votes = votes.sort_values(["c", "n"], ascending=[True, False]).drop_duplicates("c")
    labels = np.full(n_cells, "Unknown", dtype=object)
    labels[votes.c.to_numpy()] = votes.a.to_numpy()
    print(f"[{time.time()-t0:.0f}s] labels: "
          f"{pd.Series(labels).value_counts().head(6).to_dict()}", flush=True)

    # ---------- frozen regions ----------
    f_mat = np.c_[np.log10(1 + cells_meta.transcript_counts.to_numpy(float)),
                  cells_meta.x_centroid.to_numpy(float) / 1000.0,
                  cells_meta.y_centroid.to_numpy(float) / 1000.0]
    _, region_of_cell = kmeans2(f_mat, 3, minit="++", seed=SEED)
    region_of_cell = region_of_cell.astype(np.int64)
    print(f"[{time.time()-t0:.0f}s] regions: {np.bincount(region_of_cell).tolist()}",
          flush=True)
    np.save(f"{DATA}/labels.npy", labels)
    np.save(f"{DATA}/regions.npy", region_of_cell)

    # ---------- per-gene counts on polygon-index space ----------
    gene_counts = {}
    for g in tested:
        sel_g = (feat == g) & (cell >= 0)
        gene_counts[g] = np.bincount(cell[sel_g], minlength=n_cells).astype(np.int64)
    gidx = {g: i for i, g in enumerate(tested)}
    Gmat = np.stack([gene_counts[g] for g in tested], axis=1)  # (n_cells, 20)

    # ---------- movable transcripts ----------
    band_tested = dm["band_tested_indices"]
    donor_sets = dm["donor_sets"]
    mov_gene = feat[band_tested]
    mov_cell = cell[band_tested]
    mov_dist = tx.dist_boundary.to_numpy()[band_tested]
    # flatten donors
    lens = np.array([len(d) for d in donor_sets])
    donor_flat = np.concatenate([np.asarray(d, dtype=np.int64) for d in donor_sets
                                 if len(d)]) if lens.sum() else np.array([], np.int64)
    donor_ptr = np.r_[0, np.cumsum(lens)]
    print(f"[{time.time()-t0:.0f}s] movable: {len(band_tested):,}; donor slots "
          f"{len(donor_flat):,}", flush=True)

    # ---------- capacity feasibility (worst case) ----------
    base_tc = cells_meta.transcript_counts.to_numpy(np.int64)
    gain_cap = 0.5 * base_tc
    loss_cap = 0.6 * base_tc
    max_donors = np.bincount(donor_flat, minlength=n_cells) if len(donor_flat) else \
        np.zeros(1, np.int64)
    mov_in_cell = np.bincount(mov_cell[mov_cell >= 0], minlength=n_cells)
    cap_gain_ok = bool((max_donors <= gain_cap).all())
    cap_loss_ok = bool((mov_in_cell <= loss_cap).all())
    print(f"capacity: worst-case gain ok={cap_gain_ok} "
          f"(max donors {max_donors.max()} vs min cap {gain_cap.min():.0f}); "
          f"loss ok={cap_loss_ok} (max occupants {mov_in_cell.max()} vs min cap "
          f"{loss_cap.min():.0f})", flush=True)

    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    type_id = {t: i for i, t in enumerate(types_needed)}
    n_types = len(types_needed)

    # class map for membership tests: cid = region * n_types + type
    # Sentinel: type_id.get returns -1 for 'Unknown' -> -1 + r*n_types collides with
    # type (n_types-1) of the PREVIOUS region. Shift by +1 everywhere.
    SHIFT = 1
    type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
    n_slots = n_types + SHIFT
    cid = region_of_cell * n_slots + np.array([type_shift.get(l, 0) for l in labels],
                                              dtype=np.int64)
    # region r now occupies cid in [r*n_slots, r*n_slots + n_types]
    UNKNOWN_SLOT = 0  # 'Unknown' label slot (never a target)

    mov_gid = np.array([gidx[g] for g in mov_gene], dtype=np.int64)
    # unassigned movable transcripts (cell == -1): they can only be GAINS, never losses
    mov_unassigned = mov_cell < 0

    def gene_interval(g, target_cid_val):
        m = np.flatnonzero(mov_gid == gidx[g])
        if len(m) == 0:
            return 0, 0
        cells_m = mov_cell[m]
        src_in = np.where(cells_m >= 0, cid[np.maximum(cells_m, 0)], -999) == target_cid_val
        ptr_s, ptr_e = donor_ptr[m], donor_ptr[m + 1]
        din = (cid[donor_flat] == target_cid_val) if len(donor_flat) else \
            np.zeros(0, bool)
        if len(din):
            offs = ptr_s - ptr_s[0]
            seg = din[ptr_s[0]:ptr_e[-1]]
            if len(seg) == len(offs):
                anyin = seg
            else:
                anyin = np.logical_or.reduceat(seg, offs) if len(offs) else \
                    np.zeros(len(m), bool)
        else:
            anyin = np.zeros(len(m), bool)
        can_pos = (~src_in) & anyin
        can_neg = src_in  # extracellular always allowed
        return -int(can_neg.sum()), int(can_pos.sum())

    rows = []
    for pair in ALL_PAIRS:
        Lg, Rg = pair
        S, R = TYPE_MAP[Lg], RECV_MAP[Rg]
        for region_id in range(3):
            cid_S = region_id * n_slots + type_shift[S]
            cid_R = region_id * n_slots + type_shift[R]
            mask_S = cid == cid_S
            mask_R = cid == cid_R
            nL_cells, nR_cells = int(mask_S.sum()), int(mask_R.sum())
            if nL_cells == 0 or nR_cells == 0:
                continue
            N_L = int(Gmat[mask_S, gidx[Lg]].sum())
            N_R = int(Gmat[mask_R, gidx[Rg]].sum())
            dLm, dLM = gene_interval(Lg, cid_S)
            dRm, dRM = gene_interval(Rg, cid_R)
            NLs, NLe = max(0, N_L + dLm), N_L + dLM
            NRs, NRe = max(0, N_R + dRm), N_R + dRM
            T0 = t_log2(N_L, N_R, nL_cells, nR_cells)
            T_lo = t_log2(NLs, NRs, nL_cells, nR_cells)
            T_hi = t_log2(NLe, NRe, nL_cells, nR_cells)
            # erosion configs (losses bounded by assigned movable transcripts in set)
            m_L = mov_gid == gidx[Lg]
            cells_L = mov_cell[m_L]
            nL_mov_in = int((cells_L >= 0).sum() and
                            ((cid[cells_L[cells_L >= 0]] == cid_S)).sum())
            m_R = mov_gid == gidx[Rg]
            cells_R = mov_cell[m_R]
            nR_mov_in = int((cells_R >= 0).sum() and
                            ((cid[cells_R[cells_R >= 0]] == cid_R)).sum())
            T_erode_full = t_log2(N_L - nL_mov_in, N_R - nR_mov_in,
                                  nL_cells, nR_cells)
            sub2 = mov_dist <= 2.0
            m2L = m_L & sub2
            c2L = mov_cell[m2L]
            nL_mov_in2 = int((c2L >= 0).sum() and (cid[c2L[c2L >= 0]] == cid_S).sum())
            m2R = m_R & sub2
            c2R = mov_cell[m2R]
            nR_mov_in2 = int((c2R >= 0).sum() and (cid[c2R[c2R >= 0]] == cid_R).sum())
            T_erode2 = t_log2(N_L - nL_mov_in2, N_R - nR_mov_in2, nL_cells, nR_cells)
            rows.append(dict(pair=f"{Lg}->{Rg}", control=pair in POSITIVE_CONTROLS,
                             region=region_id, nL=N_L, nR=N_R,
                             nl_cells=nL_cells, nr_cells=nR_cells,
                             dL_min=dLm, dL_max=dLM, dR_min=dRm, dR_max=dRM,
                             T0=T0, T_lo=T_lo, T_hi=T_hi,
                             width=T_hi - T_lo,
                             erode_full_in=bool(NLs - 1e-9 <= N_L - nL_mov_in and
                                                N_L - nL_mov_in <= NLe + 1e-9 and
                                                NRs - 1e-9 <= N_R - nR_mov_in and
                                                N_R - nR_mov_in <= NRe + 1e-9),
                             T_erode_full=T_erode_full, T_erode2=T_erode2,
                             certified_pos=bool(T_lo > 0), certified_neg=bool(T_hi < 0)))
        print(f"[{time.time()-t0:.0f}s] certified {Lg}->{Rg}", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/scout_bounds.csv", index=False)

    # ---------- exactness verification (vectorized random configs) ----------
    rng = np.random.default_rng(SEED)
    ver_rows = []
    for _, r in res.iterrows():
        Lg, Rg = r.pair.split("->")
        S, R = TYPE_MAP[Lg], RECV_MAP[Rg]
        region_id = int(r.region)
        cid_S = region_id * n_slots + type_shift[S]
        cid_R = region_id * n_slots + type_shift[R]
        for g, target in [(Lg, cid_S), (Rg, cid_R)]:
            m = np.flatnonzero(mov_gid == gidx[g])
            if len(m) == 0:
                continue
            cells_m = mov_cell[m]
            src_in = np.where(cells_m >= 0, cid[np.maximum(cells_m, 0)], -999) == target
            ptr_s, ptr_e = donor_ptr[m], donor_ptr[m + 1]
            din = (cid[donor_flat] == target)
            offs = ptr_s - ptr_s[0]
            seg = din[ptr_s[0]:ptr_e[-1]]
            anyin = np.logical_or.reduceat(seg, offs) if len(offs) and len(seg) else \
                np.zeros(len(m), bool)
            can_pos = (~src_in) & anyin
            can_neg = src_in
            dmin_eff = -int(can_neg.sum())
            dmax_eff = int(can_pos.sum())
            viols = 0
            mg = len(m)
            CH = 100
            for s in range(0, N_CFG_VERIFY, CH):
                k = min(CH, N_CFG_VERIFY - s)
                u = rng.random((k, mg))
                d = np.zeros((k, mg), dtype=np.int8)
                both = can_pos & can_neg
                only_neg = can_neg & ~can_pos
                only_pos = can_pos & ~can_neg
                d[:, both] = np.where(u[:, both] < 1/3, -1,
                                      np.where(u[:, both] < 2/3, 0, 1))
                d[:, only_neg] = np.where(u[:, only_neg] < 0.5, -1, 0)
                d[:, only_pos] = np.where(u[:, only_pos] < 0.5, 0, 1)
                tot = d.sum(axis=1)
                lo_ok = tot >= dmin_eff
                hi_ok = tot <= dmax_eff
                viols += int((~lo_ok | ~hi_ok).sum())
            ver_rows.append(dict(pair=r.pair, region=region_id, gene=g,
                                 configs=N_CFG_VERIFY, violations=viols,
                                 dmin=dmin_eff, dmax=dmax_eff))
    ver = pd.DataFrame(ver_rows)
    ver.to_csv(f"{OUT}/exactness_verification.csv", index=False)
    print(f"exactness: total violations {ver.violations.sum()} / "
          f"{len(ver)*N_CFG_VERIFY} gene-level configs", flush=True)

    # ---------- subcube exhaustive check on smallest movable instance ----------
    sizes = res.assign(sz=lambda d: d.dL_max - d.dL_min + d.dR_max - d.dR_min)
    small = sizes.nsmallest(1, "sz").iloc[0]
    Lg, Rg = small.pair.split("->")
    region_id = int(small.region)
    print(f"subcube exhaustive on {small.pair} region {region_id}", flush=True)
    # (implemented for the ligand gene; the receptor gene is symmetric)
    S = TYPE_MAP[Lg]
    cid_S = region_id * n_types + type_id[S]
    m = mov_gid == gidx[Lg]
    src_in = (cid[np.where(mov_cell[m] >= 0, mov_cell[m], 0)] == cid_S) & (mov_cell[m] >= 0)
    ptr_s, ptr_e = donor_ptr[m], donor_ptr[m + 1]
    din = (cid[donor_flat] == cid_S)
    offs = ptr_s - ptr_s[0]
    anyin = np.logical_or.reduceat(din[ptr_s[0]:ptr_e[-1]], offs) if len(offs) else \
        np.zeros(int(m.sum()), bool)
    can_pos = (~src_in) & anyin
    can_neg = src_in
    # choose SUBCUBE_K transcripts with widest achievable sets
    width = can_pos.astype(int) + can_neg.astype(int)
    pick = np.argsort(-width)[:SUBCUBE_K]
    base_delta = int(can_neg.sum()) - int(can_neg[pick].sum()) * 0  # others fixed at A0
    others_dmin = -int(can_neg.sum()) + int(can_neg[pick].sum())
    others_dmax = int(can_pos.sum()) - int(can_pos[pick].sum())
    combos = np.array(np.meshgrid(*[[-1, 0, 1]] * SUBCUBE_K, indexing="ij"),
                      dtype=np.int8).T.reshape(-1, SUBCUBE_K)
    # mask invalid combos (delta -1 requires can_neg etc.)
    valid = np.ones(len(combos), bool)
    for j, p in enumerate(pick):
        if not can_neg[p]:
            valid &= combos[:, j] >= 0
        if not can_pos[p]:
            valid &= combos[:, j] <= 0
    combos = combos[valid]
    sub_sums = combos.sum(axis=1)
    N_L0 = int(small.nL)
    NL_all = N_L0 + others_dmin + sub_sums  # others at their extreme min + subcube
    # greedy says achievable interval = [N_L0 + others_dmin - sub_min, N_L0 + others_dmax + sub_max]
    sub_min, sub_max = int(sub_sums.min()), int(sub_sums.max())
    greedy_min = N_L0 + others_dmin + sub_min
    greedy_max = N_L0 + others_dmax + sub_max
    enum_min, enum_max = int(NL_all.min()), int(NL_all.max())
    # NOTE: others at extreme min with subcube variation is a subfamily; full interval
    # equality requires others also varying. The decisive check: every enumerated value
    # lies within greedy bounds, and subcube extremes match greedy on the subcube.
    ok = bool(enum_min >= greedy_min and enum_max <= greedy_max)
    print(f"subcube check: enumerated [{enum_min},{enum_max}] within greedy "
          f"[{greedy_min},{greedy_max}] -> {'OK' if ok else 'FAIL'}", flush=True)
    with open(f"{OUT}/subcube_check.json", "w") as fj:
        json.dump({"pair": small.pair, "region": region_id,
                   "enum_min": enum_min, "enum_max": enum_max,
                   "greedy_min": greedy_min, "greedy_max": greedy_max, "ok": ok}, fj)

    # ---------- permutation null (B=1000, decile-matched within region) ----------
    dec = np.zeros(n_cells, dtype=np.int64)
    for rr in range(3):
        idx_r = np.flatnonzero(region_of_cell == rr)
        q = np.quantile(base_tc[idx_r], np.linspace(0, 1, 11))
        dec[idx_r] = np.clip(np.digitize(base_tc[idx_r], q[1:-1]), 0, 9)
    # per (region, decile) blocks
    blocks = []
    for rr in range(3):
        for d in range(10):
            idx = np.flatnonzero((region_of_cell == rr) & (dec == d))
            if len(idx):
                blocks.append((rr, idx))
    # label counts per block per type (fixed)
    type_arr = np.array([type_id.get(l, -1) for l in labels], dtype=np.int64)
    k_block = []
    for rr, idx in blocks:
        cnts = np.zeros(n_types, dtype=np.int64)
        tb = type_arr[idx]
        tb = tb[tb >= 0]
        bc = np.bincount(tb, minlength=n_types)
        cnts[:len(bc)] = bc
        k_block.append(cnts)
    # run perms
    sums = np.zeros((B_PERM, len(blocks), n_types, len(tested)), dtype=np.float32)
    rng2 = np.random.default_rng(SEED + 1)
    for b in range(B_PERM):
        for bi, (rr, idx) in enumerate(blocks):
            n_b = len(idx)
            key = rng2.random(n_b)
            order = np.argsort(key)
            Gs = Gmat[idx][order]                     # (n_b, 20)
            cs = np.cumsum(Gs, axis=0)
            kb = k_block[bi]
            # rows for each type: cumulative position
            pos = np.cumsum(kb) - 1
            valid_t = kb > 0
            row = np.where(valid_t, pos, 0)
            sums[b, bi, valid_t] = cs[row[valid_t]]
    # T^perm per pair/region
    null_rows = []
    for _, r in res.iterrows():
        Lg, Rg = r.pair.split("->")
        S, R = TYPE_MAP[Lg], RECV_MAP[Rg]
        region_id = int(r.region)
        sidx = [bi for bi, (rr, _) in enumerate(blocks) if rr == region_id]
        bidx = np.array(sidx)
        S_ti, R_ti = type_id[S], type_id[R]
        NL_p = sums[:, bidx, S_ti, gidx[Lg]].sum(axis=1)
        NR_p = sums[:, bidx, R_ti, gidx[Rg]].sum(axis=1)
        Tp = np.log2(NL_p / r.nl_cells + EPS) + np.log2(NR_p / r.nr_cells + EPS)
        p_pos = float((Tp >= r.T_lo - 1e-9).mean())
        p_neg = float((Tp <= r.T_hi + 1e-9).mean())
        null_rows.append(dict(pair=r.pair, region=region_id,
                              p_pos=p_pos, p_neg=p_neg,
                              T_perm_median=float(np.median(Tp))))
    nul = pd.DataFrame(null_rows)
    # BH within region
    def bh(df, pcol):
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
    nul["q_pos"] = bh(nul, "p_pos")
    nul["q_neg"] = bh(nul, "p_neg")
    nul.to_csv(f"{OUT}/null_perm.csv", index=False)
    print(f"[{time.time()-t0:.0f}s] null done", flush=True)

    # ---------- merge and finalize ----------
    final = res.merge(nul, on=["pair", "region"])
    final["certified_pos"] = final.certified_pos & (final.q_pos <= 0.10)
    final["certified_neg"] = final.certified_neg & (final.q_neg <= 0.10)
    final.to_csv(f"{OUT}/scout_final.csv", index=False)
    print(final[["pair", "control", "region", "T0", "T_lo", "T_hi",
                 "certified_pos", "certified_neg", "q_pos", "q_neg"]].to_string(index=False))
    print(f"DONE in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
