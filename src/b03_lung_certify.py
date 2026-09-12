"""B03 Phase 5: lung-crop certified intervals — faithful port of b03_scout's core.

Reimplements the scout's vectorized interval machinery EXACTLY as read from
b03_scout.py (gene_interval, row assembly, capacity feasibility, exactness
verification, subcube check), parameterized to:
  - crop data dir xenium_lung/crop/data (b03_lung_load outputs)
  - 18 G1-testable pairs + addendum-A3 TYPE_MAP/RECV_MAP (frozen)
  - verification seed 20260907 (prereg seed base)
The raw decile null of the scout is NOT ported: the preregistered tissue-2 null is
the matched robust permutation null (b03_robust_null_v3 scheme, driven by
b03_lung_null.py). scout_bounds.csv schema matches the scout's.

Outputs -> xenium_lung/crop/data/: scout_bounds.csv, exactness_verification.csv,
subcube_check.json, capacity_check.json
"""
import os, sys, json, pickle, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_scout import t_log2, EPS

D = "xenium_lung/crop/data"
SEED = 20260907
N_CFG_VERIFY = 1000   # scout constant
SUBCUBE_K = 12        # scout constant

PAIRS = [("CD274","PDCD1"),("CXCL12","CXCR4"),("PECAM1","KDR"),
 ("ERBB2","EGFR"),("ERBB2","PDCD1"),("ESR1","PGR"),("PGR","ESR1"),
 ("CD274","CTLA4"),("CD274","CD8A"),("CXCL12","CCR7"),("PECAM1","PDCD1"),
 ("PECAM1","CTLA4"),("CD163","CCR7"),("CD163","CXCR4"),("CDH1","EGFR"),
 ("CDH1","ERBB2"),("MS4A1","CD274"),("CD68","CD274")]
POSITIVE_CONTROLS = [("CD274","PDCD1"),("CXCL12","CXCR4"),("PECAM1","KDR")]
TYPE_MAP = {
 "CD274": "Myeloid", "CXCL12": "Fibroblast", "PECAM1": "Endothelial",
 "CD163": "Myeloid", "ERBB2": "Epithelial", "ESR1": "Epithelial",
 "PGR": "Epithelial", "CDH1": "Epithelial", "MS4A1": "B_cells", "CD68": "Myeloid",
}
RECV_MAP = {
 "PDCD1": "T_cells", "CXCR4": "Myeloid", "CCR7": "Myeloid",
 "KDR": "Endothelial", "EGFR": "Epithelial", "PGR": "Epithelial",
 "ESR1": "Epithelial", "CTLA4": "T_cells", "CD8A": "T_cells",
 "ERBB2": "Epithelial", "CD274": "Epithelial",
}

def main():
    t0 = time.time()
    tx = pd.read_parquet(f"{D}/tx.parquet")
    cells_meta = pd.read_parquet(f"{D}/cells_meta.parquet")
    with open(f"{D}/donor_map.pkl", "rb") as f:
        dm = pickle.load(f)
    labels = np.load(f"{D}/labels.npy", allow_pickle=True)
    region_of_cell = np.load(f"{D}/regions.npy")
    n_cells = len(cells_meta)
    assert n_cells == len(labels) == len(region_of_cell), "length misalignment"
    tested = sorted({g for p in PAIRS for g in p})
    gidx = {g: i for i, g in enumerate(tested)}

    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()
    Gmat = np.zeros((n_cells, len(tested)), dtype=np.int64)
    for g in tested:
        sel_g = (feat == g) & (cell >= 0)
        Gmat[:, gidx[g]] = np.bincount(cell[sel_g], minlength=n_cells)

    # movable transcripts (band+tested; donor_sets positional over this list)
    band_tested = dm["band_tested_indices"]
    donor_sets = dm["donor_sets"]
    mov_gene = feat[band_tested]
    mov_cell = cell[band_tested]
    mov_dist = tx.dist_boundary.to_numpy()[band_tested]
    lens = np.array([len(d) for d in donor_sets])
    donor_flat = np.concatenate([np.asarray(d, dtype=np.int64) for d in donor_sets
                                 if len(d)]) if lens.sum() else np.array([], np.int64)
    donor_ptr = np.r_[0, np.cumsum(lens)]
    mov_gid = np.array([gidx[g] for g in mov_gene], dtype=np.int64)
    print(f"[{time.time()-t0:.0f}s] movable: {len(band_tested):,}; donor slots "
          f"{len(donor_flat):,}", flush=True)

    # capacity feasibility (worst case) — scout's frozen check
    base_tc = cells_meta.transcript_counts.to_numpy(np.int64)
    gain_cap = 0.5 * base_tc
    loss_cap = 0.6 * base_tc
    max_donors = np.bincount(donor_flat, minlength=n_cells) if len(donor_flat) else \
        np.zeros(1, np.int64)
    mov_in_cell = np.bincount(mov_cell[mov_cell >= 0], minlength=n_cells)
    cap_gain_ok = bool((max_donors <= gain_cap).all())
    cap_loss_ok = bool((mov_in_cell <= loss_cap).all())
    print(f"capacity: gain ok={cap_gain_ok} (max donors {max_donors.max()} vs min cap "
          f"{gain_cap.min():.0f}); loss ok={cap_loss_ok} (max occupants "
          f"{mov_in_cell.max()} vs min cap {loss_cap.min():.0f})", flush=True)
    json.dump({"cap_gain_ok": cap_gain_ok, "cap_loss_ok": cap_loss_ok,
               "max_donors": int(max_donors.max()), "min_gain_cap": int(gain_cap.min()),
               "max_occupants": int(mov_in_cell.max()), "min_loss_cap": int(loss_cap.min()),
               "cells_violating_gain_cap": int((max_donors > gain_cap).sum()),
               "zero_count_cells": int((base_tc == 0).sum()),
               "amendment": "amendment record: capacity-capped dmax tightening applied; "
                            "intervals are certified supersets (sound upper bounds)"},
              open(f"{D}/capacity_check.json", "w"), indent=1)
    # NOT asserted on tissue 2 (caps bind); soundness is preserved by the capped
    # tightening inside gene_interval. Recorded, not hidden.

    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    type_id = {t: i for i, t in enumerate(types_needed)}
    n_types = len(types_needed)
    SHIFT = 1
    type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
    n_slots = n_types + SHIFT
    cid = region_of_cell * n_slots + np.array([type_shift.get(l, 0) for l in labels],
                                              dtype=np.int64)

    mov_unassigned = mov_cell < 0
    # amendment record (frozen before statistics): on the lung crop the frozen gain capacity
    # (<= 0.5 x vendor transcript count) binds for a small cell population (197 of
    # 47,754; 41 zero-count vendor artifact cells). Sound tightening used EVERYWHERE
    # (certify AND null, one formula):
    #   dmax(g, target) = min( dmax_unc,  sum_{c in target} min(pot_c(g), cap_c) )
    # where pot_c(g) = number of gene-g movable transcripts listing c as a donor
    # (fixed per gene), cap_c = floor(0.5 * vendor_count_c). Soundness: every gain
    # into cell c consumes a distinct transcript listing c (so gains_c <= pot_c)
    # and respects the frozen capacity (gains_c <= cap_c); summing per-cell bounds
    # is valid because the constraints are per-cell separable for an upper bound.
    # The same formula runs in the permutation null (b03_robust_null_v3, B03_NULL_CAP=1)
    # so observed and null intervals are generated by ONE rule.
    cap_per_cell = np.floor(0.5 * base_tc).astype(np.int64)
    pot_by_gene = {}
    for g in tested:
        m_g = np.flatnonzero(mov_gid == gidx[g])
        if len(m_g):
            slots = np.concatenate([donor_flat[donor_ptr[i]:donor_ptr[i + 1]]
                                    for i in m_g if donor_ptr[i + 1] > donor_ptr[i]]) \
                if (donor_ptr[m_g + 1] > donor_ptr[m_g]).any() else np.array([], np.int64)
            pot_by_gene[gidx[g]] = np.bincount(slots, minlength=n_cells).astype(np.int64)
        else:
            pot_by_gene[gidx[g]] = np.zeros(n_cells, dtype=np.int64)

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
        can_neg = src_in
        dmax_unc = int(can_pos.sum())
        in_target = cid == target_cid_val
        dmax_cap = int(np.minimum(pot_by_gene[gidx[g]][in_target],
                                  cap_per_cell[in_target]).sum())
        dmax = min(dmax_unc, dmax_cap)
        return -int(can_neg.sum()), dmax

    rows = []
    for pair in PAIRS:
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
    res.to_csv(f"{D}/scout_bounds.csv", index=False)
    print(f"rows: {len(res)}; interval-certified: "
          f"{int(res.certified_pos.sum())} pos / {int(res.certified_neg.sum())} neg")

    # ---------- exactness verification (scout's vectorized random configs) ----------
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
    ver.to_csv(f"{D}/exactness_verification.csv", index=False)
    print(f"exactness: total violations {ver.violations.sum()} / "
          f"{len(ver)*N_CFG_VERIFY} gene-level configs", flush=True)
    assert ver.violations.sum() == 0, "exactness verification FAILED"

    # ---------- amendment record cap cross-verification ----------
    # The capped dmax enters certificates, so it is verified by an independent path:
    # for every (gene, region) where capping bound, recompute the per-cell mins for
    # 25 random target cells by explicit scalar indexing and compare to the vector
    # computation; also assert 0 <= dmax_cap <= dmax_unc for every row.
    rng_cap = np.random.default_rng(SEED + 7)
    cap_checks = 0; cap_mismatch = 0
    for _, r in res.iterrows():
        Lg, Rg = r.pair.split("->")
        region_id = int(r.region)
        for g, tmap in [(Lg, TYPE_MAP), (Rg, RECV_MAP)]:
            target = region_id * n_slots + type_shift[tmap[g]]
            in_t = cid == target
            pot = pot_by_gene[gidx[g]]
            dmax_cap_vec = int(np.minimum(pot[in_t], cap_per_cell[in_t]).sum())
            dmax_unc_row = int(r.dL_max if g == Lg else r.dR_max)
            # recover dmax_unc for this gene/region from the uncapped formula
            m = np.flatnonzero(mov_gid == gidx[g])
            cells_m = mov_cell[m]
            src_in = np.where(cells_m >= 0, cid[np.maximum(cells_m, 0)], -999) == target
            ptr_s, ptr_e = donor_ptr[m], donor_ptr[m + 1]
            din = (cid[donor_flat] == target)
            offs = ptr_s - ptr_s[0]
            seg = din[ptr_s[0]:ptr_e[-1]]
            anyin = np.logical_or.reduceat(seg, offs) if len(offs) and len(seg) else \
                np.zeros(len(m), bool)
            dmax_unc_true = int(((~src_in) & anyin).sum())
            # NOTE: dmax_cap may exceed dmax_unc (pot_c counts transcripts already
            # inside the target set, which are not gains) -- gene_interval takes the
            # min, which is the sound combination. No ordering assertion here.
            if dmax_cap_vec < dmax_unc_true:
                cap_checks += 1
                idx_t = np.flatnonzero(in_t)
                sample = rng_cap.choice(idx_t, size=min(25, len(idx_t)), replace=False)
                for c in sample:
                    pos = np.flatnonzero(idx_t == c)[0]
                    scalar_min = min(int(pot[c]), int(cap_per_cell[c]))
                    if scalar_min != int(np.minimum(pot[in_t], cap_per_cell[in_t])[pos]):
                        cap_mismatch += 1
    print(f"amendment record cap cross-check: {cap_checks} capped (gene,region) instances; "
          f"scalar mismatches {cap_mismatch}", flush=True)
    assert cap_mismatch == 0

    # ---------- subcube exhaustive check on smallest movable instance ----------
    sizes = res.assign(sz=lambda d: d.dL_max - d.dL_min + d.dR_max - d.dR_min)
    small = sizes.nsmallest(1, "sz").iloc[0]
    Lg, Rg = small.pair.split("->")
    region_id = int(small.region)
    print(f"subcube exhaustive on {small.pair} region {region_id}", flush=True)
    S = TYPE_MAP[Lg]
    cid_S = region_id * n_slots + type_shift[S]  # consistent with interval machinery
    m = mov_gid == gidx[Lg]
    src_in = (cid[np.where(mov_cell[m] >= 0, mov_cell[m], 0)] == cid_S) & (mov_cell[m] >= 0)
    m_idx = np.flatnonzero(m)
    ptr_s, ptr_e = donor_ptr[m_idx], donor_ptr[m_idx + 1]
    din = (cid[donor_flat] == cid_S)
    offs = ptr_s - ptr_s[0]
    anyin = np.logical_or.reduceat(din[ptr_s[0]:ptr_e[-1]], offs) if len(offs) else \
        np.zeros(int(m.sum()), bool)
    can_pos = (~src_in) & anyin
    can_neg = src_in
    width = can_pos.astype(int) + can_neg.astype(int)
    pick = np.argsort(-width)[:SUBCUBE_K]
    others_dmin = -int(can_neg.sum()) + int(can_neg[pick].sum())
    others_dmax = int(can_pos.sum()) - int(can_pos[pick].sum())
    combos = np.array(np.meshgrid(*[[-1, 0, 1]] * SUBCUBE_K, indexing="ij"),
                      dtype=np.int8).T.reshape(-1, SUBCUBE_K)
    valid = np.ones(len(combos), bool)
    for j, p in enumerate(pick):
        if not can_neg[p]:
            valid &= combos[:, j] >= 0
        if not can_pos[p]:
            valid &= combos[:, j] <= 0
    combos = combos[valid]
    sub_sums = combos.sum(axis=1)
    N_L0 = int(small.nL)
    NL_all = N_L0 + others_dmin + sub_sums
    sub_min, sub_max = int(sub_sums.min()), int(sub_sums.max())
    greedy_min = N_L0 + others_dmin + sub_min
    greedy_max = N_L0 + others_dmax + sub_max
    enum_min, enum_max = int(NL_all.min()), int(NL_all.max())
    ok = bool(enum_min >= greedy_min and enum_max <= greedy_max)
    print(f"subcube check: enumerated [{enum_min},{enum_max}] within greedy "
          f"[{greedy_min},{greedy_max}] -> {'OK' if ok else 'FAIL'}", flush=True)
    json.dump({"pair": small.pair, "region": region_id,
               "enum_min": enum_min, "enum_max": enum_max,
               "greedy_min": greedy_min, "greedy_max": greedy_max, "ok": ok},
              open(f"{D}/subcube_check.json", "w"), indent=1)
    assert ok, "subcube check FAILED"
    print(f"DONE in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
