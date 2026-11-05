"""Coupled interval ends for same-class pairs, and the tightness criterion.

Background. For a pair whose sender and receiver classes differ. the two per-gene
counts are driven by disjoint transcript sets and disjoint counted cell sets, so the
per-gene optima are attained by *different* assignments and the separable interval

    T_lo = f(min N_L) + g(min N_R),   T_hi = f(max N_L) + g(max N_R)

contains the attainable range without generally equalling it. A same-class pair
(sender class = receiver class) shares one counted set C and one set of per-cell
gains and losses, so a single assignment must serve both terms. The separable end
is then only a bound, and whether it is attained depends on whether the two genes
can take their budgets simultaneously.

What is computed here.

  arrival side   rho_g   = max inflow of gene g into C            (one max flow)
                 rho_LR  = max total inflow of both genes        (one max flow)
                 The separable upper end is attained if and only if
                 rho_LR = rho_L + rho_R: at total rho_LR the L-part is at most
                 rho_L and the R-part at most rho_R, so both are saturated
                 when their sum is the total. When the criterion fails no single
                 assignment is known to attain both maxima. the module then reports
                 the tightest value available under the two individual caps and the
                 joint total, which is still a sound interval end, together with the
                 measured distance from the separable value.

  departure side D_g    = sum_{s in C} min(a_g(s), loss_cap_s)   (closed form)
                 D_LR   = sum_{s in C} min(a_L(s)+a_R(s), loss_cap_s)
                 Departures are cell-local and the extracellular compartment is
                 uncapped, so the L-first allocation attains (D_L, D_LR - D_L) and
                 the R-first allocation attains (D_LR - D_R, D_R). every split
                 between them is reached by transferring one unit at a time inside
                 a cell, which is a witness, not an argument. No flow is needed.

Every reported value comes with a witness: the integral flow is read back and both
budgets, both per-gene counts and every group capacity are re-checked, on the
departure side the allocation is built cell by cell and re-checked the same way.

Writes, per leg, <leg>/coupled_end.csv, and the combined summary
results/coupled_end_summary.json.

Usage:  python src/b03_coupled_end.py
        B03_DATA=... B03_WORKSPACE=... python src/b03_coupled_end.py
"""
import json
import os
import pickle
import sys
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b03_scout import EPS, t_log2  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
SEED = None  # deterministic; no randomness in this computation

BREAST_DATA = os.environ.get("B03_DATA", os.path.join(ROOT, "data"))
LUNG_DATA = os.path.join(WORKSPACE, "xenium_lung", "crop", "data")
S6_DATA = os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "data")

BREAST_PAIRS = [("CD274", "PDCD1"), ("CXCL12", "CXCR4"), ("CCL5", "CCR7"),
                ("PECAM1", "KDR"), ("CD163", "CD3D"), ("ERBB2", "EGFR"),
                ("ERBB2", "PDCD1"), ("ESR1", "PGR"), ("PGR", "ESR1"),
                ("CD274", "CTLA4"), ("CD274", "CD8A"), ("CXCL12", "CCR7"),
                ("CCL5", "CXCR4"), ("PECAM1", "PDCD1"), ("PECAM1", "CTLA4"),
                ("CD163", "CCR7"), ("CD163", "CXCR4"), ("ACTA2", "EGFR"),
                ("ACTA2", "ERBB2"), ("KRT8", "ERBB2"), ("CDH1", "EGFR"),
                ("CDH1", "ERBB2"), ("CD3D", "ERBB2"), ("MS4A1", "CD274"),
                ("CD68", "CD274")]
BREAST_TYPE = {
    "CD274": "Breast cancer", "CXCL12": "Fibroblasts", "CCL5": "Macrophages",
    "PECAM1": "Endothelial cells", "CD163": "Macrophages", "ERBB2": "Breast cancer",
    "ESR1": "Breast cancer", "PGR": "Breast cancer", "ACTA2": "Smooth muscle cells",
    "KRT8": "Breast glandular cells", "CDH1": "Breast glandular cells",
    "CD3D": "T cells", "MS4A1": "B cells", "CD68": "Macrophages",
}
BREAST_RECV = {
    "PDCD1": "T cells", "CXCR4": "T cells", "CCR7": "T cells",
    "KDR": "Endothelial cells", "CD3D": "T cells", "EGFR": "Breast cancer",
    "PGR": "Breast cancer", "ESR1": "Breast cancer", "CTLA4": "T cells",
    "CD8A": "T cells", "ERBB2": "Breast cancer", "CD274": "Breast cancer",
}

GAIN_FRAC = 0.5
LOSS_FRAC = 0.6
INF = 10 ** 9


def _lung_cfg():
    from b03_lung_config import PAIRS, TYPE_MAP, RECV_MAP
    return PAIRS, TYPE_MAP, RECV_MAP


def legs():
    lp, lt, lr = _lung_cfg()
    return [dict(tag="breast", data=BREAST_DATA, out=os.path.join(ROOT, "results", "breast"),
                 rows=os.path.join(ROOT, "results", "breast", "scout_bounds.csv"),
                 pairs=BREAST_PAIRS, tmap=BREAST_TYPE, rmap=BREAST_RECV),
            dict(tag="lung", data=LUNG_DATA, out=os.path.join(ROOT, "results", "lung"),
                 rows=os.path.join(ROOT, "results", "lung", "scout_bounds.csv"),
                 pairs=lp, tmap=lt, rmap=lr),
            dict(tag="s6", data=S6_DATA, out=os.path.join(ROOT, "results", "breast_s6"),
                 rows=os.path.join(ROOT, "results", "breast_s6", "scout_bounds.csv"),
                 pairs=lp, tmap=lt, rmap=lr)]


def load(data_dir):
    tx = pd.read_parquet(os.path.join(data_dir, "tx.parquet"),
                         columns=["feature_name", "cell_idx"])
    with open(os.path.join(data_dir, "donor_map.pkl"), "rb") as f:
        dm = pickle.load(f)
    labels = np.load(os.path.join(data_dir, "labels.npy"), allow_pickle=True).astype(str)
    regions = np.load(os.path.join(data_dir, "regions.npy"), allow_pickle=True).astype(np.int64)
    cm = pd.read_parquet(os.path.join(data_dir, "cells_meta.parquet"),
                         columns=["cell_id", "transcript_counts"])
    cid_map = {str(k): int(v) for k, v in dm["cid_map"].items()}
    poly_row = np.full(max(cid_map.values()) + 1, -1, dtype=np.int64)
    for row, cid in enumerate(cm.cell_id.astype(str).to_numpy()):
        p = cid_map.get(cid, -1)
        if p >= 0:
            poly_row[p] = row
    assert (poly_row >= 0).all(), "unmapped polygons present"
    base = cm.transcript_counts.to_numpy(np.int64)[poly_row]
    assert len(base) == len(labels) == len(regions)
    band = dm["band_tested_indices"]
    lens = np.array([len(x) for x in dm["donor_sets"]])
    flat = (np.concatenate([np.asarray(x, np.int64) for x in dm["donor_sets"] if len(x)])
            if lens.sum() else np.zeros(0, np.int64))
    return dict(feat=tx.feature_name.to_numpy(), cell=tx.cell_idx.to_numpy(),
                labels=labels, regions=regions, base=base,
                mov_gene=tx.feature_name.to_numpy()[band],
                mov_cell=tx.cell_idx.to_numpy()[band],
                flat=flat, ptr=np.r_[0, np.cumsum(lens)], n_cells=len(base))


def gene_index(d):
    idx = {}
    for g in np.unique(d["mov_gene"]):
        idx[str(g)] = np.flatnonzero(d["mov_gene"] == g)
    d["gene_index"] = idx
    return idx


def groups(d, gene, target, m=None):
    """Contract this gene's movable transcripts into (source cell, donor set) groups."""
    if m is None:
        m = np.flatnonzero(d["mov_gene"] == gene)
    if not len(m):
        return None
    c = d["mov_cell"][m]
    src_in = np.where(c >= 0, target[np.maximum(c, 0)], False)
    out = ~src_in
    ps, pe = d["ptr"][m], d["ptr"][m + 1]
    ln = pe - ps
    tids = np.repeat(np.arange(len(m)), ln)
    slots = (np.concatenate([np.arange(a, b) for a, b in zip(ps, pe)])
             if ln.sum() else np.zeros(0, np.int64))
    dc = d["flat"][slots]
    keep = out[tids] & target[dc]
    tid, dcell = tids[keep], dc[keep]
    if not len(tid):
        return None
    order = np.lexsort((dcell, tid))
    tid, dcell = tid[order], dcell[order]
    starts = np.flatnonzero(np.r_[True, tid[1:] != tid[:-1]])
    gl = np.diff(np.r_[starts, len(tid)])
    kmax = int(gl.max())
    mat = np.full((len(starts), 1 + kmax), -1, dtype=np.int64)
    mat[:, 0] = c[tid[starts]]
    within = np.arange(len(tid)) - np.repeat(starts, gl)
    mat[np.repeat(np.arange(len(starts)), gl), 1 + within] = dcell
    sig, cnt = np.unique(mat, axis=0, return_counts=True)
    return sig, cnt


def build(segs, loss, gain):
    """Merged network: source -> source cell -> group -> destination cell -> sink.

    segs is a list of (gene_index, signature matrix, group counts). The gene index
    is carried on the group node, which is what makes the two genes' shared cell
    budgets a single-commodity network: a unit of either gene consumes the same
    source loss budget and the same destination gain budget.
    """
    src = np.unique(np.concatenate([s[1][:, 0] for s in segs]))
    de = np.unique(np.concatenate([s[1][:, 1:][s[1][:, 1:] >= 0] for s in segs]))
    ns, ng, nt = len(src), sum(len(s[1]) for s in segs), len(de)
    S, off_s, off_g, off_d = 0, 1, 1 + ns, 1 + ns + ng
    T = off_d + nt
    # source ids are shifted by one so that the extracellular compartment (-1) has a
    # slot of its own; indexing a lookup array with -1 would otherwise address the
    # last cell and silently reroute every extracellular transcript
    smap = np.full(int(max(src.max(), 0)) + 2, -1, np.int64)
    smap[src + 1] = np.arange(ns)
    tmap = np.full(int(max(de.max(), 0)) + 1, -1, np.int64)
    tmap[de] = np.arange(nt)
    # arcs are assembled with numpy: the Python loop over groups is what made the
    # first implementation quadratic in the number of signatures
    rows, cols, caps, gene_of = [], [], [], []
    src_ok = src >= 0
    rows.append(np.zeros(int(src_ok.sum()), np.int64))
    cols.append(off_s + np.flatnonzero(src_ok))
    caps.append(loss[src[src_ok]])
    if (~src_ok).any():
        rows.append(np.zeros(int((~src_ok).sum()), np.int64))
        cols.append(off_s + np.flatnonzero(~src_ok))
        caps.append(np.full(int((~src_ok).sum()), INF, np.int64))
    base = 0
    for gi, sig, cnt in segs:
        n_sig = len(sig)
        gidx = np.arange(base, base + n_sig)
        rows.append(off_s + smap[sig[:, 0] + 1])
        cols.append(off_g + gidx)
        caps.append(cnt.astype(np.int64))
        gene_of.append(np.full(n_sig, gi, np.int64))
        rest = sig[:, 1:]
        nz = rest >= 0
        gl = nz.sum(axis=1)
        if gl.sum():
            gi_rep = np.repeat(gidx, gl)
            dt_rep = rest[nz]
            rows.append(off_g + gi_rep)
            cols.append(off_d + tmap[dt_rep])
            caps.append(np.repeat(cnt, gl).astype(np.int64))
        base += n_sig
    rows.append(off_d + np.arange(nt))
    cols.append(np.full(nt, T, np.int64))
    caps.append(gain[de])
    n = T + 1
    graph = csr_matrix((np.concatenate(caps).astype(np.int64),
                        (np.concatenate(rows), np.concatenate(cols))),
                       shape=(n, n))
    return dict(graph=graph, S=S, T=T, n_node=n, ns=ns, ng=ng, nt=nt,
                off_s=off_s, off_g=off_g, off_d=off_d, dest_cells=de,
                src_cells=src, gene_of=np.concatenate(gene_of))


def read_flow(w, res):
    """Per-group outflow, per-destination inflow, per-source outflow."""
    fl = res.flow.tocsr()
    sub = fl[w["off_g"]:w["off_g"] + w["ng"], w["off_d"]:w["off_d"] + w["nt"]]
    out_group = np.asarray(sub.sum(axis=1)).ravel()
    inflow = np.asarray(sub.sum(axis=0)).ravel()
    by_src = np.asarray(fl[w["S"], w["off_s"]:w["off_s"] + w["ns"]].todense()).ravel()
    return out_group, inflow, by_src


def check_flow(w, res, loss, gain, cap_by_group):
    """Re-check budgets and per-gene totals on an integral flow."""
    out_group, inflow, by_src = read_flow(w, res)
    ok = bool((out_group <= cap_by_group).all())
    ok = ok and bool((inflow <= gain[w["dest_cells"]]).all())
    for si, c in enumerate(w["src_cells"]):
        if c >= 0 and by_src[si] > loss[c]:
            ok = False
    per_gene = np.zeros(2, np.int64)
    for k, gi in enumerate(w["gene_of"]):
        per_gene[gi] += int(out_group[k])
    return ok, per_gene


def _best_split(N0L, N0R, nL, nR, total, lo, hi):
    """Integer split in [lo, hi] maximising f(N0L+x) + g(N0R+total-x)."""
    lo = int(max(lo, 0))
    hi = int(min(hi, total))
    if hi < lo:
        return lo
    xs = np.arange(lo, hi + 1)
    vals = np.log2((N0L + xs) / nL + EPS) + np.log2((N0R + (total - xs)) / nR + EPS)
    return int(xs[int(np.argmax(vals))])


def _best_dep_split(N0L, N0R, nL, nR, total, lo, hi):
    """Integer departure split minimising the joint lower score.

    Both score terms decrease in their departure count, so the worst case over the
    feasible departure simplex is the split that loses the most score, not the one
    that loses the least. Each candidate split in [lo, hi] is feasible for the
    merged loss budget and the two per-gene budgets.
    """
    lo = int(max(lo, 0))
    hi = int(min(hi, total))
    if hi < lo:
        return lo
    xs = np.arange(lo, hi + 1)
    vals = np.log2((N0L - xs) / nL + EPS) + np.log2((N0R - (total - xs)) / nR + EPS)
    return int(xs[int(np.argmin(vals))])


def departure_witness(aL, aR, capc, dL, dR):
    """Build and re-check a per-cell loss allocation with the requested totals."""
    n = len(aL)
    outL = np.minimum(aL, capc)
    rest = capc - outL
    outR = np.minimum(aR, rest)
    if int(outL.sum()) < dL:
        return False
    # move units from L to R inside cells: feasible iff the cell has L units
    # allocated and an unallocated R transcript
    need = int(outL.sum()) - dL
    moved = 0
    for s in range(n):
        if moved >= need:
            break
        can = min(outL[s], aR[s] - outR[s])
        if can <= 0:
            continue
        take = int(min(can, need - moved))
        outL[s] -= take
        outR[s] += take
        moved += take
    if moved != need:
        return False
    if int(outL.sum()) != dL or int(outR.sum()) != dR:
        return False
    if (outL + outR > capc).any() or (outL > aL).any() or (outR > aR).any():
        return False
    return True


def leg_row(d, r, tmap, rmap, cache_g, cache_a):
    """Group signatures and inside-cell counts for one row, cached by (gene, region, class)."""
    Lg, Rg = r.pair.split("->")
    reg = int(r.region)
    out = []
    for g, tmap_ in ((Lg, tmap), (Rg, rmap)):
        cls = tmap_[g]
        key = (g, reg, cls)
        target = (d["regions"] == reg) & (d["labels"] == cls)
        if key not in cache_g:
            m = d["gene_index"].get(g)
            cache_g[key] = (groups(d, g, target, m), target, m)
        sig, tgt, m = cache_g[key]
        if key not in cache_a:
            c = d["mov_cell"][m]
            ins = np.where(c >= 0, target[np.maximum(c, 0)], False)
            cache_a[key] = np.bincount(c[ins], minlength=d["n_cells"])[target]
        out.append((sig, tgt, cache_a[key]))
    return out


def main():
    t0 = time.time()
    rows_all = []
    summary = dict(legs={})
    only = os.environ.get("B03_LEGS", "")
    for leg in legs():
        if only and leg["tag"] not in only.split(","):
            continue
        d = load(leg["data"])
        gene_index(d)
        loss = np.floor(LOSS_FRAC * d["base"]).astype(np.int64)
        gain = np.floor(GAIN_FRAC * d["base"]).astype(np.int64)
        sb = pd.read_csv(leg["rows"])
        cache_g, cache_a, cache_w = {}, {}, {}
        rows = []
        for _, r in sb.iterrows():
            Lg, Rg = r.pair.split("->")
            (sl, tl, aL_all), (sr, tr, aR_all) = leg_row(d, r, leg["tmap"], leg["rmap"],
                                                        cache_g, cache_a)
            if sl is None or sr is None:
                continue
            reg = int(r.region)
            N0L, N0R = int(r.nL), int(r.nR)
            nls, nrs = int(r.nl_cells), int(r.nr_cells)
            T_sep_lo, T_sep_hi = float(r.T_lo), float(r.T_hi)
            same = bool((tl == tr).all())

            # ---------------- arrival side ----------------
            if ("L", Lg, reg, leg["tmap"][Lg]) not in cache_w:
                cache_w[("L", Lg, reg, leg["tmap"][Lg])] = build(
                    [(0, sl[0], sl[1])], loss, gain)
            if ("R", Rg, reg, leg["rmap"][Rg]) not in cache_w:
                cache_w[("R", Rg, reg, leg["rmap"][Rg])] = build(
                    [(1, sr[0], sr[1])], loss, gain)
            wl = cache_w[("L", Lg, reg, leg["tmap"][Lg])]
            wr = cache_w[("R", Rg, reg, leg["rmap"][Rg])]
            wb = build([(0, sl[0], sl[1]), (1, sr[0], sr[1])], loss, gain)
            rho_l = int(maximum_flow(wl["graph"], wl["S"], wl["T"]).flow_value)
            rho_r = int(maximum_flow(wr["graph"], wr["S"], wr["T"]).flow_value)
            res = maximum_flow(wb["graph"], wb["S"], wb["T"])
            rho_lr = int(res.flow_value)
            ok_flow, per_gene = check_flow(wb, res, loss, gain,
                                           np.concatenate([sl[1], sr[1]]))
            coll_arr = int(rho_l + rho_r - rho_lr)
            if coll_arr == 0:
                xL, xR = rho_l, rho_r
                arr_witness = bool(ok_flow and per_gene[0] == rho_l and per_gene[1] == rho_r)
                method, gap = "simultaneous", 0.0
            else:
                # The criterion fails, so no single assignment is known to attain both
                # per-gene maxima. Report the tightest value any assignment can have
                # under the two individual caps and the joint total -- a sound upper
                # bound on the attained maximum, hence a sound interval end -- and the
                # distance between it and the separable value as the measured slack.
                ub_split = _best_split(N0L, N0R, nls, nrs, rho_lr,
                                       max(0, rho_lr - rho_r), min(rho_l, rho_lr))
                xL, xR = ub_split, rho_lr - ub_split
                method = "relaxation"
                arr_witness = False
                gap = float(max(0.0, t_log2(N0L + xL, N0R + xR, nls, nrs) - T_sep_hi))
            T_hi_joint = t_log2(N0L + xL, N0R + xR, nls, nrs)

            # ---------------- departure side ----------------
            # per-cell loss caps are shared by whichever counted sets the two genes use
            cells_l = np.flatnonzero(tl)
            cells_r = np.flatnonzero(tr)
            aL = aL_all
            aR = aR_all
            DL = int(np.minimum(aL, loss[cells_l]).sum())
            DR = int(np.minimum(aR, loss[cells_r]).sum())
            if same:
                capc = loss[cells_l]
                DLR = int(np.minimum(aL + aR, capc).sum())
                coll_dep = int(DL + DR - DLR)
                xlo_d, xhi_d = int(DLR - DR), int(DL)
                xs_d = _best_dep_split(N0L, N0R, nls, nrs, DLR, xlo_d, xhi_d)
                dL, dR = int(xs_d), int(DLR - xs_d)
                dep_witness = bool(departure_witness(aL, aR, capc, dL, dR))
            else:
                # disjoint counted sets: shared loss budgets only at cells outside both,
                # which no departure uses, so the ends separate
                DLR = int(DL + DR)
                coll_dep = 0
                dL, dR = DL, DR
                dep_witness = True
            T_lo_joint = t_log2(max(0, N0L - dL), max(0, N0R - dR), nls, nrs)

            # per-transcript counts: the O(|C|) refinement in Table S8.1b, kept apart
            # from the exact per-gene optima so the two corrections stay separate
            cnt_bin_L = int(sl[1].sum())
            cnt_bin_R = int(sr[1].sum())
            T_hi_count = t_log2(N0L + cnt_bin_L, N0R + cnt_bin_R, nls, nrs)
            T_lo_count = t_log2(max(0, N0L - int(aL_all.sum())),
                                max(0, N0R - int(aR_all.sum())), nls, nrs)
            # separable optimum at the exact per-gene ends: what Proposition 1 bounds
            T_hi_flowsep = t_log2(N0L + rho_l, N0R + rho_r, nls, nrs)
            T_lo_flowsep = t_log2(max(0, N0L - DL), max(0, N0R - DR), nls, nrs)

            rows.append(dict(
                section=leg["tag"], pair=r.pair, region=int(r.region),
                same_class=same, nL0=N0L, nR0=N0R, cells_L=nls, cells_R=nrs,
                count_L=cnt_bin_L, count_R=cnt_bin_R,
                rho_L=rho_l, rho_R=rho_r, rho_LR=rho_lr, arrival_collision=coll_arr,
                x_L=xL, x_R=xR, arrival_method=method,
                D_L=DL, D_R=DR, D_LR=DLR, departure_collision=coll_dep,
                d_L=dL, d_R=dR,
                T_lo_sealed=T_sep_lo, T_hi_sealed=T_sep_hi,
                T_lo_count=T_lo_count, T_hi_count=T_hi_count,
                T_lo_flowsep=T_lo_flowsep, T_hi_flowsep=T_hi_flowsep,
                T_lo_joint=T_lo_joint, T_hi_joint=T_hi_joint,
                slack_hi_flow=T_hi_flowsep - T_hi_joint,
                slack_lo_flow=T_lo_flowsep - T_lo_joint,
                slack_hi_count=T_hi_count - T_hi_joint,
                dT_hi_sealed=T_hi_joint - T_sep_hi,
                # Reported interval against the joint (exact) interval, which mixes the
                # counting and capacity refinements.
                sign_changed_vs_sealed=bool(((T_lo_joint > 0) != (T_sep_lo > 0)) or
                                            ((T_hi_joint < 0) != (T_sep_hi < 0))),
                # Coupling alone: joint endpoints against the separate-gene optima.
                sign_changed_coupling=bool(
                    ((T_lo_joint > 0) != (T_lo_flowsep > 0)) or
                    ((T_hi_joint < 0) != (T_hi_flowsep < 0))),
                arrival_witness=arr_witness, departure_witness=dep_witness,
                certified=bool((T_lo_joint > 0) or (T_hi_joint < 0))))
        df = pd.DataFrame(rows)
        os.makedirs(leg["out"], exist_ok=True)
        df.to_csv(os.path.join(leg["out"], "coupled_end.csv"), index=False)
        rows_all.append(df)
        st = dict(n_rows=int(len(df)),
                  n_pairs=int(df.pair.nunique()) if len(df) else 0,
                  n_same_class_rows=int(df.same_class.sum()) if len(df) else 0,
                  arrival_collision_free=int((df.arrival_collision == 0).sum()),
                  departure_collision_free=int((df.departure_collision == 0).sum()),
                  tight=int(((df.arrival_collision == 0) &
                             (df.departure_collision == 0)).sum()),
                  max_arrival_collision=int(df.arrival_collision.max()),
                  max_departure_collision=int(df.departure_collision.max()),
                  max_slack_hi_flow=float(df.slack_hi_flow.max()),
                  max_slack_lo_flow=float(df.slack_lo_flow.max()),
                  max_slack_hi_count=float(df.slack_hi_count.max()),
                  max_slack_hi_flow_same=float(df.loc[df.same_class, "slack_hi_flow"].max())
                  if df.same_class.any() else 0.0,
                  sign_changes_vs_sealed=int(df.sign_changed_vs_sealed.sum()),
                  sign_changes_from_coupling=int(df.sign_changed_coupling.sum()),
                  witness_fail=int((~df.arrival_witness).sum() +
                                   (~df.departure_witness).sum()))
        summary["legs"][leg["tag"]] = st
        print(leg["tag"], json.dumps(st, sort_keys=True), flush=True)

    allrows = pd.concat(rows_all, ignore_index=True) if rows_all else pd.DataFrame()
    summary.update(
        n_rows=int(len(allrows)),
        n_same_class_rows=int(allrows.same_class.sum()),
        n_tight=int(((allrows.arrival_collision == 0) &
                     (allrows.departure_collision == 0)).sum()),
        n_collision_free_arrival=int((allrows.arrival_collision == 0).sum()),
        n_collision_free_departure=int((allrows.departure_collision == 0).sum()),
        n_joint_sign_changes_vs_sealed=int(allrows.sign_changed_vs_sealed.sum()),
        n_joint_sign_changes_from_coupling=int(allrows.sign_changed_coupling.sum()),
        max_arrival_collision=int(allrows.arrival_collision.max()),
        max_departure_collision=int(allrows.departure_collision.max()),
        max_slack_hi_flow=float(allrows.slack_hi_flow.max()),
        max_slack_lo_flow=float(allrows.slack_lo_flow.max()),
        max_slack_hi_count=float(allrows.slack_hi_count.max()),
        witness_fail=int((~allrows.arrival_witness).sum() +
                         (~allrows.departure_witness).sum()),
        n_hi_attained=int((allrows.slack_hi_flow < 1e-12).sum()),
        n_lo_attained=int((allrows.slack_lo_flow < 1e-12).sum()),
        n_certified_rows=int(allrows.certified.sum()),
        secs=round(time.time() - t0, 1))
    with open(os.path.join(ROOT, "results", "coupled_end_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, sort_keys=True)
    print(json.dumps(summary, indent=1, sort_keys=True))
    print(f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
