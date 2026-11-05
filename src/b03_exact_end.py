"""Exact interval ends by max flow on the contracted cell-level network.

Three upper bounds on the positive end appear in the artefacts:

  count  dmax_bin = #{movable g-transcripts outside C with a donor in C}       O(1) per gene
  capped dmax_cap = min(dmax_bin, sum_{d in C} min(pot_d, gain_cap_d))        O(|C|) per gene
  reported       = the vectorized donor-slot count in results/breast/scout_bounds.csv

pot_d counts donor slots of gene g that point to cell d, including slots of
transcripts already in C, so the capped bound refines the count only when the
per-cell gain caps bite. On breast they do not, and the two agree on all 150
gene-ends. The reported count reads donor slots through a vectorized reduction
that can credit one transcript once for each eligible donor, so it sits above
the count bound by at most 0.722 log2 on the 75 tests. All three ignore that a
transcript that would gain also consumes its source cell's loss budget. The
exact end is

  nu+ = max flow  source -> s (cap loss_cap_s) -> d (cap n_g(s,d)) -> sink (cap gain_cap_d)

over cells. Transcripts sharing (s, d) are interchangeable,
so parallel unit arcs merge into one arc of capacity equal to their count.  The
network therefore has O(|dC|) nodes and O(d|dC|) arcs, independent of the number
of molecules.  The negative end has the closed form dmin = -sum_{s in C} min(a_s,
loss_cap_s), because dropping to the extracellular compartment is always
admissible and uncapped.

Writes results/breast/exact_end.csv and results/breast/exact_end_summary.json.

Usage:  python src/b03_exact_end.py            (DATA from ../data)
        B03_DATA=/path/to/data python src/b03_exact_end.py
"""
import json
import math
import os
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

EPS = 0.5
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("B03_DATA", os.path.join(ROOT, "data"))
OUT = os.path.join(ROOT, "results", "breast")
SEED = None  # no randomness in this check

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

GAIN_FRAC = 0.5
LOSS_FRAC = 0.6


def t_log2(nL, nR, a, b, eps=EPS):
    return math.log2(nL / a + eps) + math.log2(nR / b + eps)


def load():
    tx = pd.read_parquet(os.path.join(DATA, "tx.parquet"),
                        columns=["feature_name", "cell_idx", "dist_boundary"])
    with open(os.path.join(DATA, "donor_map.pkl"), "rb") as f:
        import pickle
        dm = pickle.load(f)
    labels = np.load(os.path.join(DATA, "labels.npy"), allow_pickle=True).astype(str)
    regions = np.load(os.path.join(DATA, "regions.npy"), allow_pickle=True).astype(np.int64)
    # cells_meta.parquet is in barcode order; every per-cell array downstream is in
    # polygon-index order, so remap as b03_scout.py does before any cap is read.
    cm = pd.read_parquet(os.path.join(DATA, "cells_meta.parquet"),
                         columns=["cell_id", "transcript_counts"])
    cid_map = {str(k): int(v) for k, v in dm["cid_map"].items()}
    poly_row = np.full(max(cid_map.values()) + 1, -1, dtype=np.int64)
    for row, cid in enumerate(cm.cell_id.astype(str).to_numpy()):
        p = cid_map.get(cid, -1)
        if p >= 0:
            poly_row[p] = row
    assert (poly_row >= 0).all(), "unmapped polygons present"
    base_tc = cm.transcript_counts.to_numpy(np.int64)[poly_row]
    assert len(base_tc) == len(labels) == len(regions)
    return tx, dm, labels, regions, base_tc


def exact_positive(mov_idx, mov_cell, donor_ptr, donor_flat, target_mask,
                   loss_cap, gain_cap, n_cells):
    """Exact max inflow into the counted set.

    mov_idx    : positions (into the band+tested arrays) of this gene's movable transcripts
    Returns (nu_plus, dmax_sum_bound, node_count, arc_count).
    """
    cell = mov_cell[mov_idx]
    src_in = np.where(cell >= 0, target_mask[np.maximum(cell, 0)], False)
    out = ~src_in
    # donor slots of this gene, as (transcript position, slot) pairs
    ps, pe = donor_ptr[mov_idx], donor_ptr[mov_idx + 1]
    lens = pe - ps
    tid_all = np.repeat(np.arange(len(mov_idx)), lens)
    slot_all = np.concatenate([np.arange(a, b) for a, b in zip(ps, pe)]) if lens.sum() \
        else np.zeros(0, np.int64)
    dcell_all = donor_flat[slot_all]
    # a transcript contributes an arc only if it starts outside C and can enter C
    keep = out[tid_all] & target_mask[dcell_all]
    tid, dcell = tid_all[keep], dcell_all[keep]

    # One unit per transcript, not one per (transcript, donor) slot: a transcript with
    # two donors in C must choose one.  Transcripts that share a source cell and an
    # identical destination set are interchangeable, so they collapse to one group
    # node of capacity equal to their count.  This is the contraction that makes the
    # network size independent of the number of molecules.
    order = np.lexsort((dcell, tid))
    tid_s, d_s = tid[order], dcell[order]
    starts = np.flatnonzero(np.r_[True, tid_s[1:] != tid_s[:-1]]) if len(tid_s) else \
        np.zeros(0, np.int64)
    if len(starts):
        grp_len = np.diff(np.r_[starts, len(tid_s)])
        kmax = int(grp_len.max())
        mat = np.full((len(starts), 1 + kmax), -1, dtype=np.int64)
        mat[:, 0] = cell[tid_s[starts]]
        within = np.arange(len(tid_s)) - np.repeat(starts, grp_len)
        mat[np.repeat(np.arange(len(starts)), grp_len), 1 + within] = d_s
        gsig, ginv, gcnt = np.unique(mat, axis=0, return_inverse=True, return_counts=True)
    else:
        gsig = np.zeros((0, 2), np.int64)
        gcnt = np.zeros(0, np.int64)

    # node numbering: 0 = source, then source cells, then groups, then destination
    # cells, then sink
    src_cells = np.unique(gsig[:, 0]) if len(gsig) else np.zeros(0, np.int64)
    dests = np.unique(d_s) if len(d_s) else np.zeros(0, np.int64)
    smap = {int(c): i for i, c in enumerate(src_cells)}
    tmap = {int(c): i for i, c in enumerate(dests)}
    ns, ng, nt = len(src_cells), len(gsig), len(dests)
    S = 0
    off_src = 1
    off_grp = off_src + ns
    off_dst = off_grp + ng
    T = off_dst + nt
    n_nodes = T + 1
    rows, cols, caps = [], [], []
    for c in src_cells:
        cap = int(loss_cap[c]) if c >= 0 else int(gcnt.sum())  # extracellular: uncapped
        rows.append(S)
        cols.append(off_src + smap[int(c)])
        caps.append(cap)
    for gi in range(ng):
        s_ = smap[int(gsig[gi, 0])]
        n_i = int(gcnt[gi])
        rows.append(off_src + s_)
        cols.append(off_grp + gi)
        caps.append(n_i)
        for d_ in gsig[gi, 1:]:
            if d_ < 0:
                break
            rows.append(off_grp + gi)
            cols.append(off_dst + tmap[int(d_)])
            caps.append(n_i)
    for c in dests:
        rows.append(off_dst + tmap[int(c)])
        cols.append(T)
        caps.append(int(gain_cap[c]))
    graph = csr_matrix((np.array(caps, np.int64), (np.array(rows), np.array(cols))),
                       shape=(n_nodes, n_nodes))
    res = maximum_flow(graph, S, T)
    nu = int(res.flow_value)
    n_arcs = len(caps)

    # Witness: an integral flow must be realisable as an assignment in U, so read the
    # group -> cell flows and check both budgets and the total. A value with no witness
    # would mean the network admits something U does not.
    flow = res.flow.tocsr()
    w_ok = True
    if ng:
        sub = flow[off_grp:off_grp + ng, off_dst:off_dst + nt]
        out_group = np.asarray(sub.sum(axis=1)).ravel()          # per group, units routed
        inflow = np.asarray(sub.sum(axis=0)).ravel()             # per destination cell
        if (out_group > gcnt).any():
            w_ok = False
        if (inflow > gain_cap[dests]).any():
            w_ok = False
        out_src = np.zeros(ns, dtype=np.int64)
        for gi in range(ng):
            out_src[smap[int(gsig[gi, 0])]] += int(out_group[gi])
        for si, c in enumerate(src_cells):
            if c >= 0 and out_src[si] > loss_cap[c]:
                w_ok = False
        if int(out_group.sum()) != nu:
            w_ok = False

    # capped count: sum over counted cells of min(donor slots, gain cap), clamped by
    # the per-transcript count; pot counts every donor slot of this gene, including
    # transcripts already in C.  On breast the cap never binds, so the clamp decides
    # and the column equals the per-transcript count on all 150 gene-ends.
    pot_all = np.bincount(dcell_all, minlength=n_cells) if len(dcell_all) else \
        np.zeros(n_cells, np.int64)
    dmax_sum = int(np.minimum(pot_all[target_mask], gain_cap[target_mask]).sum())
    dmax_bin = int(len(np.unique(tid)))
    dmax_sum = int(min(dmax_bin, dmax_sum))
    return nu, dmax_sum, dmax_bin, n_nodes, n_arcs, bool(w_ok)


def main():
    t0 = time.time()
    tx, dm, labels, regions, base_tc = load()
    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()
    n_cells = len(base_tc)
    loss_cap = np.floor(LOSS_FRAC * base_tc).astype(np.int64)
    gain_cap = np.floor(GAIN_FRAC * base_tc).astype(np.int64)

    band = dm["band_tested_indices"]
    donor_sets = dm["donor_sets"]
    lens = np.array([len(d) for d in donor_sets])
    donor_flat = np.concatenate([np.asarray(d, np.int64) for d in donor_sets if len(d)]) \
        if lens.sum() else np.zeros(0, np.int64)
    donor_ptr = np.r_[0, np.cumsum(lens)]
    mov_gene = feat[band]
    mov_cell = cell[band]
    print(f"movable transcripts: {len(band):,}; donor slots {len(donor_flat):,}", flush=True)

    rows = []
    cache = {}
    for pair in ALL_PAIRS:
        Lg, Rg = pair
        for region_id in range(3):
            for gene, tmap_ in ((Lg, TYPE_MAP), (Rg, RECV_MAP)):
                target = (regions == region_id) & (labels == tmap_[gene])
                if target.sum() == 0:
                    continue
                key = (gene, region_id, tmap_[gene])
                if key not in cache:
                    m = np.flatnonzero(mov_gene == gene)
                    t1 = time.time()
                    nu, dsum, dbin, nn, na, wok = exact_positive(
                        m, mov_cell, donor_ptr, donor_flat, target,
                        loss_cap, gain_cap, n_cells)
                    free_loss = np.full(n_cells, int(dbin) + 1, dtype=np.int64)
                    nu_nl, _, _, _, _, wok_nl = exact_positive(
                        m, mov_cell, donor_ptr, donor_flat, target,
                        free_loss, gain_cap, n_cells)
                    cache[key] = (nu, dsum, dbin, nn, na, time.time() - t1, nu_nl,
                                  bool(wok and wok_nl))
                nu, dsum, dbin, nn, na, secs, nu_nl, wok = cache[key]
                cells_m = mov_cell[np.flatnonzero(mov_gene == gene)]
                src_in = np.where(cells_m >= 0, target[np.maximum(cells_m, 0)], False)
                a_s = np.bincount(cells_m[src_in], minlength=n_cells)
                dmin_count = -int(src_in.sum())
                dmin_exact = -int(np.minimum(a_s[target], loss_cap[target]).sum())
                rows.append(dict(pair=f"{Lg}->{Rg}", region=region_id, gene=gene,
                                 side="pos", n_nodes=nn, n_arcs=na,
                                 d_count=dbin, d_sum=dsum, d_exact=nu, d_noloss=nu_nl,
                                 gap_count=dbin - nu, gap_sum=dsum - nu,
                                 witness=bool(wok), secs=round(secs, 3)))
                rows.append(dict(pair=f"{Lg}->{Rg}", region=region_id, gene=gene,
                                 side="neg", n_nodes=nn, n_arcs=na,
                                 d_count=dmin_count, d_sum=dmin_count, d_exact=dmin_exact,
                                 gap_count=dmin_count - dmin_exact,
                                 gap_sum=dmin_count - dmin_exact,
                                 secs=0.0))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "exact_end.csv"), index=False)

    # per-row interval comparison
    sealed = pd.read_csv(os.path.join(OUT, "scout_bounds.csv"))
    comp = []
    for _, r in sealed.iterrows():
        Lg, Rg = r.pair.split("->")
        sub = df[(df.pair == r.pair) & (df.region == r.region)]
        dLp = sub[(sub.gene == Lg) & (sub.side == "pos")].iloc[0]
        dLn = sub[(sub.gene == Lg) & (sub.side == "neg")].iloc[0]
        dRp = sub[(sub.gene == Rg) & (sub.side == "pos")].iloc[0]
        dRn = sub[(sub.gene == Rg) & (sub.side == "neg")].iloc[0]
        NL0, NR0 = int(r.nL), int(r.nR)
        nls, nrs = int(r.nl_cells), int(r.nr_cells)
        T0 = t_log2(NL0, NR0, nls, nrs)
        lo_s, hi_s = r.T_lo, r.T_hi
        lo_c = t_log2(max(0, NL0 + int(dLn.d_count)), max(0, NR0 + int(dRn.d_count)), nls, nrs)
        hi_c = t_log2(NL0 + int(dLp.d_count), NR0 + int(dRp.d_count), nls, nrs)
        lo_e = t_log2(max(0, NL0 + int(dLn.d_exact)), max(0, NR0 + int(dRn.d_exact)), nls, nrs)
        hi_e = t_log2(NL0 + int(dLp.d_exact), NR0 + int(dRp.d_exact), nls, nrs)
        comp.append(dict(pair=r.pair, region=int(r.region), T0=T0,
                         T_lo_sealed=lo_s, T_hi_sealed=hi_s,
                         T_lo_count=lo_c, T_hi_count=hi_c,
                         T_lo_exact=lo_e, T_hi_exact=hi_e,
                         dT_lo=lo_e - lo_s, dT_hi=hi_e - hi_s,
                         cap_dT_hi=hi_c - hi_e, cap_dT_lo=lo_c - lo_e,
                         cert_pos_sealed=bool(lo_s > 0), cert_neg_sealed=bool(hi_s < 0),
                         cert_pos_count=bool(lo_c > 0), cert_neg_count=bool(hi_c < 0),
                         cert_pos_exact=bool(lo_e > 0), cert_neg_exact=bool(hi_e < 0)))
    cdf = pd.DataFrame(comp)
    cdf.to_csv(os.path.join(OUT, "exact_end_rows.csv"), index=False)

    pos = df[df.side == "pos"]
    neg = df[df.side == "neg"]
    summ = {
        "n_gene_combos": int(len(pos)),
        "pos_exact_equal_count": int((pos.gap_count == 0).sum()),
        "pos_exact_equal_noloss": int((pos.d_count == pos.d_noloss).sum()),
        "pos_witness_ok": int(pos.witness.sum()),
        "pos_gap_capped_only": int(((pos.gap_count > 0) & (pos.d_noloss == pos.d_count)).sum()),
        "pos_exact_equal_sum": int((pos.gap_sum == 0).sum()),
        "pos_max_gap_count": int(pos.gap_count.max()),
        "pos_max_gap_sum": int(pos.gap_sum.max()),
        "pos_mean_gap_sum": float(pos.gap_sum.mean()),
        "neg_exact_equal_count": int((neg.gap_count == 0).sum()),
        "neg_max_gap_count": int(neg.gap_count.max()),
        "nodes_max": int(pos.n_nodes.max()),
        "nodes_median": float(pos.n_nodes.median()),
        "arcs_max": int(pos.n_arcs.max()),
        "arcs_median": float(pos.n_arcs.median()),
        "secs_total": float(pos.secs.sum()),
        "n_rows": int(len(cdf)),
        "rows_hi_changed": int((cdf.dT_hi.abs() > 1e-12).sum()),
        "rows_lo_changed": int((cdf.dT_lo.abs() > 1e-12).sum()),
        "max_dT_hi": float(cdf.dT_hi.abs().max()),
        "max_dT_lo": float(cdf.dT_lo.abs().max()),
        "sign_changes": int(((cdf.cert_pos_sealed != cdf.cert_pos_exact) |
                             (cdf.cert_neg_sealed != cdf.cert_neg_exact)).sum()),
        "cert_pos_exact": int(cdf.cert_pos_exact.sum()),
        "cert_neg_exact": int(cdf.cert_neg_exact.sum()),
        "cert_pos_sealed": int(cdf.cert_pos_sealed.sum()),
        "cert_neg_sealed": int(cdf.cert_neg_sealed.sum()),
        "cert_pos_count": int(cdf.cert_pos_count.sum()),
        "cert_neg_count": int(cdf.cert_neg_count.sum()),
        "cap_rows_moved": int(((cdf.cap_dT_hi.abs() > 1e-9) |
                              (cdf.cap_dT_lo.abs() > 1e-9)).sum()),
        "cap_max_dT": float(max(cdf.cap_dT_hi.abs().max(), cdf.cap_dT_lo.abs().max())),
        "cap_sign_changes": int(((cdf.cert_pos_count != cdf.cert_pos_exact) |
                                 (cdf.cert_neg_count != cdf.cert_neg_exact)).sum()),
    }
    with open(os.path.join(OUT, "exact_end_summary.json"), "w") as f:
        json.dump(summ, f, indent=1, sort_keys=True)
    print(json.dumps(summ, indent=1, sort_keys=True))
    print(f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
