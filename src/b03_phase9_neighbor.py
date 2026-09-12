"""B03 Phase 9: neighbor-effect statistic — observed side + registered gates.

Prereg: B03_PHASE9_PREREG.md (FROZEN before any compute) + the 2026-09-11
addendum (frozen before corrected compute; see amendment record in B03_DEVIATIONS.md).
amendment record amends G9-2's wording (centroid-tree vs polygon brute-force are
different geometries; the registered check is tree-vs-brute-force-tree
agreement on the SAME metric).

Design (exactness-preserving):
- Adjacency eligibility is computed ONCE from frozen vendor centroids and frozen
  scout labels (labels.npy / regions.npy, sealed Phase-1 artifacts). It is a
  function of geometry and labels only, NOT of the assignment uncertainty set,
  so denominators (nl_adj, nr_adj) are fixed numbers and every bound operation
  remains a restriction of the same bipartite graph the exactness proof covers.
- The counted-set restriction is implemented by passing the eligible cell-id
  set DIRECTLY into the corrected interval builder (addendum section 1):
  counted(sender) = {(region==rid) & (labels==S) & send_m}, counted(receiver)
  analogous. Donor sets (band geometry, kernel intervals) are untouched, so a
  band transcript on an ineligible sender legitimately contributes a +1 delta
  (it can reassign INTO an eligible cell under the uncertainty set) and a
  transcript on an eligible sender contributes its usual -1/0/+1 range.
- Interval construction: corrected per-gene-slot rule (amendment record) with the
  amendment record capacity formula (the sealed null's own rule, one formula everywhere):
    dmin(g) = -(# band transcripts of g with source in the counted set)  [exact]
    dmax(g) = min( #{out-of-set band transcripts of g with >=1 own donor slot in
                     the counted set},  sum_{c in counted} min(pot_c(g), cap_c) )
  On breast the capacity never binds (sealed capacity_check.json), so dmax is
  the attained per-gene optimum; the formula is retained for the one-rule
  discipline and is verified non-binding by assertion.
- G9-0 (artifact integrity): recomputed scout labels must equal labels.npy
  exactly, proving the frozen artifact is the sealed scout's.
- G9-1 (amended per addendum section 2): with the all-eligible counted set,
  the corrected builder must reproduce the amendment record corrected bounds
  (results/sensitivity/breast_tight_bounds_sensitivity.csv) BIT-EXACTLY on
  every (dL_min, dL_max, dR_min, dR_max), AND the sealed scout_bounds.csv
  intervals must contain the corrected intervals row-wise (superset check —
  the soundness direction).
- G9-2 (adjacency determinism, amendment record): eligibility masks recomputed by
  brute-force per-cell radius queries must equal the KD-tree masks exactly
  (explicit self-exclusion for same-type pairs).

This script computes the OBSERVED side only. p-values require the matched
worst-case permutation null (b03_phase9_null.py) under the SAME corrected
construction, registered at execution time (single-null-per-machine rule).
"""
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from b03_certify import (ALL_PAIRS, POSITIVE_CONTROLS, T_bounds_from_count_intervals,
                         assign_cell_types, load_all, make_regions)

DATA = "B03_project/data"
OUT = "B03_project/results"
SENS = f"{OUT}/sensitivity/breast_tight_bounds_sensitivity.csv"

# Frozen scout maps (sealed b03_certify main(); module-level so main() and G9-1
# share one source. Copied verbatim from the sealed file.)
TYPE_MAP = {
    "CD274": "Breast cancer", "CXCL12": "Fibroblasts", "CCL5": "Macrophages",
    "PECAM1": "Endothelial cells", "CD163": "Macrophages", "ERBB2": "Breast cancer",
    "ESR1": "Breast cancer", "PGR": "Breast cancer", "ACTA2": "Smooth muscle cells",
    "KRT8": "Breast glandular cells", "CDH1": "Breast glandular cells",
    "CD3D": "T cells", "MS4A1": "B cells", "CD68": "Macrophages",
}
RECV_MAP = {
    "PDCD1": "T cells", "CXCR4": "T cells", "CCR7": "T cells", "KDR": "Endothelial cells",
    "CD3D": "T cells", "EGFR": "Breast cancer", "PGR": "Breast cancer",
    "ESR1": "Breast cancer", "CTLA4": "T cells", "CD8A": "T cells",
    "ERBB2": "Breast cancer", "CD274": "Breast cancer",
}

D_RADIUS_UM = float(os.environ.get("B03_P9_RADIUS", "15.0"))  # registered primary
# Sensitivity radius (registered: 25 um) runs as a second pass with B03_P9_RADIUS=25.0;
# outputs carry the radius column so downstream joining is unambiguous.


def g90_label_artifact_integrity(tx, panel, labels_saved):
    tested = sorted({g for p in ALL_PAIRS for g in p})
    labels_recomputed = assign_cell_types(tx, panel, set(tested), len(labels_saved))
    if not np.array_equal(np.asarray(labels_recomputed, dtype=object),
                          np.asarray(labels_saved, dtype=object)):
        raise SystemExit("G9-0 FAIL: labels.npy is not the sealed scout artifact")
    print("G9-0 PASS: labels.npy == recomputed scout labels (exact)", flush=True)


def eligible_masks(centroids, labels, s_type, r_type, radius):
    """Sender-eligible (type S with >=1 type-R neighbor) and receiver-eligible."""
    tree = cKDTree(centroids)
    pairs = tree.query_pairs(radius, output_type="ndarray")  # i<j by index order
    is_s = labels == s_type
    is_r = labels == r_type
    send_near_r = np.zeros(len(labels), dtype=bool)
    recv_near_s = np.zeros(len(labels), dtype=bool)
    i, j = pairs[:, 0], pairs[:, 1]
    # direction i->j
    send_near_r[i[is_r[j]]] = True
    recv_near_s[j[is_s[i]]] = True
    # direction j->i
    send_near_r[j[is_r[i]]] = True
    recv_near_s[i[is_s[j]]] = True
    return send_near_r, recv_near_s


def g02_adjacency_determinism(centroids, labels, s_type, r_type, radius,
                              send_mask, recv_mask, sample_n=2000, seed=9):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(labels), size=min(sample_n, len(labels)), replace=False)
    tree = cKDTree(centroids)
    is_s = labels == s_type
    is_r = labels == r_type
    for name, mask, own_type, other_type in (
            ("senders", send_mask, is_s, is_r), ("receivers", recv_mask, is_r, is_s)):
        sample = idx[own_type[idx]]
        if len(sample) == 0:
            continue
        nbrs_all = tree.query_ball_point(centroids[sample], r=radius)
        brute = np.zeros(len(sample), dtype=bool)
        for k_i, cell_i in enumerate(sample):
            # Self is NOT a neighbor (registered semantics: adjacency is between
            # distinct cells; query_pairs in the tree mask excludes self). For
            # same-type pairs (e.g. PECAM1->KDR, both Endothelial) the ball
            # contains the cell itself, which must be excluded explicitly. The
            # neighbor test uses OTHER_TYPE (senders need a receiver-type
            # neighbor; receivers a sender-type) — the original draft hardcoded
            # is_r here, which failed G9-2 for receivers at first execution
            # (logged amendment record).
            if any(other_type[nbr] and nbr != cell_i for nbr in nbrs_all[k_i]):
                brute[k_i] = True
        if not np.array_equal(brute, mask[sample]):
            raise SystemExit(f"G9-2 FAIL: {name} mask mismatch (tree vs brute force)")
    print(f"G9-2 PASS: {s_type}->{r_type} masks reproduce brute force on sample",
          flush=True)


class CorrectedBounds:
    """amendment record corrected per-gene-slot interval builder (one rule everywhere).

    dmin(g) = -(# band transcripts of g whose source cell is in the counted set) —
    exact (extracellular destination always available).
    dmax(g) = min( binary-donor count, amendment record capacity sum ) where the binary
    count is over each transcript's OWN donor slots only (no cross-gene slot
    absorption), and the amendment record capacity sum is sum_{c in counted}
    min(pot_c(g), cap_c) with pot_c(g) fixed per gene and cap_c the frozen
    gain capacity. On breast caps never bind (asserted).
    """

    def __init__(self, tx, dm):
        feat = tx.feature_name.to_numpy()
        cell = tx.cell_idx.to_numpy()
        self.band_tested = np.asarray(dm["band_tested_indices"])
        self.donor_sets = dm["donor_sets"]
        self.mov_cell = cell[self.band_tested]
        self.mov_gene_name = feat[self.band_tested]
        lens = np.array([len(d) for d in self.donor_sets])
        self.donor_flat = np.concatenate(
            [np.asarray(d, dtype=np.int64) for d in self.donor_sets if len(d)])
        self.donor_ptr = np.r_[0, np.cumsum(lens)]
        self.gidx = {g: i for i, g in enumerate(sorted(set(self.mov_gene_name.tolist())))}
        self.mov_gid = np.array([self.gidx[g] for g in self.mov_gene_name], dtype=np.int64)
        cells_meta = pd.read_parquet(f"{DATA}/cells_meta.parquet",
                                     columns=["transcript_counts"])
        self.cap = np.floor(0.5 * cells_meta.transcript_counts.to_numpy(np.int64))
        self.n_cells = len(self.cap)
        # per-gene own-slot index (concatenated ranges ps[i]:pe[i] for this gene)
        self.own = {}
        for g in self.gidx:
            m = np.flatnonzero(self.mov_gid == self.gidx[g])
            ps, pe = self.donor_ptr[m], self.donor_ptr[m + 1]
            ok = pe > ps
            if ok.any():
                idx = np.concatenate([np.arange(a, b) for a, b in zip(ps[ok], pe[ok])])
                offs = np.r_[0, np.cumsum((pe - ps)[ok])][:-1]
            else:
                idx = np.zeros(0, np.int64)
                offs = np.zeros(0, np.int64)
            slots = self.donor_flat[idx] if len(idx) else np.zeros(0, np.int64)
            pot = np.bincount(slots, minlength=self.n_cells).astype(np.int64) \
                if len(slots) else np.zeros(self.n_cells, np.int64)
            self.own[g] = (m, self.mov_cell[m], ok, idx, offs, pot)

    def interval(self, gene, counted_ids):
        """(dmin, dmax) for gene over counted cell-id set; amendment record cap applied."""
        m, cells_m, ok, idx, offs, pot = self.own[gene]
        counted_ids = np.asarray(counted_ids, dtype=np.int64)
        src_in = np.isin(cells_m, counted_ids) & (cells_m >= 0)
        dmin = -int(src_in.sum())
        if len(idx):
            anyin = np.zeros(len(m), bool)
            anyin[ok] = np.logical_or.reduceat(
                np.isin(self.donor_flat[idx], counted_ids), offs)
            dmax_bin = int(((~src_in) & anyin).sum())
        else:
            dmax_bin = 0
        in_c = np.zeros(self.n_cells, bool)
        if pot.size:
            in_c[counted_ids[counted_ids < self.n_cells]] = True
            dmax_cap = int(np.minimum(pot[in_c], self.cap[in_c]).sum())
        else:
            dmax_cap = 0
        dmax = min(dmax_bin, dmax_cap)
        assert dmax == dmax_bin, "amendment record capacity bound on breast (must not bind)"
        return dmin, dmax


def g91_corrected_reproduction(cb, labels, region_of_cell, scout, sens):
    """Amended G9-1: corrected builder == amendment record corrected bounds bit-exactly;
    sealed intervals contain corrected intervals row-wise (superset)."""
    sens_idx = {(r.pair, int(r.region)): r for r in sens.itertuples()}
    n_ok = 0
    for _, row in scout.iterrows():
        Lg, Rg = row["pair"].split("->")
        rid = int(row["region"])
        s = sens_idx[(row["pair"], rid)]
        for G, tname, (dl_sealed, dl_corr_sens) in (
                (Lg, TYPE_MAP[Lg], ((int(row["dL_min"]), int(row["dL_max"])),
                                    (int(s.dL_min_tight) if hasattr(s, "dL_min_tight")
                                     else int(row["dL_min"]),
                                     int(s.dL_max_tight)))),
                (Rg, RECV_MAP[Rg], ((int(row["dR_min"]), int(row["dR_max"])),
                                    (int(row["dR_min"]), int(s.dR_max_tight))))):
            counted = np.flatnonzero((region_of_cell == rid) & (labels == tname))
            dmin, dmax = cb.interval(G, counted)
            assert (dmin, dmax) == (dl_corr_sens[0], dl_corr_sens[1]), \
                f"G9-1 FAIL corrected {row['pair']} r{rid} {G}: " \
                f"({dmin},{dmax}) vs {dl_corr_sens}"
            assert dl_sealed[0] <= dmin and dmax <= dl_sealed[1], \
                f"G9-1 FAIL superset {row['pair']} r{rid} {G}"
            n_ok += 1
    print(f"G9-1 PASS (amended): {n_ok} gene-side intervals reproduce amendment record "
          f"corrected bounds bit-exactly; sealed intervals contain them row-wise",
          flush=True)


def main():
    t0 = time.time()
    tx, dm, cells_meta, panel = load_all()
    labels = np.load(f"{DATA}/labels.npy", allow_pickle=True)
    region_of_cell = np.load(f"{DATA}/regions.npy")
    assert len(labels) == len(cells_meta) == len(region_of_cell)

    g90_label_artifact_integrity(tx, panel, labels)

    centroids = cells_meta[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    feat = tx.feature_name.to_numpy()
    cell = tx.cell_idx.to_numpy()
    scout = pd.read_csv(f"{OUT}/scout_bounds.csv")
    sens = pd.read_csv(SENS)

    cb = CorrectedBounds(tx, dm)
    g91_corrected_reproduction(cb, labels, region_of_cell, scout, sens)

    rows = []
    mask_cache = {}
    for (Lg, Rg) in ALL_PAIRS:
        s_type = TYPE_MAP.get(Lg)
        r_type = RECV_MAP.get(Rg)
        if s_type is None or r_type is None:
            continue
        key = (s_type, r_type)
        if key not in mask_cache:
            send_m, recv_m = eligible_masks(centroids, labels, s_type, r_type,
                                            D_RADIUS_UM)
            g02_adjacency_determinism(centroids, labels, s_type, r_type,
                                      D_RADIUS_UM, send_m, recv_m)
            mask_cache[key] = (send_m, recv_m)
        send_m, recv_m = mask_cache[key]
        for rid in range(3):
            nl_adj = int(((region_of_cell == rid) & (labels == s_type) & send_m).sum())
            nr_adj = int(((region_of_cell == rid) & (labels == r_type) & recv_m).sum())
            if nl_adj == 0 or nr_adj == 0:
                rows.append(dict(pair=f"{Lg}->{Rg}", region=rid, control=(Lg, Rg) in POSITIVE_CONTROLS,
                                 s_type=s_type, r_type=r_type, testable=False,
                                 nl_adj=nl_adj, nr_adj=nr_adj))
                continue
            # counted sets: eligible senders / receivers in this region (addendum:
            # counted set passed directly; identical to sentinel-masked labels)
            counted_S = np.flatnonzero((region_of_cell == rid) &
                                       (labels == s_type) & send_m)
            counted_R = np.flatnonzero((region_of_cell == rid) &
                                       (labels == r_type) & recv_m)
            dl = cb.interval(Lg, counted_S)
            dr = cb.interval(Rg, counted_R)
            # observed eligible counts
            cell_L = cell[feat == Lg]
            nL_adj = int(np.isin(cell_L, counted_S).sum())
            cell_R = cell[feat == Rg]
            nR_adj = int(np.isin(cell_R, counted_R).sum())
            T_lo, T_hi, _ = T_bounds_from_count_intervals(nL_adj, nR_adj,
                                                          dl, dr,
                                                          nl_adj, nr_adj)
            rows.append(dict(pair=f"{Lg}->{Rg}", region=rid,
                             control=(Lg, Rg) in POSITIVE_CONTROLS,
                             s_type=s_type, r_type=r_type, testable=True,
                             nL_adj=nL_adj, nR_adj=nR_adj,
                             nl_adj=nl_adj, nr_adj=nr_adj,
                             dL_min=dl[0], dL_max=dl[1],
                             dR_min=dr[0], dR_max=dr[1],
                             T_lo=T_lo, T_hi=T_hi,
                             certified_pos=bool(T_lo > 0), certified_neg=bool(T_hi < 0)))
    res = pd.DataFrame(rows)
    res["radius_um"] = D_RADIUS_UM
    os.makedirs(OUT, exist_ok=True)
    suffix = "" if D_RADIUS_UM == 15.0 else f"_d{int(D_RADIUS_UM)}"
    res.to_csv(f"{OUT}/phase9_neighbor_observed{suffix}.csv", index=False)
    n_t = int(res["testable"].sum())
    print(f"rows: {len(res)} ({n_t} testable) | "
          f"cert_pos {int(res.get('certified_pos', pd.Series(dtype=bool)).sum())} "
          f"cert_neg {int(res.get('certified_neg', pd.Series(dtype=bool)).sum())}",
          flush=True)
    print(f"[{time.time()-t0:.0f}s] Phase 9 observed side complete (no inference; "
          f"null is a separate registered run)", flush=True)


if __name__ == "__main__":
    main()
