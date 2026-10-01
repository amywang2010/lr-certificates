"""B03 semi-synthetic certification study runner (design frozen before compute).

Builds pseudo-sections on the real breast Rep1 molecule file. The ONLY things
interventions touch are feature-name assignments of transcript rows; molecule
coordinates, cell_idx, dist_boundary, in_band, per-cell totals, regions, band
membership, donor slots, and capacities are identical to the shipped section
in every run (donor geometry is name- and label-independent).

Representation: transcript names as integer codes (codes) over a uniq name
table. Interventions:
  (1) Arm N: seeded permutation of the 21 tested names among themselves
      (code-remap; labels, regions, geometry, tested mass all unchanged).
  (2) Arm M: seeded transcript-level REMOVAL of tested transcripts (a smaller
      panel does not record transcripts; removal is implemented by masking
      codes to a sentinel outside every panel/marker/tested list), protecting
      tested transcripts inside anchor cells (the injection substrate), toward
      an aggregate tested mass-share target.
  (3) Arms E/M: known-positive injection for (KRT8 -> ERBB2) by borrowing
      same-class tested transcripts inside the anchor cells of the pair's
      sender/receiver classes at fraction f.

The tested-name closure invariant holds for every run: the pseudo-section's
tested-name transcript population is a subset of the original one, so the
stored donor geometry (built for originally band-and-tested transcripts),
rebuilt positionally, stays exact.

Each run executes the shipped interval engine logic (b03_scout.py conventions,
in-process) and the shipped matched-null machinery (b03_robust_null_v3.py
functions, in-process; identical per-permutation seeds, blocks, tie convention,
BH, output schema), writing bounds/null/final tables plus a run manifest into
results/synthstudy/<run_id>/.

Usage (paths are repository-relative; run from anywhere):
  python src/b03_synthstudy.py --mode verify
  python src/b03_synthstudy.py --mode smoke
  python src/b03_synthstudy.py --mode run --arm N --run-id N01 --seed 910001
  python src/b03_synthstudy.py --mode run --arm E --run-id E101 --seed 911001 --f 1.0
  python src/b03_synthstudy.py --mode run --arm M --run-id M0801 --target 8.0 --seed 912001
  python src/b03_synthstudy.py --mode fleet
  python src/b03_synthstudy.py --mode analyze
"""
import argparse
import itertools
import json
import math
import multiprocessing as mp
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "results"
STUDY = OUT / "synthstudy"
sys.path.insert(0, str(SRC))  # spawned workers must re-import the null module by name

SEED = 20260905  # identical to the shipped default: permutation stream matches
EPS = 0.5
MASKED = "__REMOVED__"  # sentinel: matches no panel, marker, or tested name

ALL_PAIRS = [("CD274", "PDCD1"), ("CXCL12", "CXCR4"), ("CCL5", "CCR7"),
             ("PECAM1", "KDR"), ("CD163", "CD3D"),
             ("ERBB2", "EGFR"), ("ERBB2", "PDCD1"), ("ESR1", "PGR"), ("PGR", "ESR1"),
             ("CD274", "CTLA4"), ("CD274", "CD8A"), ("CXCL12", "CCR7"), ("CCL5", "CXCR4"),
             ("PECAM1", "PDCD1"), ("PECAM1", "CTLA4"), ("CD163", "CCR7"), ("CD163", "CXCR4"),
             ("ACTA2", "EGFR"), ("ACTA2", "ERBB2"), ("KRT8", "ERBB2"), ("CDH1", "EGFR"),
             ("CDH1", "ERBB2"), ("CD3D", "ERBB2"), ("MS4A1", "CD274"), ("CD68", "CD274")]
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
TESTED = sorted({g for p in ALL_PAIRS for g in p})

INJ_PAIR = ("KRT8", "ERBB2")
INJ_CLASSES = (TYPE_MAP[INJ_PAIR[0]], RECV_MAP[INJ_PAIR[1]])
M_TARGETS = (8.0, 4.0, 1.0, 0.5)
N_RUNS_TOTAL = 20

SYS = None  # base data cache (set by load_base())


def log(msg, run_id="study"):
    line = f"[{time.strftime('%H:%M:%S')}] [{run_id}] {msg}"
    print(line, flush=True)
    STUDY.mkdir(parents=True, exist_ok=True)
    with open(STUDY / "fleet.log", "a", encoding="utf-8") as f:
        f.write(line + "\n")


# ---------------------------------------------------------------------------
# base data (donor-map load dominates: ~30-60 s)
# ---------------------------------------------------------------------------

def load_base():
    global SYS
    if SYS is not None:
        return SYS
    t0 = time.time()
    tx = pd.read_parquet(DATA / "tx.parquet")
    cm0 = pd.read_parquet(DATA / "cells_meta.parquet")
    panel = pd.read_csv(DATA / "Xenium_FFPE_Human_Breast_Cancer_Rep1_panel.tsv", sep="\t")
    log(f"donor map loading (t={time.time()-t0:.0f}s at start)", "study")
    with open(DATA / "donor_map.pkl", "rb") as f:
        dm = pickle.load(f)
    # polygon idx -> cells_meta row (shipped scout alignment)
    cid_map = {str(k): v for k, v in dm["cid_map"].items()}
    cm_ids = cm0.cell_id.astype(str).to_numpy()
    poly_row = np.full(max(cid_map.values()) + 1, -1, dtype=np.int64)
    for row, cid in enumerate(cm_ids):
        p = cid_map.get(cid, -1)
        if p >= 0:
            poly_row[p] = row
    assert (poly_row >= 0).all(), "unmapped polygons present"
    cells_aligned = cm0.iloc[poly_row].reset_index(drop=True)  # row == polygon idx
    SYS = dict(tx=tx, panel=panel, dm=dm, cells_aligned=cells_aligned,
               labels=np.load(DATA / "labels.npy", allow_pickle=True))
    log(f"base data loaded in {time.time()-t0:.0f}s", "study")
    return SYS


# ---------------------------------------------------------------------------
# shared geometry helpers
# ---------------------------------------------------------------------------

def proximity_mask(anchor_seed, k=300, near=100):
    """Anchor set: k seeded cells plus their `near` nearest neighbours per
    anchor (centroid space, squared Euclidean). Deterministic in seed."""
    cm = SYS["cells_aligned"]
    xy = np.c_[cm.x_centroid.to_numpy(float), cm.y_centroid.to_numpy(float)]
    rng = np.random.default_rng(anchor_seed)
    centers = xy[rng.choice(len(xy), size=k, replace=False)]
    d2 = ((xy[:, None, :] - centers[None, :, :]) ** 2).sum(-1)
    idx = np.argpartition(d2, near, axis=0)[:near]          # (near, k)
    mask = np.zeros(len(xy), dtype=bool)
    mask[idx.ravel()] = True
    return mask


def tx_in_anchor(mask):
    cell = SYS["tx"].cell_idx.to_numpy()
    in_a = np.zeros(len(cell), dtype=bool)
    ap = cell >= 0
    in_a[ap] = mask[np.clip(cell[ap], 0, len(mask) - 1)]
    return in_a


def tested_class_pool(cls, exclude):
    """Tested genes of a given annotation class (borrowing sources)."""
    pool = [g for g, c in list(TYPE_MAP.items()) + list(RECV_MAP.items())
            if c == cls and g not in exclude]
    return sorted(set(pool))


def tx_class_int(labels):
    """Per-transcript sender/receiver class id (int; -1 unassigned/out of cell),
    restricted to the two injection classes (everything else -> -1)."""
    cell = SYS["tx"].cell_idx.to_numpy()
    cls_ids = {c: i for i, c in enumerate(INJ_CLASSES)}
    lab_int = np.array([cls_ids.get(l, -1) for l in labels], dtype=np.int8)
    out = np.full(len(cell), -1, dtype=np.int8)
    ok = cell >= 0
    out[ok] = lab_int[np.clip(cell[ok], 0, len(lab_int) - 1)]
    return out


# ---------------------------------------------------------------------------
# arms (spec S2)
# ---------------------------------------------------------------------------

def arm_null(seed):
    rng = np.random.default_rng(seed)
    return dict(pi_seed=int(seed), injection=None,
                spec=dict(arm="N", seed=int(seed),
                          note="pi = random permutation of the 21 tested names"))


def arm_effect(anchor_seed, f):
    """Arm E, amended (spec S7.7): base-null permutation FIRST (tested names
    scrambled, every pair null by construction), THEN the known-positive
    injection. Recovery of a certified positive on KRT8->ERBB2 therefore
    measures detection of an injected effect, not recovery of a pair that
    already certifies on the real section."""
    rng = np.random.default_rng(anchor_seed * 10 + int(f * 10))
    return dict(pi_seed=int(anchor_seed) + 100000,
                injection=dict(anchor_seed=int(anchor_seed), f=float(f),
                               pair=list(INJ_PAIR), classes=list(INJ_CLASSES),
                               rng_seed=int(rng.integers(0, 2**31 - 1))),
                spec=dict(arm="E", seed=int(anchor_seed), f=float(f),
                          note="base-null pi (seed anchor+100000) then injection"))


def arm_mass(target_pct, seed):
    rng = np.random.default_rng(seed)
    return dict(pi_seed=None,
                injection=dict(anchor_seed=911001, f=1.0,
                               pair=list(INJ_PAIR), classes=list(INJ_CLASSES),
                               rng_seed=int(rng.integers(0, 2**31 - 1))),
                spec=dict(arm="M", seed=int(seed), target_pct=float(target_pct),
                          note="seeded transcript-level removal to target share "
                               "(anchor substrate protected), then injection "
                               "at f=1.0 of surviving borrowable mass"))


# ---------------------------------------------------------------------------
# pseudo-section construction (codes + uniq representation)
# ---------------------------------------------------------------------------

def build_pseudo_section(run):
    """Returns (codes2, uniq2, info). codes2: int64 per transcript into uniq2
    (object name array; MASKED appended when removal ran)."""
    tx = SYS["tx"]
    codes, uniq = pd.factorize(tx.feature_name, sort=False)
    codes = np.asarray(codes, dtype=np.int64)
    uniq = np.asarray(uniq, dtype=object)
    info = {}
    name_code = {g: i for i, g in enumerate(uniq.tolist())}
    missing = [g for g in TESTED if g not in name_code]
    assert not missing, f"tested genes absent from tx: {missing}"
    tested_codes = np.array([name_code[g] for g in TESTED], dtype=np.int64)

    # (1) Arm N: permute tested names among themselves via code remap
    if run["pi_seed"] is not None:
        rng = np.random.default_rng(run["pi_seed"])
        perm = rng.permutation(len(tested_codes))
        remap = np.arange(len(uniq), dtype=np.int64)
        remap[tested_codes] = tested_codes[perm]
        codes = remap[codes]
        info["pi_applied"] = True

    # (2) Arm M: transcript-level removal toward target share
    removed = None
    if run["spec"]["arm"] == "M":
        codes, uniq, removed = remove_tested(codes, uniq, tested_codes,
                                             run["spec"]["target_pct"],
                                             run["spec"]["seed"])
        info["removed"] = removed

    # (3) injection (Arms E and M)
    if run["injection"] is not None:
        codes, uniq = inject(codes, uniq, name_code, run["injection"])
        info["injection"] = {k: v for k, v in run["injection"].items()
                             if k != "pair" and k != "classes"}
    return codes, uniq, info


def remove_tested(codes, uniq, tested_codes, target_pct, seed):
    """Uniform transcript-level thinning (spec S7.3): a seeded uniform subset of
    ALL tested transcripts is removed toward the aggregate share target, with
    no anchor protection. The certified-bound construction is invariant to
    which molecules are removed (interior molecules do not move under U), and
    the anchor injection is re-applied post-thinning at f = 1.0 of the
    surviving borrowable mass, so the injected contrast persists in the
    anchor cells. Removal = masking codes to a sentinel."""
    rng = np.random.default_rng(seed)
    uniq2 = np.append(uniq, MASKED)
    masked_code = len(uniq2) - 1
    n_total = len(codes)
    in_tested = np.isin(codes, tested_codes)
    tested_mass = int(in_tested.sum())
    to_remove = max(0, tested_mass - int(round(target_pct / 100.0 * n_total)))
    removed = {}
    if to_remove > 0:
        cand = np.flatnonzero(in_tested)
        take = rng.choice(cand, size=to_remove, replace=False)
        names_removed, counts_removed = np.unique(uniq[codes[take]],
                                                  return_counts=True)
        removed = {str(g): int(c) for g, c in zip(names_removed, counts_removed)}
        codes[take] = masked_code
    return codes, uniq2, removed


def inject(codes, uniq, name_code, inj):
    """Known-positive: relabel same-class tested transcripts inside the anchor
    ring (nearest-neighbour ring around seeded anchors; ring radius scales
    with f so injected MASS is the dose axis — f <= 1 saturates the fixed
    anchor pool, so dose must grow the ring) to the pair names, taking the
    full borrowable mass of the ring's sender/receiver cells."""
    rng = np.random.default_rng(inj["rng_seed"])
    labels = np.load(DATA / "labels.npy", allow_pickle=True)
    cls_int = tx_class_int(labels)
    ring = max(10, int(round(100 * inj["f"])))
    anchor_tx = tx_in_anchor(proximity_mask(inj["anchor_seed"], near=ring))
    Lg, Rg = inj["pair"]
    Lcls, Rcls = inj["classes"]
    pool_L = [g for g in tested_class_pool(Lcls, exclude={Lg})]
    pool_R = [g for g in tested_class_pool(Rcls, exclude={Rg})]
    if not pool_L or not pool_R:
        raise RuntimeError("borrow pool empty")
    code_L, code_R = name_code[Lg], name_code[Rg]
    sel_L = anchor_tx & (cls_int == INJ_CLASSES.index(Lcls))
    sel_R = anchor_tx & (cls_int == INJ_CLASSES.index(Rcls))
    idx_L = np.flatnonzero(sel_L & np.isin(codes, [name_code[g] for g in pool_L]))
    idx_R = np.flatnonzero(sel_R & np.isin(codes, [name_code[g] for g in pool_R]))
    take_l = idx_L
    take_r = idx_R
    codes[take_l] = code_L
    codes[take_r] = code_R
    inj["swapped_L"] = int(len(take_l))
    inj["swapped_R"] = int(len(take_r))
    inj["ring_nn"] = int(ring)
    inj["borrow_pools"] = dict(L=pool_L, R=pool_R)
    return codes, uniq


# ---------------------------------------------------------------------------
# labels (shipped marker-vote rule on the pseudo-section; asserted identical)
# ---------------------------------------------------------------------------

def compute_labels(codes, uniq):
    tx = SYS["tx"]
    cell = tx.cell_idx.to_numpy()
    n_cells = len(SYS["cells_aligned"])
    ann = dict(zip(SYS["panel"].Name, SYS["panel"].Annotation))
    ann_by_code = np.array([ann.get(g, "NA") if g != MASKED else "NA"
                            for g in uniq.tolist()] + ["NA"], dtype=object)
    ann_per_tx = ann_by_code[codes]
    tested_set = set(TESTED)
    tested_by_code = np.array([g in tested_set if g != MASKED else False
                               for g in uniq.tolist()] + [False], dtype=bool)
    sel = (ann_per_tx != "NA") & (~tested_by_code[codes]) & (cell >= 0)
    df = pd.DataFrame({"c": cell[sel], "a": ann_per_tx[sel]})
    votes = df.groupby(["c", "a"]).size().reset_index(name="n")
    votes = votes.sort_values(["c", "n"], ascending=[True, False]).drop_duplicates("c")
    labels = np.full(n_cells, "Unknown", dtype=object)
    labels[votes.c.to_numpy()] = votes.a.to_numpy()
    ref = np.load(DATA / "labels.npy", allow_pickle=True)
    assert (labels == ref).all(), "labels diverged from shipped section"
    return labels


# ---------------------------------------------------------------------------
# interval engine (shipped scout conventions, in-process)
# ---------------------------------------------------------------------------

def t_log2(nL, nR, a, b):
    return math.log2(nL / a + EPS) + math.log2(nR / b + EPS)


def run_intervals(codes, uniq, run_dir: Path):
    tx = SYS["tx"]
    dm = SYS["dm"]
    cells_meta = SYS["cells_aligned"]
    n_cells = len(cells_meta)
    cell = tx.cell_idx.to_numpy()
    # spec S7.1c: labels are invariant across arms (tested-name closure: voting
    # reads only non-tested marker transcripts); computed once per fleet and
    # asserted against the shipped file in --mode verify
    labels = SYS["labels"]
    region_of_cell = np.load(DATA / "regions.npy")

    name_code = {g: i for i, g in enumerate(uniq.tolist())}
    gidx = {g: i for i, g in enumerate(TESTED)}
    tested_codes = np.array([name_code[g] for g in TESTED], dtype=np.int64)
    Gmat = np.zeros((n_cells, len(TESTED)), dtype=np.int64)
    for g in TESTED:
        sel_g = (codes == name_code[g]) & (cell >= 0)
        Gmat[:, gidx[g]] = np.bincount(cell[sel_g], minlength=n_cells)

    # band-tested geometry rebuilt positionally over the STORED band-tested set
    band_tested_full = dm["band_tested_indices"]
    donor_sets_full = dm["donor_sets"]
    keep = np.isin(codes[band_tested_full], tested_codes)
    band_tested = band_tested_full[keep]
    donor_sets = [donor_sets_full[i] for i in np.flatnonzero(keep)]
    lens = np.array([len(d) for d in donor_sets])
    donor_flat = (np.concatenate([np.asarray(d, dtype=np.int64) for d in donor_sets
                                  if len(d)]) if lens.sum() else np.array([], np.int64))
    donor_ptr = np.r_[0, np.cumsum(lens)]
    mov_gid = np.searchsorted(TESTED, uniq[codes[band_tested]])
    mov_cell = cell[band_tested]

    # capacity worst case (relabeling/removal changes per-cell totals never;
    # donor slots and band membership are untouched)
    base_tc = cells_meta.transcript_counts.to_numpy(np.int64)
    max_donors = np.bincount(donor_flat, minlength=n_cells) if len(donor_flat) \
        else np.zeros(1, np.int64)
    mov_in_cell = np.bincount(mov_cell[mov_cell >= 0], minlength=n_cells)
    cap_gain_ok = bool((max_donors <= 0.5 * base_tc).all())
    cap_loss_ok = bool((mov_in_cell <= 0.6 * base_tc).all())

    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    n_types = len(types_needed)
    SHIFT = 1
    type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
    n_slots = n_types + SHIFT
    cid = region_of_cell * n_slots + np.array(
        [type_shift.get(l, 0) for l in labels], dtype=np.int16)
    # ONE donor-class gather per run (perf: was ~300 gathers of 27.8M elements)
    cid_donor = cid[donor_flat] if len(donor_flat) else None

    def gene_interval(g, target):
        m = np.flatnonzero(mov_gid == gidx[g])
        if len(m) == 0:
            return 0, 0
        cells_m = mov_cell[m]
        src_in = np.where(cells_m >= 0, cid[np.maximum(cells_m, 0)], -999) == target
        ptr_s, ptr_e = donor_ptr[m], donor_ptr[m + 1]
        if cid_donor is not None:
            din = cid_donor == target
            offs = ptr_s - ptr_s[0]
            seg = din[ptr_s[0]:ptr_e[-1]]
            anyin = (np.logical_or.reduceat(seg, offs)
                     if len(offs) and len(seg) >= len(offs)
                     else np.zeros(len(m), bool))
        else:
            anyin = np.zeros(len(m), bool)
        return -int(src_in.sum()), int(((~src_in) & anyin).sum())

    controls = {("CD274", "PDCD1"), ("CXCL12", "CXCR4"), ("CCL5", "CCR7"),
                ("PECAM1", "KDR"), ("CD163", "CD3D")}
    rows = []
    for pair in ALL_PAIRS:
        Lg, Rg = pair
        S, R = TYPE_MAP[Lg], RECV_MAP[Rg]  # frozen type tables, all arms
        for region_id in range(3):
            cid_S = region_id * n_slots + type_shift[S]
            cid_R = region_id * n_slots + type_shift[R]
            mask_S, mask_R = cid == cid_S, cid == cid_R
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
            rows.append(dict(pair=f"{Lg}->{Rg}", control=pair in controls,
                             region=region_id, nL=N_L, nR=N_R,
                             nl_cells=nL_cells, nr_cells=nR_cells,
                             dL_min=dLm, dL_max=dLM, dR_min=dRm, dR_max=dRM,
                             T0=T0, T_lo=T_lo, T_hi=T_hi, width=T_hi - T_lo,
                             certified_pos=bool(T_lo > 0), certified_neg=bool(T_hi < 0)))
    res = pd.DataFrame(rows)
    res.to_csv(run_dir / "scout_bounds.csv", index=False)
    np.save(run_dir / "labels.npy", labels)  # per-run copy for the null state

    # randomized feasible-assignment verification (shipped convention:
    # sampled per-transcript deltas must land inside the certified interval)
    rng = np.random.default_rng(SEED)
    cid_donor = cid[donor_flat].astype(np.int32) if len(donor_flat) else None  # once
    ver_rows = []
    for _, r in res.iterrows():
        Lg, Rg = r.pair.split("->")
        region_id = int(r.region)
        for g, cls in [(Lg, TYPE_MAP[Lg]), (Rg, RECV_MAP[Rg])]:
            target = region_id * n_slots + type_shift[cls]
            m = np.flatnonzero(mov_gid == gidx[g])
            if len(m) == 0:
                continue
            cells_m = mov_cell[m]
            src_in = np.where(cells_m >= 0, cid[np.maximum(cells_m, 0)], -999) == target
            ptr_s, ptr_e = donor_ptr[m], donor_ptr[m + 1]
            if cid_donor is not None:
                din = cid_donor == target
                offs = ptr_s - ptr_s[0]
                seg = din[ptr_s[0]:ptr_e[-1]]
                anyin = (np.logical_or.reduceat(seg, offs)
                         if len(offs) and len(seg) >= len(offs)
                         else np.zeros(len(m), bool))
            else:
                anyin = np.zeros(len(m), bool)
            lo = np.where(src_in, -1, 0).astype(np.int8)
            hi = np.where((~src_in) & anyin, 1, 0).astype(np.int8)
            dmin_eff, dmax_eff = int(lo.sum()), int(hi.sum())
            # shipped convention, CHUNKED like b03_scout.py (CH=100): sample
            # per-transcript achievable deltas from the EXACT per-transcript
            # ranges and check totals against the interval; batching preserves
            # the per-config draws while removing per-config Python overhead.
            # Attainability evidence is the shipped subcube enumeration and the
            # independent re-derivation artifacts, not this loop.
            span = (hi - lo + 1).astype(np.int16)
            viol = 0
            CH = 100
            # 200 configs per gene-row (spec S7.9): this check catches range-
            # computation bugs; it is arithmetically incapable of failing for
            # a correct [lo, hi] (deltas drawn within per-transcript ranges
            # sum into the aggregate range), so 200 draws retain the full
            # bug-catch value at a fifth of the cost. Attainability evidence
            # remains the shipped subcube enumeration and independent
            # re-derivation artifacts.
            for _s in range(0, 200, CH):
                u = rng.random((CH, len(m)))
                d = (lo[None, :] + (u * span[None, :]).astype(np.int16)).astype(np.int16)
                tot = d.sum(axis=1)
                viol += int(((tot < dmin_eff) | (tot > dmax_eff)).sum())
            ver_rows.append(dict(pair=r.pair, region=region_id, gene=g, configs=200,
                                 violations=viol, dmin=dmin_eff, dmax=dmax_eff))
    ver = pd.DataFrame(ver_rows)
    ver.to_csv(run_dir / "exactness_verification.csv", index=False)

    in_tested = np.isin(codes, tested_codes)
    # Capacity semantics: the aggregate boolean is NOT a binding test (it counts
    # cells outside every counted set). The per-target audit on the shipped
    # section found 42/78 gain and 12/45 loss targets tightenable, 0 of 75
    # verdict changes; identical donor geometry and counted sets transfer this
    # finding verbatim to every pseudo-section, so intervals here are certified
    # supersets under the same uncapped construction as the shipped tables.
    diag = dict(capacity_note="uncapped construction (certified supersets); "
                              "per-target capacity audit transfers: geometry and "
                              "counted sets identical to shipped section",
                band_tested=int(len(band_tested)),
                violations=int(ver.violations.sum()),
                configs=int(len(ver) * 200), n_rows=int(len(res)),
                tested_mass=int(in_tested.sum()),
                tested_share=100 * float(in_tested.mean()),
                label_top=pd.Series(labels).value_counts().head(8).to_dict())
    return res, diag


# ---------------------------------------------------------------------------
# matched null: shipped b03_robust_null_v3 semantics, single-gather fast path
# (spec S7.1b). The fast path restructures ONLY the data movement: one gather
# of lab[donor_flat] per permutation and per-gene reduceat over the gathered
# array (bitwise-compatible with the shipped logical_or.reduceat, including
# its behavior on zero-donor rows). Seeds, block order, rng consumption,
# tie conventions, N computation, and output schema are identical to the
# shipped functions. Equivalence is enforced three ways: a fleet-gate
# validation mode (--mode validate-null, 20 named permutations on the REAL
# section vs the shipped functions, bitwise), a per-run perm-0 cross-check
# against the shipped function, and the shipped observed-interval assertion.
# ---------------------------------------------------------------------------

_W2 = {}


def worker_init_fast(Gmat, mov_cell, gene_slices_list, donor_flat, blocks,
                     obs_arr, cid):
    """Shipped-form state in, derived structures built worker-side once
    (spec S7.1b: spawn payload stays small; offsets never pickled).
    gene_slices_list entries: (gi, mv, absolute start offsets, [span_end])."""
    _W2["Gmat"] = Gmat
    _W2["mov_cell"] = mov_cell
    gene_rows = {}
    for gi, mv, ps, pe in gene_slices_list:
        offs = (ps - ps[0]).astype(np.int64) if len(mv) else np.zeros(0, np.int64)
        ps0 = int(ps[0]) if len(mv) else 0
        pe_last = int(pe[-1]) if len(mv) else 0
        gene_rows[gi] = (mv, offs, ps0, pe_last, mov_cell[mv])
    _W2["gene_rows"] = gene_rows
    _W2["donor_flat"] = donor_flat
    _W2["blocks"] = blocks
    _W2["obs"] = obs_arr
    _W2["cid"] = cid


def permute_and_score_fast(p):
    import b03_robust_null_v3 as NV  # for SEED_BASE only; sampling is local
    cid = _W2["cid"]
    lab_perm = cid.copy()
    rng = np.random.default_rng(NV.SEED_BASE + p)
    for idx in _W2["blocks"]:
        pm = rng.permutation(len(idx))
        lab_perm[idx] = lab_perm[idx[pm]]
    d = lab_perm[_W2["donor_flat"]]  # the single big gather
    Gmat = _W2["Gmat"]
    obs_arr = _W2["obs"]
    # N per (unique target, gene) via one boolean matvec per target
    # (integer-exact; replaces ~150 Gmat row-gathers per permutation)
    targets_unique = sorted({t for o in obs_arr for t in (o[2], o[3])})
    gsum = {}
    for t in targets_unique:
        in_t = lab_perm == t
        gsum[t] = in_t.astype(np.int64) @ Gmat if in_t.any() else \
            np.zeros(Gmat.shape[1], dtype=np.int64)
    n_keys = len(obs_arr)
    Tlo = np.empty(n_keys)
    Thi = np.empty(n_keys)
    # per-(gene, target) bounds are SHARED across keys (same geometry, same
    # counted set); compute once per unique pair instead of once per key
    # (~150 keys -> <= 42 unique pairs; the seg comparison over the gene's
    # donor segment is the dominant per-perm cost)
    cache = {}
    for k in range(n_keys):
        giL, giR, tS, tR, nLc, nRc, T_lo_o, T_hi_o = obs_arr[k]
        bounds = []
        for gi, target in ((giL, tS), (giR, tR)):
            key = (gi, target)
            if key not in cache:
                mv, offs, ps0, pe_last, cells_m = _W2["gene_rows"][gi]
                # unassigned movers (cell == -1) are never in a counted set
                src_in = np.where(cells_m >= 0,
                                  lab_perm[np.maximum(cells_m, 0)],
                                  -999) == target
                # boolean mask FIRST, then OR-reduceat (bitwise quirk parity
                # with the shipped logical_or.reduceat path)
                seg = d[ps0:pe_last] == target
                anyin = (np.logical_or.reduceat(seg, offs)
                         if len(offs) and len(seg) >= len(offs)
                         else np.zeros(len(mv), bool))
                cache[key] = (int(src_in.sum()),
                              int(((~src_in) & anyin).sum()))
            dmin, dmax = cache[key]
            N = int(gsum[target][gi])
            bounds.append((N, dmin, dmax))
        (NL, dLmin, dLmax), (NR, dRmin, dRmax) = bounds
        Tlo[k] = (math.log2(max(0, NL - dLmin) / nLc + EPS) +
                  math.log2(max(0, NR - dRmin) / nRc + EPS))
        Thi[k] = (math.log2((NL + dLmax) / nLc + EPS) +
                  math.log2((NR + dRmax) / nRc + EPS))
    return p, Tlo, Thi


def build_null_state(codes, uniq, labels, region_of_cell, res):
    """Assemble the null computation state for one pseudo-section (identical
    construction to the shipped main; obs cross-checked vs scout_bounds)."""
    tx = SYS["tx"]
    dm = SYS["dm"]
    cells_meta = SYS["cells_aligned"]
    n_cells = len(cells_meta)
    gidx = {g: i for i, g in enumerate(TESTED)}
    name_code = {g: i for i, g in enumerate(uniq.tolist())}
    cell = tx.cell_idx.to_numpy()
    Gmat = np.zeros((n_cells, len(TESTED)), dtype=np.int64)
    for g in TESTED:
        sel_g = (codes == name_code[g]) & (cell >= 0)
        Gmat[:, gidx[g]] = np.bincount(cell[sel_g], minlength=n_cells)
    types_needed = sorted(set(TYPE_MAP.values()) | set(RECV_MAP.values()))
    n_types = len(types_needed)
    SHIFT = 1
    type_shift = {t: i + SHIFT for i, t in enumerate(types_needed)}
    n_slots = n_types + SHIFT
    cid = region_of_cell * n_slots + np.array(
        [type_shift.get(l, 0) for l in labels], dtype=np.int64)
    tested_codes = np.array([name_code[g] for g in TESTED], dtype=np.int64)
    band_tested_full = dm["band_tested_indices"]
    keep = np.isin(codes[band_tested_full], tested_codes)
    band_tested = band_tested_full[keep]
    donor_sets = [dm["donor_sets"][i] for i in np.flatnonzero(keep)]
    mov_cell = cell[band_tested]
    lens = np.array([len(d) for d in donor_sets])
    donor_flat = np.concatenate([np.asarray(d, dtype=np.int64) for d in donor_sets
                                 if len(d)])
    donor_ptr = np.r_[0, np.cumsum(lens)]
    mov_gid = np.searchsorted(TESTED, uniq[codes[band_tested]])
    gene_rows = {}
    for gi in range(len(TESTED)):
        mv = np.flatnonzero(mov_gid == gi)
        ps, pe = donor_ptr[mv], donor_ptr[mv + 1]
        offs = ps - ps[0] if len(mv) else np.zeros(0, np.int64)
        ps0 = ps[0] if len(mv) else 0
        pe_last = pe[-1] if len(mv) else 0
        gene_rows[gi] = (mv, offs, ps0, pe_last, mov_cell[mv])
    base_tc = cells_meta.transcript_counts.to_numpy(np.int64)
    dec = np.zeros(n_cells, dtype=np.int64)
    for rr in range(3):
        idx_r = np.flatnonzero(region_of_cell == rr)
        q = np.quantile(base_tc[idx_r], np.linspace(0, 1, 11))
        dec[idx_r] = np.clip(np.digitize(base_tc[idx_r], q[1:-1]), 0, 9)
    blocks = []
    for rr in range(3):
        for dd in range(10):
            idx = np.flatnonzero((region_of_cell == rr) & (dec == dd))
            if len(idx) > 1:
                blocks.append(idx)

    # observed intervals, single code path (the shipped interval_for),
    # cross-checked against scout_bounds before any permutation runs
    keys, obs_arr = [], []
    for _, r in res.iterrows():
        Lg, Rg = r.pair.split("->")
        region_id = int(r.region)
        tS = region_id * n_slots + type_shift[TYPE_MAP[Lg]]
        tR = region_id * n_slots + type_shift[RECV_MAP[Rg]]
        NL, dLmin, dLmax = _interval_obs(Gmat, mov_cell, gene_rows, donor_flat,
                                         gidx[Lg], tS, cid)
        NR, dRmin, dRmax = _interval_obs(Gmat, mov_cell, gene_rows, donor_flat,
                                         gidx[Rg], tR, cid)
        nLc, nRc = int(r.nl_cells), int(r.nr_cells)
        T_lo = (math.log2(max(0, NL - dLmin) / nLc + EPS) +
                math.log2(max(0, NR - dRmin) / nRc + EPS))
        T_hi = (math.log2((NL + dLmax) / nLc + EPS) +
                math.log2((NR + dRmax) / nRc + EPS))
        worst = max(abs(T_lo - float(r.T_lo)), abs(T_hi - float(r.T_hi)))
        assert worst < 1e-6, f"observed interval mismatch {r.pair} r{region_id}: {worst}"
        keys.append((r.pair, region_id))
        obs_arr.append((gidx[Lg], gidx[Rg], tS, tR, nLc, nRc,
                        float(r.T_lo), float(r.T_hi)))
    # dtype diet (spec S7.1b): donor indices int32, labels int16 — the gather
    # is bandwidth-bound, this halves/quarters its traffic; values bounded
    # (slot indices < 2^31; cid < 3*n_slots+1 << 2^15). Slice form for the
    # worker: (gi, mv, absolute start offsets, [span_end]).
    donor_flat_32 = donor_flat.astype(np.int32)
    cid16 = cid.astype(np.int16)
    slices_list = []
    for gi in range(len(TESTED)):
        mv, offs, ps0, pe_last, _cells = gene_rows[gi]
        slices_list.append((gi, mv, (ps0 + offs).astype(np.int32),
                            np.array([pe_last], dtype=np.int64)))
    return (Gmat, mov_cell, slices_list, donor_flat_32, blocks, obs_arr, cid16,
            keys)


def _interval_obs(Gmat, mov_cell, gene_rows, donor_flat, gi, target, lab):
    mv, offs, ps0, pe_last, cells_m = gene_rows[gi]
    # unassigned movers (cell == -1) are never in a counted set (shipped guard)
    src_in = (np.where(cells_m >= 0, lab[np.maximum(cells_m, 0)], -999) == target) \
        if len(mv) else np.zeros(0, bool)
    if len(mv) and len(donor_flat):
        seg = (lab[donor_flat] == target)[ps0:pe_last]
        anyin = np.logical_or.reduceat(seg, offs) if len(offs) else \
            np.zeros(len(mv), bool)
    else:
        anyin = np.zeros(len(mv), bool)
    dmax = int(((~src_in) & anyin).sum())
    dmin = int(src_in.sum())
    N = int(Gmat[np.flatnonzero(lab == target), gi].sum())
    return N, dmin, dmax


def run_null_inproc(run_dir: Path, B: int, codes, uniq):
    """Matched null via the validated single-gather fast path (spec S7.1b).
    Semantics identical to shipped: per-perm seeds rng(NV.SEED_BASE + p), block
    order, rng consumption, tie conventions, N computation, count outputs,
    BH within region via the shipped implementation. Equivalence enforcement:
    fleet-gate validation on the REAL section (mode validate-null) plus a
    per-run bitwise perm-0 cross-check against the shipped function."""
    import b03_robust_null_v3 as NV  # shipped functions: cross-check + BH

    labels = np.load(run_dir / "labels.npy", allow_pickle=True)
    region_of_cell = np.load(DATA / "regions.npy")
    res = pd.read_csv(run_dir / "scout_bounds.csv")
    (Gmat, mov_cell, gene_rows, donor_flat, blocks, obs_arr, cid, keys) = \
        build_null_state(codes, uniq, labels, region_of_cell, res)

    # shipped-path state (gene_slices form) for the perm-0 cross-check
    dm = SYS["dm"]
    name_code = {g: i for i, g in enumerate(uniq.tolist())}
    tested_codes = np.array([name_code[g] for g in TESTED], dtype=np.int64)
    bti = dm["band_tested_indices"]
    keep2 = np.isin(codes[bti], tested_codes)
    dsets2 = [dm["donor_sets"][i] for i in np.flatnonzero(keep2)]
    lens2 = np.array([len(d) for d in dsets2])
    donor_ptr2 = np.r_[0, np.cumsum(lens2)]
    mov_cell2 = SYS["tx"].cell_idx.to_numpy()[bti[keep2]]
    mg2 = np.searchsorted(TESTED, uniq[codes[bti[keep2]]])
    gene_slices = {gi: (np.flatnonzero(mg2 == gi),
                        donor_ptr2[np.flatnonzero(mg2 == gi)],
                        donor_ptr2[np.flatnonzero(mg2 == gi) + 1])
                   for gi in range(len(TESTED))}

    # pre-flight: shipped vs fast path BITWISE on perm 0 (per-run), then the
    # single-gather sweep (shipped R2 seeds; R4 cross-check strengthened)
    ctx = mp.get_context("spawn")
    fast_args = (Gmat, mov_cell, gene_rows, donor_flat, blocks, obs_arr, cid)
    shipped_args = (Gmat, mov_cell2, mg2, donor_flat, gene_slices, blocks, obs_arr,
                    cid, None)  # caps: breast verified non-binding (shipped default)
    pool = ctx.Pool(int(os.environ.get("B03_STUDY_WORKERS", "4")),
                    initializer=worker_init_fast, initargs=fast_args)
    NV.worker_init(*shipped_args)
    worker_init_fast(*fast_args)
    ref_shipped = NV.permute_and_score(0)
    ref_fast = permute_and_score_fast(0)
    assert np.array_equal(ref_shipped[1], ref_fast[1]) and \
        np.array_equal(ref_shipped[2], ref_fast[2]), \
        "fast-path cross-check failed vs shipped implementation (perm 0)"
    w1 = pool.apply(permute_and_score_fast, (0,))
    assert np.array_equal(ref_fast[1], w1[1]) and np.array_equal(ref_fast[2], w1[2]), \
        "pre-flight determinism failed (parent vs worker)"
    Tlo_all = np.full((B, len(keys)), np.nan)
    Thi_all = np.full((B, len(keys)), np.nan)
    Tlo_all[0], Thi_all[0] = ref_fast[1], ref_fast[2]
    done = {0}
    t_last = time.time()
    for p, Tlo, Thi in pool.imap_unordered(permute_and_score_fast,
                                           range(1, B), chunksize=4):
        Tlo_all[p], Thi_all[p] = Tlo, Thi
        done.add(p)
        if len(done) % 25 == 0:
            rate = (time.time() - t_last) / 25
            t_last = time.time()
            log(f"  null {len(done)}/{B} perms | {rate:.1f}s/perm eff | "
                f"ETA {rate*(B-len(done))/60:.0f} min", run_dir.name)
    pool.close()
    pool.join()
    assert not np.isnan(Tlo_all).any(), "missing permutations"

    T_lo_obs = np.array([o[6] for o in obs_arr])
    T_hi_obs = np.array([o[7] for o in obs_arr])
    p_pos_acc = (Tlo_all >= T_lo_obs[None, :] - 1e-9).sum(axis=0)
    p_neg_acc = (Thi_all <= T_hi_obs[None, :] + 1e-9).sum(axis=0)
    rows = [dict(pair=k[0], region=k[1],
                 p_pos=p_pos_acc[i] / B, p_neg=p_neg_acc[i] / B,
                 T_lo_perm_med=float(np.median(Tlo_all[:, i])),
                 T_hi_perm_med=float(np.median(Thi_all[:, i])))
            for i, k in enumerate(keys)]
    nul = pd.DataFrame(rows)
    nul["q_pos"] = NV.bh_within_region(nul, "p_pos")
    nul["q_neg"] = NV.bh_within_region(nul, "p_neg")
    nul.to_csv(run_dir / "null_robust_matched.csv", index=False)
    (run_dir / "scout_summary_matched.json").write_text(json.dumps(dict(
        method="synthstudy in-process reuse of b03_robust_null_v3 (identical "
               "per-perm seeds, blocks, tie convention)",
        B=B, n_workers=int(os.environ.get("B03_STUDY_WORKERS", "4"))), indent=2))


# ---------------------------------------------------------------------------
# adjudication (shipped watchdog convention: add-one + BH within region)
# ---------------------------------------------------------------------------

def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def adjudicate(run_dir: Path):
    res = pd.read_csv(run_dir / "scout_bounds.csv")
    nul = pd.read_csv(run_dir / "null_robust_matched.csv")
    B = json.loads((run_dir / "scout_summary_matched.json").read_text())["B"]
    acc_pos = np.round(nul.p_pos * B).astype(int)
    acc_neg = np.round(nul.p_neg * B).astype(int)
    assert (np.abs(acc_pos - nul.p_pos * B) < 1e-6).all()
    assert (np.abs(acc_neg - nul.p_neg * B) < 1e-6).all()
    nul["p_pos_raw"] = (1 + acc_pos) / (B + 1)
    nul["p_neg_raw"] = (1 + acc_neg) / (B + 1)
    nul["q_pos"], nul["q_neg"] = np.nan, np.nan
    for reg in sorted(nul.region.unique()):
        sel = nul.region == reg
        nul.loc[sel, "q_pos"] = bh(nul.loc[sel, "p_pos_raw"])
        nul.loc[sel, "q_neg"] = bh(nul.loc[sel, "p_neg_raw"])
    final = res.merge(nul[["pair", "region", "p_pos_raw", "p_neg_raw",
                           "q_pos", "q_neg"]], on=["pair", "region"], how="inner")
    final["certified_pos_matched"] = (final.T_lo > 0) & (final.q_pos <= 0.10)
    final["certified_neg_matched"] = (final.T_hi < 0) & (final.q_neg <= 0.10)
    final.to_csv(run_dir / "scout_final_matched.csv", index=False)
    return final


# ---------------------------------------------------------------------------
# modes
# ---------------------------------------------------------------------------

def do_validate_null(n_perms=20):
    """Fleet gate (spec S7.1b): on the REAL section (base data, no pi), the
    fast path must reproduce the shipped functions BITWISE on named
    permutations. Writes results/null_path_validation.csv."""
    load_base()
    import b03_robust_null_v3 as NV
    labels = np.load(DATA / "labels.npy", allow_pickle=True)
    region_of_cell = np.load(DATA / "regions.npy")
    res = pd.read_csv(OUT / "scout_bounds.csv")
    tx = SYS["tx"]
    codes, uniq = pd.factorize(tx.feature_name, sort=False)
    codes = np.asarray(codes, dtype=np.int64)
    uniq = np.asarray(uniq, dtype=object)
    (Gmat, mov_cell, gene_rows, donor_flat, blocks, obs_arr, cid, keys) = \
        build_null_state(codes, uniq, labels, region_of_cell, res)
    # shipped state
    dm = SYS["dm"]
    name_code = {g: i for i, g in enumerate(uniq.tolist())}
    tested_codes = np.array([name_code[g] for g in TESTED], dtype=np.int64)
    bti = dm["band_tested_indices"]
    keep2 = np.isin(codes[bti], tested_codes)
    dsets2 = [dm["donor_sets"][i] for i in np.flatnonzero(keep2)]
    lens2 = np.array([len(d) for d in dsets2])
    donor_ptr2 = np.r_[0, np.cumsum(lens2)]
    mov_cell2 = tx.cell_idx.to_numpy()[bti[keep2]]
    mg2 = np.searchsorted(TESTED, uniq[codes[bti[keep2]]])
    gene_slices = {gi: (np.flatnonzero(mg2 == gi),
                        donor_ptr2[np.flatnonzero(mg2 == gi)],
                        donor_ptr2[np.flatnonzero(mg2 == gi) + 1])
                   for gi in range(len(TESTED))}
    NV.worker_init(Gmat, mov_cell2, mg2, donor_flat, gene_slices, blocks,
                   obs_arr, cid, None)
    worker_init_fast(Gmat, mov_cell, gene_rows, donor_flat, blocks, obs_arr, cid)
    rows = []
    for p in range(n_perms):
        a = NV.permute_and_score(p)
        b = permute_and_score_fast(p)
        eq = bool(np.array_equal(a[1], b[1]) and np.array_equal(a[2], b[2]))
        rows.append(dict(perm=p, shipped_Tlo_min=float(a[1].min()),
                         fast_Tlo_min=float(b[1].min()),
                         shipped_Thi_max=float(a[2].max()),
                         fast_Thi_max=float(b[2].max()),
                         bitwise_equal=eq))
        assert eq, f"fast path diverges from shipped at perm {p}"
    pd.DataFrame(rows).to_csv(OUT / "null_path_validation.csv", index=False)
    print(f"VALIDATE-NULL: {n_perms}/{n_perms} permutations bitwise identical "
          f"(shipped vs single-gather fast path). Fleet gate PASSED.")


def do_verify():
    load_base()
    feat = SYS["tx"].feature_name.to_numpy()
    gcount = pd.Series(feat).value_counts()
    tmass = int(gcount.reindex(TESTED).fillna(0).sum())
    print(f"VERIFY baseline tested-gene mass share: {100*tmass/len(feat):.2f}% "
          f"({tmass:,} of {len(feat):,})")
    # one-time label-vote recomputation (spec S3.5, run once per fleet)
    codes_v, uniq_v = pd.factorize(SYS["tx"].feature_name, sort=False)
    labels_re = compute_labels(np.asarray(codes_v, dtype=np.int64),
                               np.asarray(uniq_v, dtype=object))
    assert (labels_re == SYS["labels"]).all(), "fleet labels diverge from shipped"
    print("VERIFY label-vote recomputation matches shipped labels (one-time)")
    anchor_tx = tx_in_anchor(proximity_mask(911001))
    for cls in INJ_CLASSES:
        pool = tested_class_pool(cls, exclude=dict(zip(INJ_CLASSES, INJ_PAIR)).get(cls, set()))
        print(f"VERIFY class {cls!r}: borrow pool {pool}")
    dm = SYS["dm"]
    assert len(dm["band_tested_indices"]) == len(dm["donor_sets"])
    print(f"VERIFY donor map: {len(dm['band_tested_indices']):,} band-tested "
          f"transcripts; geometry rebuilt per run over the stored set")
    print("VERIFY OK")


def do_run(run_id, arm, target_pct=None, seed=None, f=None, B=250):
    load_base()
    if arm == "N":
        run = arm_null(seed)
    elif arm == "E":
        run = arm_effect(seed, f)
    elif arm == "M":
        run = arm_mass(target_pct, seed)
    else:
        raise ValueError(arm)
    run_dir = STUDY / run_id
    if run_dir.exists():
        raise RuntimeError(f"run dir exists: {run_dir}")
    run_dir.mkdir(parents=True)
    t0 = time.time()
    log("building pseudo-section", run_id)
    codes, uniq, info = build_pseudo_section(run)
    log(f"pseudo-section built in {time.time()-t0:.0f}s "
        f"(tested share {info.get('tested_share', 'n/a') if 'tested_share' in info else ''})"
        , run_id)
    res, diag = run_intervals(codes, uniq, run_dir)
    log(f"intervals done: {diag['n_rows']} rows; violations "
        f"{diag['violations']}/{diag['configs']}; band_tested {diag['band_tested']:,}; "
        f"tested share {diag['tested_share']:.2f}%", run_id)
    run_null_inproc(run_dir, B, codes, uniq)
    log(f"null done (B={B}); adjudicating", run_id)
    final = adjudicate(run_dir)
    npos = int(final.certified_pos_matched.sum())
    nneg = int(final.certified_neg_matched.sum())
    inj_rows = None
    if run["injection"]:
        key = f"{run['injection']['pair'][0]}->{run['injection']['pair'][1]}"
        inj_rows = final[final.pair == key].to_dict("records")
    manifest = dict(run_id=run_id, spec=run["spec"], info=info, diag=diag, B=B,
                    n_cert_pos=npos, n_cert_neg=nneg, inj_pair_rows=inj_rows)
    (run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    log(f"RUN COMPLETE {run_id} in {(time.time()-t0)/60:.1f} min: "
        f"{npos} cert-pos / {nneg} cert-neg / {len(final)-npos-nneg} non-id", run_id)


def fleet_list():
    # spec S7.4/S7.8 final trim: E (f in {1.0, 4.0} x 3 anchors, dose-response
    # design) -> M (4 targets x 2 seeds) -> N (6 replicates; FCR saturated)
    runs = []
    for s in (911001, 911002, 911003):
        for f in (1.0, 4.0):
            runs.append(dict(run_id=f"E{s % 1000:03d}F{int(f*10)}", arm="E",
                             seed=s, f=f))
    mseed = itertools.count(912001)
    for t in M_TARGETS:
        for _k in (1, 2):
            runs.append(dict(run_id=f"M{int(t*10):03d}{_k}", arm="M",
                             target_pct=t, seed=next(mseed)))
    runs += [dict(run_id=f"N{i:02d}", arm="N", seed=910000 + i)
             for i in range(1, 7)]
    return runs


def do_fleet():
    runs = fleet_list()
    assert len(runs) == N_RUNS_TOTAL
    t_start = time.time()
    log(f"fleet: {len(runs)} runs queued (B=250 uniform)")
    for i, r in enumerate(runs):
        if (STUDY / r["run_id"] / "run_manifest.json").exists():
            log(f"skip {r['run_id']} (manifest exists)")
            continue
        if (r["arm"] == "N" and r["run_id"] > "N05" and len(runs) > 17):
            elapsed = time.time() - t_start
            done_frac = i / len(runs)
            projected = elapsed / max(done_frac, 1e-9)
            if projected > 3.5 * 3600:
                log(f"TRIM RULE: projected {projected/3600:.1f}h > 3.5h; "
                    f"skipping Arm-N replicate {r['run_id']} (5 replicates "
                    f"saturate FCR reporting)")
                continue
        log(f"fleet launch {i+1}/{len(runs)}: {r}")
        try:
            do_run(r["run_id"], r["arm"], target_pct=r.get("target_pct"),
                   seed=r.get("seed"), f=r.get("f"), B=250)
        except Exception as e:
            log(f"RUN FAILED {r['run_id']}: {e}; fleet continues")


def do_analyze():
    rows = []
    for mf in sorted(STUDY.glob("*/run_manifest.json")):
        m = json.loads(mf.read_text())
        spec = m["spec"]
        diag = m["diag"]
        row = dict(run_id=m["run_id"], arm=spec["arm"], B=m["B"],
                   tested_share=round(diag.get("tested_share", float("nan")), 3),
                   n_rows=diag["n_rows"],
                   n_cert_pos=m["n_cert_pos"], n_cert_neg=m["n_cert_neg"],
                   n_nonid=diag["n_rows"] - m["n_cert_pos"] - m["n_cert_neg"],
                   violations=diag["violations"])
        if spec["arm"] == "M":
            row["target_pct"] = spec["target_pct"]
        if m.get("inj_pair_rows"):
            r0 = m["inj_pair_rows"][0]
            row["inj_cert_any"] = bool(r0["certified_pos_matched"])
            row["inj_T_lo"] = float(r0["T_lo"])
            row["inj_q_pos"] = float(r0["q_pos"])
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(STUDY / "fleet_summary.csv", index=False)
    print(out.to_string(index=False))
    n_arm = out[out.arm == "N"]
    if len(n_arm):
        fcr_run = (n_arm.n_cert_pos + n_arm.n_cert_neg) / n_arm.n_rows
        print(f"\nV1 known-null FCR: per-run max {fcr_run.max():.4f} (gate <= 0.15); "
              f"pooled {n_arm.n_cert_pos.sum() + n_arm.n_cert_neg.sum()}/"
              f"{n_arm.n_rows.sum()} = "
              f"{(n_arm.n_cert_pos.sum() + n_arm.n_cert_neg.sum()) / n_arm.n_rows.sum():.4f}")
    m_arm = out[out.arm == "M"]
    if len(m_arm):
        print("\nV3 power curve (injected pair by target mass share):")
        print(m_arm.groupby("target_pct")[["inj_cert_any", "inj_T_lo"]].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["verify", "validate-null", "smoke", "run",
                             "fleet", "analyze"])
    ap.add_argument("--run-id")
    ap.add_argument("--arm")
    ap.add_argument("--target", type=float)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--f", type=float)
    ap.add_argument("--B", type=int, default=200)
    a = ap.parse_args()
    if a.mode == "verify":
        do_verify()
    elif a.mode == "validate-null":
        do_validate_null()
    elif a.mode == "smoke":
        do_run("SMOKE_N01", "N", seed=910001, B=10)
    elif a.mode == "run":
        do_run(a.run_id, a.arm, target_pct=a.target, seed=a.seed, f=a.f, B=a.B)
    elif a.mode == "fleet":
        do_fleet()
    elif a.mode == "analyze":
        do_analyze()


if __name__ == "__main__":
    main()
