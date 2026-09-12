"""B03 Phase 5 G3: vendor-assignment concordance — per B03_PHASE5_PREREG.md gate G3.

Verifies that per-cell per-gene counts derived from transcripts.parquet (cell_id
barcodes) match the vendor's cell_feature_matrix.h5 on a random 30-gene subset.

XOA 3.0 schema facts (verified from the file, not assumed — amendment record):
  - cell_id is a STRING BARCODE; unassigned = 'UNASSIGNED' (17.5% of rows,
    consistent with metrics fraction_transcripts_assigned = 0.8201).
  - h5 shape is [features, cells] = [10029, 278328]; gene row i -> counts H[i, :].
  - h5 barcodes[] gives the column order.

Memory architecture (v2): the naive single-pass to_pandas() hit ArrowMemoryError
(177M string rows). v2 iterates PARQUET ROW GROUPS with constant memory (~2 GB):
per chunk, numpy masks + one vectorized barcode->column remap + bincount into a
(30 x ncells) accumulator. Disk streams once (~2.4 GB x read amplification).
"""
import time

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.sparse import csc_matrix

TX = "xenium_lung/extracted/transcripts.parquet"
H5 = "xenium_lung/extracted/cell_feature_matrix.h5"
OUT = "B03_project/results/lung_concordance.csv"
SEED = 20260907
N_TEST_GENES = 30
UNASSIGNED = "UNASSIGNED"

t0 = time.time()


def log(msg):
    print(f"[{time.time()-t0:.0f}s] {msg}", flush=True)


def main():
    pf = pq.ParquetFile(TX)
    n_rowgroups = pf.metadata.num_row_groups
    log(f"parquet: {pf.metadata.num_rows:,} rows in {n_rowgroups} row groups")

    with h5py.File(H5, "r") as f:
        mat = f["matrix"]
        shape = mat["shape"][:]  # [features, cells]
        barcodes = np.array([b.decode() for b in mat["barcodes"][:]])
        features = np.array([b.decode() for b in mat["features"]["name"][:]])
        H = csc_matrix((mat["data"][:], mat["indices"][:], mat["indptr"][:]),
                       shape=tuple(shape))
    ncells = int(shape[1])
    log(f"h5 [features, cells]: {shape}, nnz {H.nnz:,}")

    barcode_pos = pd.Series(np.arange(ncells, dtype=np.int64), index=barcodes)
    is_ctrl = np.array([g.startswith(("Blank", "NegControl", "Unassigned",
                                       "antisense", "DeprecatedCodeword",
                                       "GenomicControl"))
                        or "Codeword" in g
                        for g in features])
    gene_pool = features[~is_ctrl]
    rng = np.random.default_rng(SEED)
    test_genes = np.sort(rng.choice(gene_pool, N_TEST_GENES, replace=False))
    gene_pos = {g: int(np.flatnonzero(features == g)[0]) for g in test_genes}
    gene_slot = {g: k for k, g in enumerate(test_genes)}
    acc = np.zeros((N_TEST_GENES, ncells), dtype=np.int64)
    totals = dict(assigned=0, unassigned=0)

    for rg in range(n_rowgroups):
        chunk = pf.read_row_group(rg, columns=["feature_name", "cell_id", "qv",
                                               "is_gene"])
        qv = chunk.column("qv").to_numpy()
        is_gene = np.asarray(chunk.column("is_gene").to_numpy(zero_copy_only=False),
                             dtype=bool)
        cid = chunk.column("cell_id").to_numpy(zero_copy_only=False)
        feat = chunk.column("feature_name").to_numpy(zero_copy_only=False)
        sel = (qv >= 20) & is_gene & (cid != UNASSIGNED)
        totals["assigned"] += int(sel.sum())
        totals["unassigned"] += int(((qv >= 20) & is_gene
                                     & (cid == UNASSIGNED)).sum())
        if not sel.any():
            continue
        ccols = barcode_pos.reindex(cid[sel]).to_numpy()  # -1 for unknown barcode
        known = ccols >= 0
        ccols = ccols[known]
        fsub = feat[sel][known]
        for k, g in enumerate(test_genes):
            m = fsub == g
            if m.any():
                acc[k] += np.bincount(ccols[m], minlength=ncells)
        if (rg + 1) % 25 == 0 or rg == n_rowgroups - 1:
            log(f"row group {rg+1}/{n_rowgroups} | assigned so far "
                f"{totals['assigned']:,}")

    log(f"QV20 assigned: {totals['assigned']:,}; unassigned: "
        f"{totals['unassigned']:,} "
        f"({totals['unassigned']/(totals['assigned']+totals['unassigned'])*100:.1f}%)")

    rows = []
    for k, g in enumerate(test_genes):
        counts_tx = acc[k]
        counts_h5 = np.asarray(H[gene_pos[g], :].todense()).ravel()
        tot_tx, tot_h5 = int(counts_tx.sum()), int(counts_h5.sum())
        absdiff = int(np.abs(counts_tx - counts_h5).sum())
        rows.append(dict(gene=g, tx_total=tot_tx, h5_total=tot_h5,
                         abs_diff=absdiff,
                         rel_diff=(tot_tx - tot_h5) / max(tot_h5, 1)))
    res = pd.DataFrame(rows)
    res.to_csv(OUT, index=False)
    print(res.to_string(), flush=True)
    mism = res[~((res.rel_diff.abs() < 0.01)
                 & (res.abs_diff <= np.maximum(3, 0.001 * res.h5_total)))]
    verdict = "PASS" if len(mism) == 0 else "CHECK: " + ",".join(mism.gene)
    print(f"G3 CONCORDANCE VERDICT: {len(res)-len(mism)}/{N_TEST_GENES} genes "
          f"concordant -> {verdict}", flush=True)


if __name__ == "__main__":
    main()
