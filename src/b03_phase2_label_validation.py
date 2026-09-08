"""Label-transfer accuracy validation (2-fold CV on vendor cells).

Question (referee-proofing the Phase-2 methods): how accurate is nearest-centroid
cosine label transfer in the frozen marker space? Estimated by 2-fold CV on VENDOR
cells, whose labels are the scout's marker-vote labels (the same source the centroids
were built from). Reports overall + per-type accuracy and the confusion structure.

Uses the frozen transform constants from vendor_markers.npz (scale, marker panel,
centroids) so the validated classifier is EXACTLY the one applied to Proseg cells.
The h5 read repeats the layout logic of b03_phase2_vendor_prep.py (deterministic).
"""

import numpy as np
import pandas as pd
import h5py
import pickle
import sys
import time

sys.path.insert(0, "B03_project/src")
from b03_scout import ALL_PAIRS

DATA = "B03_project/data"
t0 = time.time()
tested = sorted({g for p in ALL_PAIRS for g in p})

vm = np.load(f"{DATA}/vendor_markers.npz", allow_pickle=True)
marker_genes = [str(g) for g in vm["marker_genes"]]
type_names = [str(t) for t in vm["type_names"]]
scale = float(vm["scale"])

# rebuild vendor marker matrix (same read as vendor_prep, deterministic)
labels_poly = np.load(f"{DATA}/labels.npy", allow_pickle=True)
with open(f"{DATA}/donor_map.pkl", "rb") as f:
    dm = pickle.load(f)
cells_meta = pd.read_parquet(f"{DATA}/cells_meta.parquet")
cid_map = {str(k): v for k, v in dm["cid_map"].items()}
cm_ids = cells_meta.cell_id.astype(str).to_numpy()
poly_row = np.full(max(cid_map.values()) + 1, -1, dtype=np.int64)
for row, cid in enumerate(cm_ids):
    p = cid_map.get(cid, -1)
    if p >= 0:
        poly_row[p] = row
labels_row = labels_poly[poly_row]

h5_path = f"{DATA}/Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5"
with h5py.File(h5_path, "r") as h5:
    features = [b.decode() if isinstance(b, bytes) else str(b)
                for b in h5["matrix/features/name"][:]]
    indptr = h5["matrix/indptr"][:]
    indices = h5["matrix/indices"][:]
    data = h5["matrix/data"][:].astype(np.int64, copy=False)
n_cells = len(cells_meta)
assert len(indptr) == n_cells + 1
n_markers = len(marker_genes)
M = np.zeros((n_cells, n_markers), dtype=np.int32)
feat_to_col = np.full(len(features), -1, dtype=np.int64)
for j, g in enumerate(marker_genes):
    feat_to_col[features.index(g)] = j
col_ids = np.repeat(np.arange(n_cells, dtype=np.int64), np.diff(indptr))
j = feat_to_col[indices]
sel = j >= 0
M[col_ids[sel], j[sel]] = data[sel]
print(f"[{time.time()-t0:.0f}s] vendor marker matrix rebuilt: {M.shape}", flush=True)

C = M.sum(axis=1).astype(np.int64)
X = np.log1p((M / np.maximum(C, 1)[:, None]) * scale).astype(np.float64)
Xn = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)

# evaluation set: cells whose vendor label is one of the 15 transfer targets
lab_arr = np.asarray(labels_row, dtype=object)
eval_mask = np.isin(lab_arr, type_names)
y = np.array([type_names.index(l) for l in lab_arr[eval_mask]])
Xe = Xn[eval_mask]
print(f"[{time.time()-t0:.0f}s] eval set: {len(y)} cells across {len(type_names)} types",
      flush=True)

# 2-fold CV (deterministic split)
rng = np.random.RandomState(20260906)
perm = rng.permutation(len(y))
folds = [perm[: len(y) // 2], perm[len(y) // 2:]]
accs = []
conf = np.zeros((len(type_names), len(type_names)), dtype=np.int64)
for k in range(2):
    tr, te = folds[1 - k], folds[k]
    cent = np.stack([Xe[tr][y[tr] == t].mean(axis=0) for t in range(len(type_names))])
    cent = cent / np.maximum(np.linalg.norm(cent, axis=1, keepdims=True), 1e-12)
    pred = (Xe[te] @ cent.T).argmax(axis=1)
    accs.append((pred == y[te]).mean())
    for a, b in zip(y[te], pred):
        conf[a, b] += 1
acc = float(np.mean(accs))
print(f"[{time.time()-t0:.0f}s] 2-fold CV accuracy: {acc*100:.1f}% "
      f"(folds: {accs[0]*100:.1f}%, {accs[1]*100:.1f}%)", flush=True)
print("per-type recall (vendor label -> transferred):")
for i, t in enumerate(type_names):
    n_i = conf[i].sum()
    if n_i:
        print(f"  {t:26s} {conf[i, i]}/{n_i} = {conf[i, i]/n_i*100:.1f}%")
top_conf = []
for a in range(len(type_names)):
    for b in range(len(type_names)):
        if a != b and conf[a, b] > 0.05 * max(conf[a].sum(), 1):
            top_conf.append((type_names[a], type_names[b],
                             conf[a, b], int(conf[a].sum())))
top_conf.sort(key=lambda x: -x[2])
print("largest confusions (true -> predicted, count, row total):")
for a, b, c, n in top_conf[:8]:
    print(f"  {a} -> {b}: {c}/{n}")

out = pd.DataFrame(conf, index=type_names, columns=type_names)
out.to_csv("B03_project/results/phase2_label_transfer_cv_confusion.csv")
with open("B03_project/results/phase2_label_transfer_cv.json", "w") as f:
    import json
    json.dump(dict(accuracy_2fold_cv=acc, n_eval=int(eval_mask.sum()),
                   n_types=len(type_names), transform="frozen from vendor_markers.npz"),
              f, indent=2)
print("LABEL VALIDATION COMPLETE", flush=True)
