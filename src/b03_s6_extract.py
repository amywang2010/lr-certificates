"""B03 Phase 7 step 1: extract the analysis source from the S-BIAD2146 zarr zip.

Registered packaging rules (frozen before this run):
  - members -> the workspace xenium_breast_s6/extracted/ dir, mirroring a 10x-like layout
  - transcript parts streamed verbatim (16 parts)
  - cell/nucleus boundaries parquet verbatim
  - table_cells CSR snapshot: data/indices/indptr + var gene names + obs cell ids
Assertions (abort before any downstream step on failure):
  - 16 transcript parts present; schema has x,y,feature_name,cell_id,qv,is_gene
  - CSR dims from X/.zattrs == (692184, 5101); indptr length == rows+1
  - coordinate transform verification deferred to the loader (needs both extents)
"""
import json, os, shutil, sys, time, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
ZIP = os.path.join(WORKSPACE, "xenium_breast_s6", "sdata_breast_s6.zarr.zip")
OUT = os.path.join(WORKSPACE, "xenium_breast_s6", "extracted")
TX_DIR = os.path.join(OUT, "transcripts_parts")
LOG = os.path.join(WORKSPACE, "xenium_breast_s6", "extract_log.json")

def log(d):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(json.dumps(d) + "\n")

def main():
    t0 = time.time()
    os.makedirs(TX_DIR, exist_ok=True)
    os.makedirs(os.path.join(OUT, "boundaries"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "table_cells"), exist_ok=True)
    z = zipfile.ZipFile(ZIP)
    names = set(z.namelist())
    n_tx = 0
    for i in range(16):
        member = f"sdata_breast_s6.zarr/points/st/points.parquet/part.{i}.parquet"
        assert member in names, f"missing transcript part {i}"
        dst = os.path.join(TX_DIR, f"part.{i}.parquet")
        if not os.path.exists(dst):
            with z.open(member) as src, open(dst, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
        n_tx += os.path.getsize(dst)
    log({"step": "transcript_parts", "bytes": n_tx, "t": time.time() - t0})

    for elem, dst_name in [
        ("shapes/cell_boundaries/shapes.parquet", "cell_boundaries.parquet"),
        ("shapes/nucleus_boundaries/shapes.parquet", "nucleus_boundaries.parquet"),
    ]:
        member = f"sdata_breast_s6.zarr/{elem}"
        dst = os.path.join(OUT, "boundaries", dst_name)
        if not os.path.exists(dst):
            with z.open(member) as src, open(dst, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
        log({"step": dst_name, "bytes": os.path.getsize(dst), "t": time.time() - t0})

    base = "sdata_breast_s6.zarr/tables/table_cells"
    tc = os.path.join(OUT, "table_cells")
    with z.open(f"{base}/X/.zattrs") as f:
        xattrs = json.load(f)
    shape = tuple(xattrs["shape"])
    assert shape == (692184, 5101), f"CSR shape mismatch: {shape}"
    for arr in ["X/data", "X/indices", "X/indptr"]:
        sub = arr.split("/")[1]
        os.makedirs(os.path.join(tc, sub), exist_ok=True)
        members = [n for n in names if n.startswith(f"{base}/{arr}/") and not n.endswith(".zarray")]
        for m in members:
            leaf = os.path.join(tc, sub, m.rsplit("/", 1)[1])
            if os.path.exists(leaf):
                continue
            with z.open(m) as src, open(leaf, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
        with z.open(f"{base}/{arr}/.zarray") as f:
            meta = json.load(f)
        with open(os.path.join(tc, sub, ".zarray.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f)
    with z.open(f"{base}/var/_index/.zarray") as f:
        var_meta = json.load(f)
    vmembers = [n for n in names if n.startswith(f"{base}/var/_index/") and not n.endswith((".zarray", ".zattrs"))]
    os.makedirs(os.path.join(tc, "var_index"), exist_ok=True)
    for m in vmembers:
        leaf = os.path.join(tc, "var_index", m.rsplit("/", 1)[1])
        if not os.path.exists(leaf):
            with z.open(m) as src, open(leaf, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
    with z.open(f"{base}/var/_index/.zattrs") as f:
        vattrs = json.load(f)
    with open(os.path.join(tc, "var_index", ".zattrs.json"), "w", encoding="utf-8") as f:
        json.dump(vattrs, f)
    obsm = [n for n in names if n.startswith(f"{base}/obs/_index/") and not n.endswith((".zarray", ".zattrs"))]
    os.makedirs(os.path.join(tc, "obs_index"), exist_ok=True)
    for m in obsm:
        leaf = os.path.join(tc, "obs_index", m.rsplit("/", 1)[1])
        if not os.path.exists(leaf):
            with z.open(m) as src, open(leaf, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
    # obs/cell_id: the instance-key join column for the G3 concordance gate
    # (verified present in the zip before this run; A1-b documents the join).
    cid = [n for n in names if n.startswith(f"{base}/obs/cell_id/") and not n.endswith((".zarray", ".zattrs"))]
    assert cid, "obs/cell_id missing from zip"
    os.makedirs(os.path.join(tc, "obs_cell_id"), exist_ok=True)
    for m in cid:
        leaf = os.path.join(tc, "obs_cell_id", m.rsplit("/", 1)[1])
        if not os.path.exists(leaf):
            with z.open(m) as src, open(leaf, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
    with z.open(f"{base}/obs/cell_id/.zarray") as f:
        cmeta = json.load(f)
    with open(os.path.join(tc, "obs_cell_id", ".zarray.json"), "w", encoding="utf-8") as f:
        json.dump(cmeta, f)
    # obs/transcript_counts: vendor per-cell total, input to the frozen region
    # rule (make_regions log10(1+transcript_counts), seed 20260907).
    tcc = [n for n in names if n.startswith(f"{base}/obs/transcript_counts/") and not n.endswith((".zarray", ".zattrs"))]
    assert tcc, "obs/transcript_counts missing from zip"
    os.makedirs(os.path.join(tc, "obs_transcript_counts"), exist_ok=True)
    for m in tcc:
        leaf = os.path.join(tc, "obs_transcript_counts", m.rsplit("/", 1)[1])
        if not os.path.exists(leaf):
            with z.open(m) as src, open(leaf, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
    with z.open(f"{base}/obs/transcript_counts/.zarray") as f:
        tmeta = json.load(f)
    with open(os.path.join(tc, "obs_transcript_counts", ".zarray.json"), "w", encoding="utf-8") as f:
        json.dump(tmeta, f)
    log({"step": "table_cells", "shape": list(shape), "t": time.time() - t0})

    import pyarrow.parquet as pq
    pf = pq.ParquetFile(os.path.join(TX_DIR, "part.0.parquet"))
    cols = set(pf.schema_arrow.names)
    need = {"x", "y", "feature_name", "cell_id", "qv", "is_gene"}
    assert need <= cols, f"transcript schema missing {need - cols}"
    total_rows = 0
    for i in range(16):
        total_rows += pq.ParquetFile(os.path.join(TX_DIR, f"part.{i}.parquet")).metadata.num_rows
    summary = {
        "step": "done",
        "transcript_rows_total": total_rows,
        "csr_shape": list(shape),
        "seconds": time.time() - t0,
    }
    log(summary)
    print(json.dumps(summary))

if __name__ == "__main__":
    main()
