"""B03 Phase 5: proseg completion watcher — validates outputs the moment proseg exits.

Polls every 30 s. On proseg exit: verifies the amendment record naming (cells/genes/counts
under xenium_lung/crop/proseg_out/), row-count sanity vs the crop, and writes
proseg_done.json. Does NOT run analyses (E2 runs as its own audited step).
"""
import json, os, subprocess, sys, time

OUT = "xenium_lung/crop/proseg_out"
DONE = "xenium_lung/crop/proseg_done.json"

def proseg_alive():
    r = subprocess.run(["powershell", "-NoProfile", "-Command",
                        "if(Get-Process proseg -ErrorAction SilentlyContinue){'1'}else{'0'}"],
                       capture_output=True, text=True, timeout=60)
    return r.stdout.strip() == "1"

def validate():
    need = ["cells", "genes", "counts"]
    info = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"), "files": {}, "ok": False}
    for n in need:
        p = os.path.join(OUT, n)
        info["files"][n] = os.path.getsize(p) if os.path.exists(p) else None
    try:
        import pyarrow.parquet as pq
        cpath = os.path.join(OUT, "cells")
        if os.path.exists(cpath):
            md = pq.ParquetFile(cpath).metadata
            info["proseg_cells"] = md.num_rows
            cm = pq.ParquetFile("xenium_lung/crop/cells.parquet").metadata
            info["vendor_cells"] = cm.num_rows
            info["cell_count_delta"] = md.num_rows - cm.num_rows
        gpath = os.path.join(OUT, "genes")
        if os.path.exists(gpath):
            info["proseg_genes"] = pq.ParquetFile(gpath).metadata.num_rows
        info["ok"] = all(info["files"][n] for n in need)
    except Exception as e:  # noqa: BLE001
        info["error"] = f"{type(e).__name__}: {e}"
    return info

def main():
    print("watcher armed: polling proseg every 30 s", flush=True)
    while proseg_alive():
        time.sleep(30)
    print("proseg exited; validating outputs", flush=True)
    info = validate()
    with open(DONE, "w") as f:
        json.dump(info, f, indent=1)
    print(json.dumps(info, indent=1), flush=True)
    sys.exit(0 if info["ok"] else 2)

if __name__ == "__main__":
    main()
