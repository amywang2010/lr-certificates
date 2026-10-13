"""B03 Phase 7: proseg completion watcher — validates outputs the moment proseg exits.

Polls every 60 s for the proseg process. On exit, verifies the DEV amendment 010 output
naming (cells/, genes/, counts/ under the workspace xenium_breast_s6/crop/proseg_out/),
row-count sanity vs the crop, and writes proseg_done.json. Does NOT run
analyses: E2/E3/E4 run as their own audited step after this gate.
"""
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
OUT = os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "proseg_out")
DONE = os.path.join(WORKSPACE, "xenium_breast_s6", "crop", "proseg_done.json")
N_CROP_CELLS = 65964
N_CROP_TX = 4912129  # is_gene & QV>=20 (the only rows proseg consumes)

def proseg_alive():
    out = subprocess.run(["powershell", "-NoProfile", "-Command",
                          "(Get-Process proseg -ErrorAction SilentlyContinue) -ne $null"],
                         capture_output=True, text=True, timeout=60)
    return out.stdout.strip().lower() == "true"

def main():
    t0 = time.time()
    while proseg_alive():
        time.sleep(60)
    print(f"[{time.time()-t0:.0f}s] proseg exited; validating outputs", flush=True)
    res = {"proseg_exited": time.strftime("%Y-%m-%d %H:%M:%S"),
           "elapsed_s": round(time.time() - t0)}
    try:
        cells = os.listdir(f"{OUT}/cells")
        genes = os.listdir(f"{OUT}/genes")
        counts = os.listdir(f"{OUT}/counts")
        res["outputs_present"] = {"cells": cells, "genes": genes, "counts": counts}
        n_seg = None
        for f in cells:
            if f.endswith(".parquet"):
                import pyarrow.parquet as pq
                n_seg = pq.ParquetFile(f"{OUT}/cells/{f}").metadata.num_rows
        res["n_proseg_cells"] = n_seg
        res["sanity"] = {
            "proseg_cells_ge_crop": bool(n_seg and n_seg >= 0.5 * N_CROP_CELLS),
            "structure_complete": bool(cells and genes and counts),
        }
        res["status"] = "OK" if all(res["sanity"].values()) else "SUSPECT"
    except Exception as e:  # noqa: BLE001 - watcher records, never crashes silently
        res["status"] = "ERROR"
        res["error"] = f"{type(e).__name__}: {e}"
    with open(DONE, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1), flush=True)
    sys.exit(0 if res.get("status") == "OK" else 2)

if __name__ == "__main__":
    main()
