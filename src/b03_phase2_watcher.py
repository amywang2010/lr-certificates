"""Phase-2 completion watcher: fires the (already fixture-validated) analysis the
moment proseg run 4 exits, with completion-vs-crash discrimination.

Normal exit: counts + cells.parquet + genes.parquet all present AND counts mtx
header/body consistent -> run b03_phase2_analysis.py, log result.
Crash/abort: outputs missing/incomplete -> log loudly, DO NOT run analysis.
Exactly-once: exits after firing (or after logging a crash); never re-runs.
"""

import os
import subprocess
import time
import gzip
import sys

LOG = "B03_project/logs/phase2_watcher.log"
OUT = "B03_project/data/proseg_out"
POLL = 60.0

def w(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)

def proseg_alive():
    out = subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         "(Get-Process proseg -ErrorAction SilentlyContinue) -ne $null"],
        capture_output=True, text=True, timeout=60).stdout.strip()
    return out == "True"

def counts_complete():
    """Header nnz == body lines (validates the mtx fully flushed to disk)."""
    try:
        with gzip.open(f"{OUT}/counts", "rt") as f:
            hdr = [next(f) for _ in range(2)]
            nnz = int(hdr[1].split()[2])
            n = 0
            for _ in f:
                n += 1
        return n == nnz, f"nnz {n}/{nnz}"
    except Exception as e:  # noqa: BLE001
        return False, f"{type(e).__name__}: {e}"

w("watcher started: poll 60s for proseg exit")
while True:
    if not proseg_alive():
        w("proseg process gone - checking completion")
        files_ok = all(os.path.exists(f"{OUT}/{f}") for f in
                       ("counts", "cells.parquet", "genes.parquet"))
        mtx_ok, detail = counts_complete() if files_ok else (False, "files missing")
        w(f"outputs: files_complete={files_ok} mtx_complete={mtx_ok} ({detail})")
        if files_ok and mtx_ok:
            w("COMPLETION VERIFIED - launching analysis")
            r = subprocess.run(["python", "B03_project/src/b03_phase2_analysis.py"],
                               capture_output=True, text=True, timeout=3600)
            tail = "\n".join((r.stdout + r.stderr).splitlines()[-12:])
            w(f"analysis exit={r.returncode}; tail:\n{tail}")
            w("WATCHER DONE (analysis fired)")
        else:
            w("CRASH/ABORT detected - analysis NOT run (outputs unusable)")
            w("WATCHER DONE (crash path)")
        sys.exit(0)
    time.sleep(POLL)
