"""B03 Phase 6 continuous monitor (non-invasive; reads only).

Appends one JSON line every 300 s to results/phase6/monitor.log (repository-relative):
  ts, done (perms), dl_bytes, battery_status, charge_pct, free_phys_mb,
  and a "stall" flag when the null's done count has not advanced for
  >= 15 min (900 s) — the slip-catcher the polling discipline needs.
Exit codes: runs until killed. Reads checkpoint read-only; np.load does not
lock the writer's rename cycle (checkpoint writes are os.replace-atomic in v3).
"""
import datetime, os
import json
import subprocess
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CKPT = os.path.join(ROOT, "results", "b10k_breast", "null_v3_checkpoint.npz")
WORKSPACE = os.environ.get("B03_WORKSPACE", os.path.dirname(ROOT))
DL = os.path.join(WORKSPACE, "xenium_breast_s6", "sdata_breast_s6.zarr.zip")
LOG = os.path.join(ROOT, "results", "phase6", "monitor.log")
STALL_S = 900.0
POLL_S = 300.0


def read_done():
    try:
        return int(np.load(CKPT)["done"].max())
    except Exception:
        return -1


def dl_bytes():
    try:
        import os
        final = os.path.getsize(DL)
        return final
    except Exception:
        pass
    # v6 segmented phase: final file absent; count accepted chunk bytes
    try:
        import os
        cdir = os.path.join(WORKSPACE, "xenium_breast_s6", "chunks")
        return sum(os.path.getsize(os.path.join(cdir, f))
                   for f in os.listdir(cdir) if f.endswith(".bin"))
    except Exception:
        return -1


def power():
    try:
        out = subprocess.run(
            ["powershell", "-NoProfile", "-Command",
             "$b=Get-CimInstance Win32_Battery; \"$($b.BatteryStatus),$($b.EstimatedChargeRemaining)\""],
            capture_output=True, text=True, timeout=30).stdout.strip()
        st, pct = out.split(",")
        return int(st), int(pct)
    except Exception:
        return -1, -1


def free_mb():
    try:
        out = subprocess.run(
            ["wmic", "OS", "get", "FreePhysicalMemory", "/format:list"],
            capture_output=True, text=True, timeout=30).stdout
        for line in out.splitlines():
            if line.startswith("FreePhysicalMemory="):
                return int(line.split("=")[1].strip()) // 1024
        return -1
    except Exception:
        return -1


def main():
    last_done = read_done()
    last_change = time.time()
    while True:
        time.sleep(POLL_S)
        done = read_done()
        now = time.time()
        if done != last_done:
            last_done, last_change = done, now
        st, pct = power()
        rec = {
            "ts": datetime.datetime.now().isoformat(timespec="seconds"),
            "done": done,
            "dl_bytes": dl_bytes(),
            "battery_status": st,
            "charge_pct": pct,
            "free_phys_mb": free_mb(),
            "stall_min": round((now - last_change) / 60.0, 1),
        }
        if now - last_change >= STALL_S:
            rec["stall"] = True
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
