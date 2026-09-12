"""Proseg run guard v2: exception-proof resource logger + hard abort thresholds.

Defect fixed from v1: a poll failure (e.g. powershell spawn timeout under memory
starvation) raised an uncaught exception and the guard died silently because output
was redirected to /dev/null. v2 never crashes on a poll: every exception is caught,
logged to guard.log, and the loop continues. Only two exits exist: proseg gone
(normal), or an abort threshold tripped (logged loudly).

Triggers (documented before any abort decision):
  - free virtual memory < 3.0 GB (commit exhaustion danger), OR
  - proseg private memory > 14.0 GB (would leave <1.5 GB commit headroom)
"""

import subprocess
import time
import sys

LOG = "B03_project/logs/guard.log"
POLL = 60.0
FREE_VIRT_ABORT_GB = 3.0
PROSEG_PRIV_ABORT_GB = 14.0
PS_CMD = (
    '$os=Get-CimInstance Win32_OperatingSystem;'
    '$p=Get-Process proseg -ErrorAction SilentlyContinue;'
    '"{0:N3}|{1:N3}|{2:N3}|{3:N0}" -f ($os.FreePhysicalMemory/1MB),'
    '($os.FreeVirtualMemory/1MB),'
    '($(if($p){$p.PrivateMemorySize64/1GB}else{-1})),'
    '($(if($p){$p.TotalProcessorTime.TotalSeconds}else{-1}))'
)

def w(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)

def get_ps():
    """Return measurement dict, or None on any failure (never raises)."""
    for attempt in range(3):
        try:
            out = subprocess.run(["powershell", "-NoProfile", "-Command", PS_CMD],
                                 capture_output=True, text=True, timeout=90)
            parts = [p.replace(",", "") for p in out.stdout.strip().split("|")]
            if len(parts) != 4:
                w(f"poll parse failure (attempt {attempt+1}): {out.stdout!r} "
                  f"{out.stderr[:200]!r}")
            else:
                return dict(free_ram=float(parts[0]), free_virt=float(parts[1]),
                            priv_gb=float(parts[2]), cpu_s=float(parts[3]))
        except Exception as e:  # noqa: BLE001 - guard must never die on a poll
            w(f"poll exception (attempt {attempt+1}): {type(e).__name__}: {e}")
        time.sleep(5)
    return None

w("guard v2 started: poll 60s; ABORT if free_virtual<3.0GB or proseg_private>14.0GB")
prev_cpu = None
while True:
    s = get_ps()
    if s is None:
        w("WARNING: 3 poll attempts failed; continuing (no abort decision without data)")
        time.sleep(POLL)
        continue
    if s["priv_gb"] < 0:
        w("proseg process GONE (finished or crashed) - guard exiting")
        sys.exit(0)
    rate = ""
    if prev_cpu is not None:
        rate = f" cpu_rate={(s['cpu_s'] - prev_cpu) / POLL:.2f}x"
    prev_cpu = s["cpu_s"]
    w(f"free_ram={s['free_ram']:.2f}GB free_virt={s['free_virt']:.2f}GB "
      f"proseg_priv={s['priv_gb']:.2f}GB cpu_s={s['cpu_s']:.0f}{rate}")
    if s["free_virt"] < FREE_VIRT_ABORT_GB or s["priv_gb"] > PROSEG_PRIV_ABORT_GB:
        w(f"ABORT TRIGGER: free_virt={s['free_virt']:.2f}GB "
          f"priv={s['priv_gb']:.2f}GB - killing proseg to protect the machine")
        # dynamic PID resolution (v3 fix: v2 hardcoded a stale PID)
        pid_out = subprocess.run(
            ["powershell", "-NoProfile", "-Command",
             "(Get-Process proseg -ErrorAction SilentlyContinue).Id"],
            capture_output=True, text=True, timeout=60).stdout.strip()
        if pid_out:
            subprocess.run(["taskkill", "/PID", pid_out, "/F"], capture_output=True)
            w(f"proseg killed (PID {pid_out}); guard exiting with ABORT status")
        else:
            w("abort triggered but proseg already gone; guard exiting")
        sys.exit(2)
    time.sleep(POLL)
