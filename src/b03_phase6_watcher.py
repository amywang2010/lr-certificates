"""B03 Phase 6 finish-line watcher (automation of the registered plan only).

The chain's registered job ends at the determinism gate: it runs both B=10k
nulls, compares per-permutation Tlo/Thi bitwise against the sealed B=1k
checkpoints, writes determinism_gate.json, and exits. B03_PHASE6_PREREG.md
then requires the verdict-stability analysis (b03_phase6_analysis.py) on gate
PASS. Until now that step was manual; this watcher closes the gap so an
~04:00 completion is not idle until morning.

Behavior:
  - Poll every 60 s for determinism_gate.json.
  - overall == "PASS" -> run the analysis ONCE (atomic marker claim so a
    watcher restart cannot double-run), then verify the two registered
    analysis outputs (b10k/analysis/b10k_null.csv under each tissue's
    out_dir) exist with mtime after the fire time.
  - overall != "PASS" -> touch nothing and exit 2 (registered anomaly path:
    outputs quarantined, human code audit required).
  - Idempotent: only if BOTH registered per-tissue outputs already exist.
    (DEV-031: the original any-match check was satisfied by lung's Sep-8
    output alone, causing a false "already present" skip of the breast
    analysis. Presence is now verified per tissue.)
  - A gate file that fails to parse (caught mid-write) is retried, not fatal.
"""
import json
import os
import subprocess
import sys
import time

GATE = "B03_project/results/phase6/determinism_gate.json"
MARK = "B03_project/results/phase6/analysis_fired.json"
SCRIPT = "B03_project/src/b03_phase6_analysis.py"

# Registered analysis outputs, per tissue (must ALL exist post-fire).
EXPECTED_OUTPUTS = [
    "xenium_lung/crop/data/b10k/analysis/b10k_null.csv",
    "B03_project/results/b10k_breast/analysis/b10k_null.csv",
]


def find_outputs(after_ts):
    """Registered outputs with mtime strictly after the fire time."""
    return [p for p in EXPECTED_OUTPUTS
            if os.path.exists(p) and os.path.getmtime(p) > after_ts]


def all_outputs_present():
    return all(os.path.exists(p) for p in EXPECTED_OUTPUTS)


def main():
    while not os.path.exists(GATE):
        time.sleep(60)
    try:
        g = json.load(open(GATE))
    except json.JSONDecodeError:
        time.sleep(5)
        g = json.load(open(GATE))

    if g.get("overall") != "PASS":
        print(f"GATE overall={g.get('overall')!r}; NOT running analysis "
              f"(registered anomaly path: quarantine + human audit)",
              flush=True)
        return 2
    if all_outputs_present():
        print("analysis outputs already present; nothing to do", flush=True)
        return 0
    if os.path.exists(MARK):
        print("marker exists but no outputs found; a previous fire failed - "
              "inspect analysis_fired.json", flush=True)
        return 1

    fired = time.time()
    tmp = MARK + ".tmp"
    fd = os.open(tmp, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, json.dumps({"fired_at": fired}).encode())
    os.close(fd)
    os.replace(tmp, MARK)

    print(f"gate PASS; firing registered analysis at "
          f"{time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    r = subprocess.run([sys.executable, SCRIPT],
                       capture_output=True, text=True)
    record = {
        "fired_at": fired,
        "rc": r.returncode,
        "wall_s": time.time() - fired,
        "stdout_tail": r.stdout[-2000:],
        "stderr_tail": r.stderr[-2000:],
    }
    with open(MARK, "w") as fh:
        json.dump(record, fh, indent=1)

    outs = find_outputs(fired)
    ok = r.returncode == 0 and len(outs) >= len(EXPECTED_OUTPUTS)
    print(f"{'ANALYSIS COMPLETE' if ok else 'ANALYSIS FAILED'} "
          f"rc={r.returncode} wall={record['wall_s']:.0f}s "
          f"fresh_outputs={outs} (of {len(EXPECTED_OUTPUTS)} required)", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
