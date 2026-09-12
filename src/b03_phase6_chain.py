"""B03 Phase 6 chain: lung B=10k -> breast B=10k -> determinism gate (B03_PHASE6_PREREG.md).

Execution model matches the sealed Phase 1/5 runs: each wrapper is invoked as a
subprocess (python b03_phase6_<tissue>.py), so the wrapper's __main__ guard fires
exactly once and Windows-spawn workers re-import it as __mp_main__ without
re-executing the run. Importing the wrappers directly would be silently a no-op
(guard never true on import) or, if unguarded, would recurse the entire run on
every spawned worker (the Phase 5 E3 failure class).

The determinism gate is the registered pre-flight: per-permutation Tlo/Thi for
perms 0..999 must match the sealed B=1000 checkpoints EXACTLY (max |delta| = 0.0)
on both tissues before the phase is reported. Any mismatch exits 2 and quarantines
the outputs (registered anomaly path: code audit, no output used). Sequential
execution keeps the memory envelope identical to the sealed runs.
"""
import json
import os
import subprocess
import sys
import time

t0 = time.time()

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))  # Research/ (sealed runs resolve
os.chdir(ROOT)  #  all paths relative to this cwd; src/ is two levels down)
PY = sys.executable

for tag in ("lung", "breast"):
    print(f"[{time.time()-t0:.0f}s] PHASE 6 {tag.upper()}: B=10000 launching subprocess", flush=True)
    r = subprocess.run([PY, os.path.join(HERE, f"b03_phase6_{tag}.py")])
    if r.returncode != 0:
        print(f"PHASE 6 {tag} FAILED rc={r.returncode}; aborting chain before gate.", flush=True)
        sys.exit(1)
    print(f"[{time.time()-t0:.0f}s] PHASE 6 {tag.upper()} null complete", flush=True)

# ------------------------------------------------- determinism gate ------
import numpy as np  # noqa: E402

GATE = {}

for tag, new_ckpt, sealed_ckpt, n_perm_sealed in [
    ("breast", "B03_project/results/b10k_breast/null_v3_checkpoint.npz",
     "B03_project/results/null_v3_checkpoint.npz", 1000),
    ("lung", "xenium_lung/crop/data/b10k/null_v3_checkpoint.npz",
     "xenium_lung/crop/data/null_v3_checkpoint.npz", 1000),
]:
    new = np.load(new_ckpt)
    sealed = np.load(sealed_ckpt)
    done_new, done_sealed = new["done"], sealed["done"]
    assert len(done_new) >= n_perm_sealed, f"{tag}: fewer perms than sealed run"
    assert len(done_sealed) == n_perm_sealed, f"{tag}: sealed ckpt unexpected size"
    # Both checkpoints store rows in completion order with done[] = perm ids.
    # Sort both by perm id and compare the first 1000.
    pos_new = np.argsort(done_new)
    pos_sealed = np.argsort(done_sealed)
    Tlo_new = new["Tlo"][pos_new][:n_perm_sealed]
    Thi_new = new["Thi"][pos_new][:n_perm_sealed]
    Tlo_sealed = sealed["Tlo"][pos_sealed][:n_perm_sealed]
    Thi_sealed = sealed["Thi"][pos_sealed][:n_perm_sealed]
    dlo = float(np.max(np.abs(Tlo_new - Tlo_sealed)))
    dhi = float(np.max(np.abs(Thi_new - Thi_sealed)))
    GATE[tag] = {
        "shared_perms": n_perm_sealed,
        "max_abs_delta_Tlo": dlo,
        "max_abs_delta_Thi": dhi,
        "pass": bool(dlo == 0.0 and dhi == 0.0),
    }
    print(f"[{time.time()-t0:.0f}s] GATE {tag}: shared={n_perm_sealed} "
          f"max|dTlo|={dlo:.3g} max|dThi|={dhi:.3g} -> "
          f"{'PASS' if GATE[tag]['pass'] else 'FAIL'}", flush=True)

overall = all(g["pass"] for g in GATE.values())
report = {
    "phase": "B03 Phase 6 (B=10000 null resolution)",
    "prereg": "B03_PHASE6_PREREG.md",
    "determinism_gate": GATE,
    "overall": "PASS" if overall else "FAIL",
    "wall_seconds": time.time() - t0,
}
os.makedirs("B03_project/results/phase6", exist_ok=True)
with open("B03_project/results/phase6/determinism_gate.json", "w") as fh:
    json.dump(report, fh, indent=1)

if not overall:
    print("PHASE 6 GATE FAILED: outputs quarantined, code audit required "
          "(registered anomaly path).", flush=True)
    sys.exit(2)

print(f"[{time.time()-t0:.0f}s] PHASE 6 COMPLETE: gate PASS on both tissues. "
      f"b10k null CSVs are in B03_project/results/b10k_breast and "
      f"xenium_lung/crop/data/b10k.", flush=True)
