"""B03 Phase 5 watcher: fires the arrival audit + extraction the moment the
download completes. Polls every 30 s; completion = curl process gone AND zip size
matches the Content-Length (22,618,312,339 B). Then: CRC verify -> G1/G2 audit ->
extract needed members to xenium_lung/extracted/. Stops there: the sender/receiver
map addendum freeze is the next human-gated step (prereg discipline), then compute.

Runs once, exactly once (exit after the chain). All output logged to
xenium_lung/watcher.log. Robust to transient stat errors (download in flight).
"""
import os
import subprocess
import sys
import time

sys.path.insert(0, "B03_project/src")

ZIP = "xenium_lung/Xenium_Prime_Human_Lung_Cancer_FFPE_outs.zip"
EXPECTED = 22618312339
LOG = "xenium_lung/watcher.log"
AUDIT_OUT = "B03_project/logs/lung_audit.json"
EXTRACT_DIR = "xenium_lung/extracted"
EXTRACT_MEMBERS = [
    "gene_panel.json",
    "metrics_summary.csv",
    "cells.parquet",
    "cell_feature_matrix.h5",
    "cell_boundaries.parquet",
    "nucleus_boundaries.parquet",
    "transcripts.csv.gz",
    "transcripts.parquet",
]


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def curl_alive():
    r = subprocess.run(["tasklist", "/FI", "IMAGENAME eq curl.exe"],
                       capture_output=True, text=True)
    return "curl.exe" in r.stdout


def size():
    try:
        return os.path.getsize(ZIP)
    except OSError:
        return 0


def main():
    log(f"watcher armed: waiting for {ZIP} to reach {EXPECTED} B")
    while True:
        s = size()
        alive = curl_alive()
        if s >= EXPECTED and not alive:
            log(f"download complete: {s} B")
            break
        if not alive and s < EXPECTED:
            log(f"WARNING: curl exited early at {s} B - restarting download "
                f"(resume-capable)")
            subprocess.Popen(
                ["curl", "-sS", "-L", "--retry", "8", "--retry-delay", "5", "-C", "-",
                 "-o", ZIP,
                 "https://cf.10xgenomics.com/samples/xenium/3.0.0/"
                 "Xenium_Prime_Human_Lung_Cancer_FFPE/"
                 "Xenium_Prime_Human_Lung_Cancer_FFPE_outs.zip"],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        time.sleep(30)

    # integrity + audit (audit includes full CRC walk)
    log("running arrival audit (CRC walk over 22.6 GB - several minutes)...")
    import b03_lung_audit as A
    audit = A.main(ZIP, AUDIT_OUT)
    log(f"audit complete: {len(audit['testable_pairs'])}/25 pairs testable")

    # extract needed members
    import zipfile
    os.makedirs(EXTRACT_DIR, exist_ok=True)
    zf = zipfile.ZipFile(ZIP)
    for name, member in audit["needed"].items():
        if not member:
            continue
        target = os.path.join(EXTRACT_DIR, name)
        if name == "transcripts.csv.gz" and audit["needed"].get("transcripts.parquet"):
            continue  # parquet preferred when present
        with zf.open(member) as src, open(target, "wb") as dst:
            while True:
                chunk = src.read(1 << 24)
                if not chunk:
                    break
                dst.write(chunk)
        log(f"extracted {name} -> {target} ({os.path.getsize(target):,} B)")
    log("WATCHER CHAIN COMPLETE - pipeline paused at the map-freeze gate")


if __name__ == "__main__":
    main()
