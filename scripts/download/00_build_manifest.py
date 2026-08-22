#!/usr/bin/env python3
# =============================================================================
# 00_build_manifest.py
# -----------------------------------------------------------------------------
# PURPOSE:
#   Build the download manifest that drives the Slurm array in
#   01_download_array.sbatch. One row per ENA run that actually has FASTQ,
#   carrying the original sample code so downloads can be renamed later.
#
# INPUTS:
#   alicia_bundle/run_metadata.tsv   (563 runs; sample_title == original code)
#   sample_info .csv                 (Sample_code -> ERS accession; cross-check)
#
# OUTPUTS:
#   manifests/download_list.tsv      tab-separated, NO header, array-indexed:
#                                    idx  run  sample_code  status  prep  gb
#   manifests/skipped_bam_only.tsv   runs with no ENA FASTQ (byte-range route)
#
# CONDA ENV: cancer  (stdlib only, any python3 works)
#
# RUN:
#   python3 00_build_manifest.py
#
# VERSION: 1.0  (2026-08-12)
# =============================================================================
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                       # .../Alicia_nature_data — everything is under here

RUN_META = ROOT / "alicia_bundle" / "run_metadata.tsv"
SAMPLE_INFO = ROOT / "sample_info .csv"   # note: the space in the filename is real
MANIFEST_DIR = ROOT / "manifests"


def main() -> int:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)

    with RUN_META.open() as fh:
        runs = list(csv.DictReader(fh, delimiter="\t"))

    # sample_info maps the original sample code to its ERS sample accession.
    # run_metadata already carries the code as `sample_title`; we use sample_info
    # purely as a consistency check so a silent metadata drift can't rename files wrong.
    with SAMPLE_INFO.open() as fh:
        info = {r["Sample_code"]: r["SAMPLE"] for r in csv.DictReader(fh)}

    titles = {r["sample_title"] for r in runs}
    missing = sorted(set(info) - titles)
    if missing:
        print(f"ERROR: {len(missing)} sample_info codes absent from run_metadata: "
              f"{missing[:5]}", file=sys.stderr)
        return 1
    print(f"cross-check OK: all {len(info)} sample_info codes present in run_metadata")

    have, bam_only = [], []
    for r in runs:
        (have if r["fastq_available"] == "true" else bam_only).append(r)

    # Largest first: long-pole runs start early, so the array's tail is short runs
    # and the whole job finishes sooner under a fixed concurrency throttle.
    have.sort(key=lambda r: float(r["fastq_gb"] or 0), reverse=True)

    out = MANIFEST_DIR / "download_list.tsv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for i, r in enumerate(have, start=1):
            w.writerow([i, r["run_accession"], r["sample_title"],
                        r["disease_status"], r["dna_prep"], r["fastq_gb"]])

    skipped = MANIFEST_DIR / "skipped_bam_only.tsv"
    with skipped.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run_accession", "sample_code", "disease_status", "dna_prep",
                    "submitted_format"])
        for r in bam_only:
            w.writerow([r["run_accession"], r["sample_title"], r["disease_status"],
                        r["dna_prep"], r["submitted_format"]])

    total_gb = sum(float(r["fastq_gb"] or 0) for r in have)
    print(f"\nmanifest : {out}")
    print(f"  runs with FASTQ : {len(have)}")
    print(f"  total size      : {total_gb:.1f} GB ({total_gb/1024:.2f} TB)")
    print(f"  array range     : 1-{len(have)}")
    print(f"\nBAM-only (no ENA FASTQ, need golden_set byte-range route): "
          f"{len(bam_only)} -> {skipped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
