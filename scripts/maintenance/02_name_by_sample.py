#!/usr/bin/env python3
# =============================================================================
# 02_name_by_sample.py
# -----------------------------------------------------------------------------
# PURPOSE:
#   Give every downloaded FASTQ its ORIGINAL sample name (e.g. ENCE17_4531H-wga)
#   instead of the ENA run accession (ERR10673701).
#
#   Default mode is `link`: the canonical files stay as <RUN>_N.fastq.gz in
#   fastq/ and by_sample/ gets symlinks named <SAMPLE_CODE>_N.fastq.gz. That way
#   you get the readable names AND every bundled pipeline script — which all
#   expect <ACC>_1.fastq.gz — keeps working untouched. Symlinks cost no disk.
#
#   `--mode rename` physically renames the files instead. Only use it if you have
#   decided the accession-named layout is not needed, because it breaks
#   profile_microbiome.sh / validate_taxa.sh / the Nextflow inputs.
#
# INPUTS:
#   ../manifests/download_list.tsv   (idx run sample_code status prep gb)
#   ../fastq/<RUN>_{1,2}.fastq.gz
#
# OUTPUTS:
#   ../by_sample/<SAMPLE_CODE>_{1,2}.fastq.gz    (symlinks, or real files if renamed)
#   ../manifests/sample_name_map.csv             full run <-> sample-code mapping
#
# CONDA ENV: cancer  (stdlib only)
#
# RUN:
#   python3 02_name_by_sample.py               # symlink farm (safe, default)
#   python3 02_name_by_sample.py --mode rename # physically rename
#
# VERSION: 1.0  (2026-08-12)
# =============================================================================
import argparse
import csv
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "manifests" / "download_list.tsv"
NAME_MAP = ROOT / "manifests" / "sample_name_map.csv"

# The cohort lives on one node's node-local /scratch, which has two different
# paths depending on where you run this: /scratch on the node itself, and
# /mnt/<node> from the login node. Resolve whichever is visible.
# Override with the COCKLE_DATA environment variable.
#
# NODE is whichever node 01_download_array.sbatch currently targets (--nodelist=).
# It has moved before (slurm-404 -> slurm-406, 2026-08-14, see ../README.md).
# Note /mnt/<node> exists ONLY on the login node; compute nodes see /scratch.
NODE = "slurm-406"
USER = os.environ.get("USER", "ahakimzadeh")
if os.environ.get("COCKLE_DATA"):
    DATA_ROOT = Path(os.environ["COCKLE_DATA"])
elif Path(f"/scratch/{USER}/Alicia_nature_data").is_dir():
    DATA_ROOT = Path(f"/scratch/{USER}/Alicia_nature_data")      # on the node
else:
    DATA_ROOT = Path(f"/mnt/{NODE}/{USER}/Alicia_nature_data")   # on the login node

if not DATA_ROOT.is_dir():
    sys.exit(f"ERROR: {DATA_ROOT} does not exist.\n"
             f"  NODE is hardcoded to '{NODE}' above — if the download moved to a\n"
             f"  different node, update NODE (or set COCKLE_DATA=/path/to/Alicia_nature_data).")

FASTQ_DIR = DATA_ROOT / "fastq"
BY_SAMPLE = DATA_ROOT / "by_sample"


def load_manifest():
    rows = []
    with MANIFEST.open() as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if not r:
                continue
            rows.append({"idx": r[0], "run": r[1], "sample_code": r[2],
                         "status": r[3], "prep": r[4], "gb": r[5]})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["link", "rename"], default="link",
                    help="link = symlinks in by_sample/ (default); "
                         "rename = physically rename files in fastq/")
    args = ap.parse_args()

    rows = load_manifest()
    BY_SAMPLE.mkdir(parents=True, exist_ok=True)

    # Sample codes must be unique or the by_sample/ layout would silently collide.
    codes = [r["sample_code"] for r in rows]
    dupes = {c for c in codes if codes.count(c) > 1}
    if dupes:
        print(f"ERROR: sample codes are not unique: {sorted(dupes)[:5]}", file=sys.stderr)
        return 1

    done = missing = 0
    map_rows = []

    for r in rows:
        run, code = r["run"], r["sample_code"]
        verified = (FASTQ_DIR / f"{run}.md5.ok").exists()
        for mate in ("1", "2"):
            src = FASTQ_DIR / f"{run}_{mate}.fastq.gz"
            dst_name = f"{code}_{mate}.fastq.gz"
            if not src.exists():
                missing += 1
                map_rows.append([run, code, r["status"], r["prep"], r["gb"],
                                 mate, src.name, dst_name,
                                 "verified" if verified else "not_downloaded"])
                continue

            if args.mode == "link":
                dst = BY_SAMPLE / dst_name
                # Relative target keeps the tree movable as a whole.
                target = os.path.relpath(src, BY_SAMPLE)
                if dst.is_symlink() or dst.exists():
                    dst.unlink()
                dst.symlink_to(target)
            else:
                dst = FASTQ_DIR / dst_name
                if dst.exists():
                    dst.unlink()
                src.rename(dst)
            done += 1
            map_rows.append([run, code, r["status"], r["prep"], r["gb"],
                             mate, src.name, dst_name,
                             "verified" if verified else "present_unverified"])

    with NAME_MAP.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_accession", "sample_code", "disease_status", "dna_prep",
                    "fastq_gb", "mate", "accession_filename", "sample_filename", "state"])
        w.writerows(map_rows)

    verb = "symlinked" if args.mode == "link" else "renamed"
    print(f"{verb}: {done} file(s)")
    print(f"not yet downloaded: {missing} file(s)")
    print(f"mapping written to: {NAME_MAP}")
    if args.mode == "link":
        print(f"readable names in : {BY_SAMPLE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
