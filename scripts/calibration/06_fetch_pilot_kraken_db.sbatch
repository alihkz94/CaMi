#!/usr/bin/env bash
# =============================================================================
# 06_fetch_pilot_kraken_db.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Fetch k2_standard_16gb_20250402 — the EXACT Kraken2 database the pilot used.
#   We need it to calibrate, not to analyse: the shared PlusPF database is far
#   larger, and with --confidence 0.1 it reported 26x fewer bacteria than the
#   pilot on the same run (ERR10680552). Reproducing the pilot baseline tells us
#   whether that gap is the database, the threshold, or something else.
#
#   Note: do NOT run this with `sbatch --wrap`. --wrap executes the body with
#   /bin/sh (dash), where `set -o pipefail` is an illegal option and the job dies
#   immediately. A real script file with a bash shebang is used instead.
#
# OUT: kraken2_db/standard_16gb_20250402/{hash,opts,taxo}.k2d + kmer_distrib files
#
# RUN: sbatch 06_fetch_pilot_kraken_db.sbatch
# VERSION: 1.0 (2026-08-14)
# =============================================================================
#SBATCH --job-name=k2_16gb
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/k2_16gb_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/k2_16gb_%j.out

set -euo pipefail

export PATH=/home/ahakimzadeh/miniforge3/envs/cancer/bin:$PATH
export TMPDIR=/scratch/ahakimzadeh/Alicia_nature_data/tmp
case "$TMPDIR" in
    /scratch/*) ;;
    *) echo "ERROR: TMPDIR must be under /scratch. Got: $TMPDIR" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"

DB_DIR="/scratch/ahakimzadeh/Alicia_nature_data/kraken2_db/standard_16gb_20250402"
URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20250402.tar.gz"

mkdir -p "$DB_DIR"
cd "$DB_DIR"

if [[ -f hash.k2d && -f taxo.k2d && -f opts.k2d ]]; then
    echo "[$(date '+%F %T')] database already present — nothing to do"
    exit 0
fi

echo "[$(date '+%F %T')] downloading $URL"
# The site proxy flaps (README trap 12). aria2c's own retries plus an outer wait.
for attempt in 1 2 3 4; do
    if aria2c -x8 -s8 -c --console-log-level=warn --retry-wait=15 --max-tries=8 \
              -o db.tar.gz "$URL"; then
        break
    fi
    if (( attempt == 4 )); then
        echo "ERROR: download failed after 4 attempts" >&2
        exit 1
    fi
    wait_s=$(( attempt * 300 ))
    echo "[$(date '+%F %T')] attempt $attempt failed; waiting ${wait_s}s for the proxy"
    sleep "$wait_s"
done

echo "[$(date '+%F %T')] extracting"
tar -xzf db.tar.gz
rm -f db.tar.gz

echo "[$(date '+%F %T')] DONE"
ls -la
n=$(ls -1 ./*.k2d 2>/dev/null | wc -l)
echo "k2d files: $n"
(( n >= 3 )) || { echo "ERROR: expected 3 .k2d files, found $n" >&2; exit 1; }
