#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — the MPI Bremen launcher for CaMi
# -----------------------------------------------------------------------------
# THIS IS A SITE WRAPPER, NOT THE PIPELINE INTERFACE.
#
# Everywhere else, run CaMi the ordinary way and skip this file entirely:
#     nextflow run alihkz94/CaMi -profile singularity,slurm \
#         --dataRoot /path/to/cohort --step PROFILE -resume
#
# This wrapper exists because one cluster needs three things that no command line
# should have to remember, and getting any of them wrong costs days:
#
#   1. The head process must run ON slurm-406. The data is on that node's local
#      disk, and the same bytes have a different path from the login node
#      (/scratch/... vs /mnt/slurm-406/...). A head started on the wrong side
#      hands every task a path that does not exist there — and, because -resume
#      keys on the path, would also split the cache and repeat finished work.
#
#   2. TMPDIR must be on /scratch. Filling /tmp on a shared node breaks the node
#      for every user; the administrators have asked that this not happen again.
#
#   3. Only one Nextflow head at a time. Two runs in this directory share
#      .nextflow/cache and the history file, and running them together corrupts
#      the resume cache.
#
# It always adds -resume, so a stopped or failed run continues where it stopped.
#
# USE:
#   ./run_pipeline.sh --step REFERENCE
#   ./run_pipeline.sh --step FETCH_BAM
#   ./run_pipeline.sh --step PROFILE --samples PACE17_428F,PACE17_456F --subsample 200000
#   ./run_pipeline.sh --step STATS
#   ./run_pipeline.sh                      # REFERENCE, FETCH_BAM and PROFILE, in order
#
# WATCH:  tail -f ../slurm_logs/nf_<jobid>.out
# =============================================================================
#SBATCH --job-name=cami_nf
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# 14 days, not 7. MEASURED 2026-08-16: host removal averages ~2 h per sample at
# 16 CPUs, and the node fits ~7 concurrently, so 498 samples need ~5-6 days for
# that stage alone. A 7-day head job left about 13 hours of margin — and the run
# submitted with that limit did in fact reach it. The partition allows 120 days
# (MaxTime), so 7 was a self-imposed limit, not a cluster one.
#
# A USER CANNOT RAISE THE LIMIT ON A RUNNING JOB (scontrol denies it). This has
# to be right at submission time, or the run dies at the wall clock with no
# recourse beyond resubmitting with -resume.
#SBATCH --time=14-00:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/CaMi
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/nf_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/nf_%j.out

set -euo pipefail

# CODE root (this repository) and DATA root. They are separate on purpose: the
# repository is a few megabytes of text, the data is about 4 TB.
CODE="/scratch/ahakimzadeh/CaMi"
DATA="/scratch/ahakimzadeh/Alicia_nature_data"

# Only one Nextflow head at a time — see reason 3 in the header.
#
# No grep in a pipeline here: with `set -e` and `pipefail`, a grep that matches
# nothing exits 1, the command substitution below fails, and the script dies
# silently in the ordinary case where no other job is running.
running_heads() {
    local ids out=""
    ids=$(squeue -u "$USER" -h -n cami_nf -o '%A' 2>/dev/null || true)
    for id in $ids; do
        [[ "$id" == "${SLURM_JOB_ID:-none}" ]] && continue
        out+="$id "
    done
    printf '%s' "$out"
}

OTHER="$(running_heads)"
if [[ -n "${OTHER// /}" ]]; then
    echo "ERROR: another Nextflow head job is already running: ${OTHER}" >&2
    echo "  Two runs would share .nextflow/cache and break -resume." >&2
    echo "  Wait for it, or stop it with: scancel ${OTHER}" >&2
    exit 1
fi

# Not in a Slurm job yet: submit this same file and stop.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submitting the Nextflow head job to slurm-406..."
    exec sbatch "$0" "$@"
fi

# --- /tmp guard — see reason 2 in the header ---------------------------------
export TMPDIR="$DATA/tmp"
case "$TMPDIR" in
    /scratch/*) ;;
    *) echo "ERROR: TMPDIR must be under /scratch. Got: $TMPDIR" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR" "$DATA/slurm_logs" "$DATA/results"

# --- tools -------------------------------------------------------------------
# This cluster has no Singularity, Apptainer or Docker, so software comes from a
# conda environment. conf/site/mpi_bremen.config puts the same directory on PATH
# inside every task, because a Slurm child job does not inherit this activation.
eval "$(/home/ahakimzadeh/miniforge3/bin/conda shell.bash hook)"
conda activate cancer

# Keep Nextflow's own state inside the project, on /scratch.
export NXF_HOME="$CODE/.nextflow"
export NXF_WORK="$DATA/work"
export NXF_OPTS='-Xms1g -Xmx4g'

cd "$CODE"

echo "=============================================="
echo " pipeline  : CaMi"
echo " node      : $(hostname)"
echo " job       : ${SLURM_JOB_ID}"
echo " code      : ${CODE}"
echo " data      : ${DATA}"
echo " work dir  : ${NXF_WORK}"
echo " TMPDIR    : ${TMPDIR}"
echo " arguments : $*"
echo "=============================================="

exec nextflow run main.nf -profile mpi_bremen -resume "$@"
