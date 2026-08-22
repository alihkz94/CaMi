#!/usr/bin/env bash
# =============================================================================
# 24_run_stats_isolated.sbatch — run --step STATS beside a live PROFILE run
# -----------------------------------------------------------------------------
# PURPOSE:
#   Run the downstream statistics in their OWN Nextflow session, so they can be
#   run (and re-run) while the 563-sample PROFILE run is still going.
#
#   Once PROFILE has finished, prefer the normal path instead:
#       ./run_pipeline.sh --step STATS
#   which merges everything into one resume cache.
#
# WHY A SEPARATE SESSION:
#   run_pipeline.sh deliberately refuses to start a second Nextflow head while
#   one is already running (job name "cami_nf"). That guard protects the
#   shared resume cache and must not be worked around casually.
#
#   Isolation here comes from two things, both real:
#     1. A separate LAUNCH DIRECTORY, so `.nextflow/cache` and `.nextflow/history`
#        are this session's own. (NXF_HOME does NOT do this — verified 2026-08-15;
#        the cache lives in the launch directory, not in NXF_HOME.)
#     2. A separate work directory (-w), so no task files can ever collide.
#
#   The launch directory must be under /scratch: lib/helpers.nf checkLocation()
#   refuses to run from the /mnt path, because tasks on slurm-406 resolve /scratch
#   locally and would otherwise be handed paths that do not exist for them.
#
# SAFE WHILE PROFILE RUNS:
#   The statistics only READ what PROFILE has published. The integrity gate in
#   aggregate_counts.py quarantines any sample whose files are still half-written,
#   so a partially-finished cohort produces a smaller correct table, never a
#   corrupt one.
#
# RESOURCE FOOTPRINT: tiny. 2 CPUs for the head, and the tasks are 1-2 CPUs on
# tables rather than reads. It does not compete with host removal.
#
# USE:
#   sbatch scripts/analysis/24_run_stats_isolated.sbatch
#   sbatch scripts/analysis/24_run_stats_isolated.sbatch --stats_max_tier 4
#
# VERSION: 1.0 (2026-08-20)
# =============================================================================
#SBATCH --job-name=cockle_stats
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/stats_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/stats_%j.out

set -euo pipefail

DATA="/scratch/ahakimzadeh/Alicia_nature_data"
# CODE is derived from this script's own location, not hardcoded. The repository
# was renamed to CaMi and moved out of the data root; a script that assumes either
# would break silently the next time something moves.
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LAUNCH="$DATA/.stats_session"          # this session's own .nextflow/
WORKDIR="$DATA/work_stats"             # this session's own task files

export TMPDIR="$DATA/tmp"
case "$TMPDIR" in
    /scratch/*) ;;
    *) echo "ERROR: TMPDIR must be under /scratch. Got: $TMPDIR" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR" "$DATA/slurm_logs" "$LAUNCH" "$WORKDIR"

eval "$(/home/ahakimzadeh/miniforge3/bin/conda shell.bash hook)"
conda activate cancer

export NXF_OPTS='-Xms512m -Xmx2g'

# The statistics venv must exist. Build it once with setup_stats_env.sh.
STATS_PY="$DATA/envs/stats/bin/python"
if [[ ! -x "$STATS_PY" ]]; then
    echo "ERROR: statistics environment missing at $STATS_PY" >&2
    echo "       Create it with: bash $CODE/scripts/analysis/setup_stats_env.sh" >&2
    exit 1
fi

cd "$LAUNCH"

echo "=============================================="
echo " node      : $(hostname)"
echo " job       : ${SLURM_JOB_ID}"
echo " launch    : ${LAUNCH}   (isolated .nextflow/)"
echo " workDir   : ${WORKDIR}  (isolated task files)"
echo " stats env : ${STATS_PY}"
echo " arguments : --step STATS $*"
echo "=============================================="

# Reports go to their own directory. nextflow.config points report/timeline/
# trace/dag at results/pipeline_info/ for every session, and overwrite=true, so
# without this a second session would clobber the live run's diagnostics. No
# result data is at risk either way — only the trace files.
INFO="${DATA}/results/pipeline_info_stats"
mkdir -p "$INFO"

# -profile mpi_bremen is REQUIRED here. The site
# settings — the node pinning, the conda environment on PATH, the /scratch guard
# and the database paths — used to be the unconditional defaults in
# nextflow.config. They now live in conf/site/mpi_bremen.config so that the
# pipeline is usable elsewhere. Without the profile this runs with portable
# defaults and every task lands on the wrong node with no software.
exec nextflow run "$CODE/main.nf" -profile mpi_bremen -resume -w "$WORKDIR" --step STATS \
    -with-report   "${INFO}/report.html" \
    -with-timeline "${INFO}/timeline.html" \
    -with-trace    "${INFO}/trace.txt" \
    "$@"
