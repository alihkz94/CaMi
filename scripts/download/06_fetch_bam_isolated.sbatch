#!/usr/bin/env bash
# =============================================================================
# 06_fetch_bam_isolated.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Run --step FETCH_BAM in its OWN Nextflow session, separate from the main
#   PROFILE session. Use this ONLY while a PROFILE run is already active on
#   slurm-406. Once PROFILE finishes, use run_pipeline.sh normally instead —
#   that is the standard path and it merges everything into one resume cache.
#
#   WHY A SEPARATE SESSION:
#     run_pipeline.sh refuses to start a second Nextflow head while one is
#     already running (job name "cami_nf"), on purpose.
#
#   WHAT ACTUALLY PROVIDES THE ISOLATION — verified 2026-08-15, and NOT what an
#   earlier version of this comment claimed:
#     Nextflow keeps its resume cache and history in `.nextflow/` inside the
#     LAUNCH DIRECTORY, not inside NXF_HOME. Setting NXF_HOME here only moves
#     plugins/secrets/tmp; `.nextflow_bamfetch/` holds no cache at all. Both this
#     job and the main PROFILE job launch from the same directory and therefore
#     SHARE `.nextflow/cache/` and `.nextflow/history`.
#
#     They stay safe for two real reasons:
#       1. Different SESSION IDs -> different `.nextflow/cache/<uuid>/` subdirs,
#          so neither can corrupt the other's resume state.
#       2. Separate `-w` work directories, so no task files ever collide.
#
#     The one genuinely shared file is the append-only `.nextflow/history` log.
#     Concurrent appends there are metadata only and cannot damage results.
#
#   DO NOT run `nextflow clean` while a PROFILE run is live. `-resume` keeps one
#   session ID across runs, so the old --subsample test runs and the production
#   run share a session: cleaning "the old test session" deletes the LIVE run's
#   work directory.
#
#   RESOURCE FOOTPRINT: small on purpose. FETCH_UNMAPPED_BAM asks for 4 CPUs,
#   runs one sample at a time (maxForks 1), and does mostly network I/O. This
#   does not compete meaningfully with the PROFILE run for the node's 128 CPUs.
#
# USE:
#   sbatch 06_fetch_bam_isolated.sbatch --samples EECE18_949H-wga,EPCE18_851H
#   sbatch 06_fetch_bam_isolated.sbatch                # all 65, one at a time
#
# VERSION: 1.0 (2026-08-15)
# =============================================================================
#SBATCH --job-name=cockle_bam
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
# 10 days, not 3. Each fetch is only ~6 minutes of real work, but while a PROFILE
# run holds the node these tasks wait hours for a free CPU slot: 8 of 65 completed
# in the first 20 hours. The limit must cover the WAITING, not the work.
#SBATCH --time=10-00:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/CaMi
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/bamfetch_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/bamfetch_%j.out

set -euo pipefail

DATA="/scratch/ahakimzadeh/Alicia_nature_data"
# CODE is derived from this script's own location, not hardcoded. The repository
# was renamed to CaMi and moved out of the data root; a script that assumes either
# would break silently the next time something moves.
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

export TMPDIR="$DATA/tmp"
case "$TMPDIR" in
    /scratch/*) ;;
    *) echo "ERROR: TMPDIR must be under /scratch. Got: $TMPDIR" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR" "$DATA/slurm_logs"

eval "$(/home/ahakimzadeh/miniforge3/bin/conda shell.bash hook)"
conda activate cancer

# Isolation: distinct NXF_HOME (resume cache + history) AND distinct work
# directory. Neither is shared with the live PROFILE session.
export NXF_HOME="$CODE/.nextflow_bamfetch"
export NXF_OPTS='-Xms512m -Xmx2g'

cd "$CODE"

echo "=============================================="
echo " node      : $(hostname)"
echo " job       : ${SLURM_JOB_ID}"
echo " NXF_HOME  : ${NXF_HOME}   (isolated from the main run)"
echo " workDir   : ${DATA}/work_bamfetch  (isolated from the main run)"
echo " arguments : --step FETCH_BAM $*"
echo "=============================================="

# nextflow.config points report/timeline/trace/dag at the SAME results/pipeline_info/
# regardless of which session runs — NXF_HOME and workDir isolate the resume cache
# and the task work directories, but not this. A second session starting here
# overwrites (trace.overwrite=true) whatever the main PROFILE session had already
# written. It does not touch any real result — every per-sample file is still
# correct and independent — only the diagnostic trace/report get clobbered.
# Point this session's reports somewhere of its own.
BAMINFO="${DATA}/results/pipeline_info_bamfetch"
mkdir -p "$BAMINFO"

# -profile mpi_bremen is REQUIRED here. The site
# settings — the node pinning, the conda environment on PATH, the /scratch guard
# and the database paths — used to be the unconditional defaults in
# nextflow.config. They now live in conf/site/mpi_bremen.config so that the
# pipeline is usable elsewhere. Without the profile this runs with portable
# defaults and every task lands on the wrong node with no software.
exec nextflow run "$CODE/main.nf" -profile mpi_bremen -resume -w "${DATA}/work_bamfetch" --step FETCH_BAM \
    -with-report   "${BAMINFO}/report.html" \
    -with-timeline "${BAMINFO}/timeline.html" \
    -with-trace    "${BAMINFO}/trace.txt" \
    "$@"
