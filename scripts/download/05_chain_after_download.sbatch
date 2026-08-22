#!/usr/bin/env bash
# =============================================================================
# 05_chain_after_download.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Start the Nextflow pipeline automatically as soon as the download array has
#   finished, so no time is lost waiting for someone to look at the queue.
#
#   It is submitted with a Slurm dependency on the download job, so it wakes up
#   the moment that job ends (for any reason).
#
#   It does NOT simply launch the pipeline. First it counts what is verified:
#
#     all 498 verified            -> launch the pipeline
#     some missing, rounds left   -> resubmit the download array, re-arm itself
#     some missing, no rounds left-> launch the pipeline with what exists, and say so
#
#   The middle case is the important one. On 2026-08-14 the site web proxy blinked
#   and 337 runs failed together in under a minute. A round of retries recovers
#   from that without anybody being awake for it.
#
# IN : manifests/download_list.tsv , fastq/<ACC>.md5.ok
# OUT: submits 01_download_array.sbatch and/or CaMi/run_pipeline.sh
#
# RUN (arm it against the running download job):
#   cd /mnt/slurm-406/ahakimzadeh/Alicia_nature_data/scripts
#   sbatch --dependency=afterany:<DOWNLOAD_JOBID> 05_chain_after_download.sbatch
#
#   Check it is armed:  squeue -u $USER -n cockle_chain
#   Cancel it:          scancel -n cockle_chain
#
# VERSION: 1.0 (2026-08-14)
# =============================================================================
#SBATCH --job-name=cockle_chain
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/chain_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/chain_%j.out

set -uo pipefail

# Data root, code repository, and the download scripts inside it. These are
# absolute on purpose: sbatch runs this file from a spool copy, so $0 and
# BASH_SOURCE do not point at the repository.
PROJECT="/scratch/ahakimzadeh/Alicia_nature_data"          # data root
# CODE is derived from this script's own location, not hardcoded. The repository
# was renamed to CaMi and moved out of the data root; a script that assumes either
# would break silently the next time something moves.
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPTS="$CODE/scripts/download"
MANIFEST="$PROJECT/manifests/download_list.tsv"

ROUND="${1:-1}"        # which retry round this is
MAX_ROUNDS=4           # after this, stop retrying and process what we have

# Which pipeline step to launch. PROFILE, not 'all', on purpose: 'all' would also
# run FETCH_BAM, and the byte-range BAM fetch is not yet proven through this site
# proxy — on 2026-08-14 the proxy injected HTTP 503 into the ranges and produced a
# full-size BAM with a CRC32 mismatch. Do not start unproven work unattended.
# Run FETCH_BAM by hand once it is verified. Override: sbatch ... 05_...sbatch 1 all
STEP="${2:-PROFILE}"

log() { echo "[$(date '+%F %T')] $*"; }

log "chain round $ROUND/$MAX_ROUNDS starting on $(hostname), will launch --step $STEP"

TOTAL=$(wc -l < "$MANIFEST")
DONE=$(find "$PROJECT/fastq" -maxdepth 1 -name '*.md5.ok' 2>/dev/null | wc -l)
LEFT=$(( TOTAL - DONE ))

log "verified $DONE / $TOTAL runs, $LEFT outstanding"

# --- everything is here: go ---------------------------------------------------
if (( LEFT <= 0 )); then
    log "download complete. Launching the pipeline."
# Launch the pipeline with sbatch, NOT by calling ./run_pipeline.sh directly.
# run_pipeline.sh submits itself only when SLURM_JOB_ID is empty. This chain IS a
# Slurm job, so that variable is already set, and a direct call would run the
# Nextflow head INSIDE this job — 1 CPU, 2 GB, 30 minutes. The 498-sample run
# would be killed half an hour in, overnight, with nobody watching.
    NF_JOB=$(sbatch --parsable "$CODE/run_pipeline.sh" --step "$STEP")
    if [[ -z "$NF_JOB" ]]; then
        log "ERROR: could not submit the pipeline. Launch it by hand:"
        log "  cd $CODE && ./run_pipeline.sh --step $STEP"
        exit 1
    fi
    log "pipeline submitted as job $NF_JOB"
    exit 0
fi

# --- still missing runs, and we have rounds left: retry the download ----------
if (( ROUND < MAX_ROUNDS )); then
    log "$LEFT run(s) still missing. Resubmitting the download array (round $((ROUND+1)))."
    NEW_DL=$(sbatch --parsable "$SCRIPTS/01_download_array.sbatch")
    if [[ -z "$NEW_DL" ]]; then
        log "ERROR: could not resubmit the download array. Not re-arming."
        exit 1
    fi
    log "download resubmitted as job $NEW_DL"

    NEW_CHAIN=$(sbatch --parsable --dependency=afterany:"$NEW_DL" \
                       "$SCRIPTS/05_chain_after_download.sbatch" $((ROUND+1)) "$STEP")
    log "chain re-armed as job $NEW_CHAIN, waiting on $NEW_DL"
    exit 0
fi

# --- out of rounds: process what we have, and be explicit about it ------------
log "WARNING: $LEFT run(s) never completed after $MAX_ROUNDS rounds."
log "Launching the pipeline with the $DONE verified runs."
log "The missing runs are listed by: bash $SCRIPTS/03_status.sh"
log "The pipeline only reads runs that have an .md5.ok marker, so nothing unverified is used."
NF_JOB=$(sbatch --parsable "$CODE/run_pipeline.sh" --step "$STEP")
if [[ -z "$NF_JOB" ]]; then
    log "ERROR: could not submit the pipeline. Launch it by hand:"
    log "  cd $CODE && ./run_pipeline.sh --step $STEP"
    exit 1
fi
log "pipeline submitted as job $NF_JOB"
