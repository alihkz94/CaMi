#!/usr/bin/env bash
# =============================================================================
# 04_migrate_to_marinecancers.sbatch  (v1.0, 2026-08-13)
#
# Move the project to the group directory and delete the originals.
# rsync -> md5 every FASTQ source vs dest -> only then rm the source.
# A single mismatch aborts the delete. Re-runnable (rsync skips what's copied).
#
# Runs as a Slurm job because login sessions may still carry the old group set;
# compute nodes already resolve `marinecancers`.
#
# BLOCKED as of 2026-08-13: /misc writes hang (mkdir ok, file data hangs).
# Check with a small write before submitting this. See ../README.md.
#
# Edit SRC/DEST below if the paths changed.
# =============================================================================
#SBATCH --job-name=cockle_migrate
#SBATCH --partition=slurm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --nodelist=slurm-406          # data lives on this node's /scratch
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/migrate_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/migrate_%j.out

set -uo pipefail

# SAFETY GATE ---------------------------------------------------------------
# This script ends with `rm -rf "$SRC"`. $SRC is the WHOLE project, and since
# 2026-08-15 that includes the git repository itself.
# Push the repository to GitHub before you migrate, or move it out first.
#
# The script refuses to run without an explicit confirmation word:
#     sbatch 04_migrate_to_marinecancers.sbatch CONFIRM
if [[ "${1:-}" != "CONFIRM" ]]; then
    echo "REFUSING TO RUN." >&2
    echo "  This job copies the project to marine-cancers and then DELETES the source." >&2
    echo "  The source now contains the git repository." >&2
    echo "  Re-run with:  sbatch 04_migrate_to_marinecancers.sbatch CONFIRM" >&2
    exit 1
fi

SRC=/scratch/ahakimzadeh/Alicia_nature_data
DEST=/misc/marine-cancers/05_SEQUENCING-DATA/0501_raw/Ali/Alicia_nature_data

log() { echo "[$(date '+%F %T')] $*"; }

log "migration starting on $(hostname)"
log "groups here: $(id)"
log "source: $SRC"
log "dest  : $DEST"

# --- sanity: we must actually have the group before touching anything ----------
if ! id -nG | tr ' ' '\n' | grep -qx marinecancers; then
  log "ABORT: this node does not resolve the marinecancers group. Nothing done."
  exit 1
fi

mkdir -p "$DEST" || { log "ABORT: cannot create $DEST"; exit 1; }

# --- 1) copy ------------------------------------------------------------------
# rsync rather than cp: resumable, preserves times, and skips what is already
# there so this job can be safely re-run after an interruption.
log "copying (rsync)..."
rsync -a --info=progress2 --no-inc-recursive "$SRC/" "$DEST/"
RC=$?
if [[ $RC -ne 0 ]]; then
  log "ABORT: rsync exited $RC. Source untouched."
  exit 1
fi
log "copy finished"

# --- 2) verify every FASTQ byte-for-byte --------------------------------------
# The source FASTQs already passed ENA's published MD5 (that is what .md5.ok
# records), so proving source == dest preserves that guarantee end to end.
log "verifying FASTQ checksums source vs dest..."
FAIL=0; NCHECK=0
for f in "$SRC"/fastq/*.fastq.gz; do
  [[ -e "$f" ]] || continue
  base=$(basename "$f")
  if [[ ! -f "$DEST/fastq/$base" ]]; then
    log "  MISSING at dest: $base"; FAIL=$((FAIL+1)); continue
  fi
  s=$(md5sum "$f" | awk '{print $1}')
  d=$(md5sum "$DEST/fastq/$base" | awk '{print $1}')
  NCHECK=$((NCHECK+1))
  if [[ "$s" == "$d" ]]; then
    log "  OK   $base"
  else
    log "  FAIL $base (src=$s dst=$d)"; FAIL=$((FAIL+1))
  fi
done

# checkpoint markers must survive the move too, or the download array would
# re-download everything it already has
SRC_OK=$(find "$SRC/fastq" -maxdepth 1 -name '*.md5.ok' | wc -l)
DST_OK=$(find "$DEST/fastq" -maxdepth 1 -name '*.md5.ok' | wc -l)
log "checkpoints: source=$SRC_OK dest=$DST_OK"
[[ "$SRC_OK" != "$DST_OK" ]] && FAIL=$((FAIL+1))

{
  echo "migration report  $(date '+%F %T')"
  echo "from : $SRC"
  echo "to   : $DEST"
  echo "fastq files verified : $NCHECK"
  echo "failures             : $FAIL"
  echo "checkpoints          : src=$SRC_OK dst=$DST_OK"
} > "$DEST/migration_report.txt"

# --- 3) delete the source, but only if everything checked out -----------------
if [[ $FAIL -eq 0 && $NCHECK -gt 0 ]]; then
  log "all $NCHECK FASTQ verified. Removing source tree."
  rm -rf "$SRC"
  log "source removed: $SRC"
  echo "source deleted: yes" >> "$DEST/migration_report.txt"
  log "MIGRATION COMPLETE"
else
  log "NOT deleting source: $FAIL problem(s) found across $NCHECK file(s)."
  echo "source deleted: NO ($FAIL failures)" >> "$DEST/migration_report.txt"
  exit 1
fi
