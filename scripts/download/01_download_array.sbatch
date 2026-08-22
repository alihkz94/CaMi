#!/usr/bin/env bash
# =============================================================================
# 01_download_array.sbatch
#
# Download all 498 PRJEB58149 runs that have ENA FASTQ (5.24 TB), one array task
# per run, all on slurm-406 so the cohort stays on one disk.
# (Migrated from slurm-404 2026-08-14: 404 sat in a maintenance drain while its
# hardware twins 401/403/406 were already back up; 9 verified runs were rsynced
# across and re-MD5'd clean against ENA before cutover.)
# Each file is MD5-checked against ENA before the .md5.ok marker is written.
# Idempotent: re-run the same sbatch to resume.
#
# in : manifests/download_list.tsv   (idx run sample_code status prep gb)
# out: fastq/<RUN>_[12].fastq.gz + <RUN>.md5.ok , logs/<RUN>.log
#
# Notes (details in ../README.md):
#  - writes to /scratch only. NFS silently corrupts aria2c output; /tmp is the
#    root fs and filling it takes the node down.
#  - MD5 is not optional: corrupt files have appeared at the exact right size.
#  - flock guards the space check; a plain df races between concurrent tasks.
# =============================================================================
#SBATCH --job-name=cockle_dl
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --array=1-498%12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=08:00:00
# Without this, slurmstepd cannot chdir to the submit directory (that path only
# exists on the login node) and drops every task into /tmp — the filesystem that
# took slurm-001 down. Nothing here writes relative paths, but do not rely on it.
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/dl_%A_%a.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/dl_%A_%a.out

set -uo pipefail

# fall back to the literal: `set -u` would abort every task if USER were unset
OWNER="${USER:-ahakimzadeh}"
PROJECT="/scratch/$OWNER/Alicia_nature_data"          # scripts, manifests, logs
DATA_ROOT="/scratch/$OWNER/Alicia_nature_data"        # the cohort itself
MANIFEST="$PROJECT/manifests/download_list.tsv"
FASTQ_DIR="$DATA_ROOT/fastq"
LOG_DIR="$PROJECT/logs"
ARIA2C=/home/ahakimzadeh/miniforge3/envs/cancer/bin/aria2c
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"

CONN_PER_FILE=8      # aria2c -x/-s
FILES_AT_ONCE=2      # aria2c -j  (a run is 2 files: _1 and _2)
MAX_ATTEMPTS=4       # full re-download attempts before giving up on a run
RETRY_BACKOFF_S=300  # wait attempt*this between attempts (5, 10, 15 min).
                     # A proxy/network outage lasts minutes; instant retries just
                     # burn all the attempts (2026-08-14: 337 runs lost that way).
MIN_FREE_GB=2000     # leave headroom on a scratch disk other users also share

mkdir -p "$FASTQ_DIR" "$LOG_DIR" || { echo "ERROR: cannot create output dirs" >&2; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:?must be run as a Slurm array job (sbatch, not bash)}"
ROW=$(awk -F'\t' -v i="$IDX" '$1==i {print; exit}' "$MANIFEST")
if [[ -z "$ROW" ]]; then
  echo "ERROR: no manifest row with idx=$IDX in $MANIFEST" >&2
  exit 1
fi
IFS=$'\t' read -r _ ACC SAMPLE_CODE STATUS PREP GB <<< "$ROW"

LOG="$LOG_DIR/${ACC}.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

log "TASK $IDX  run=$ACC  sample=$SAMPLE_CODE  status=$STATUS  prep=$PREP  size=${GB}GB  node=$(hostname)"

# --- idempotency: a verified run is never re-downloaded --------------------------
if [[ -f "$FASTQ_DIR/${ACC}.md5.ok" ]]; then
  log "SKIP $ACC — already downloaded and MD5-verified"
  exit 0
fi

# --- space reservation (flock: a bare df check races between tasks) ---------------
NEED_GB=$(awk -v g="$GB" 'BEGIN{printf "%d", g+1}')
RESERVE_LOCK="/scratch/$OWNER/.reserve.lock"
touch "$RESERVE_LOCK" 2>/dev/null
GOT_SPACE=0
for wait_round in $(seq 1 120); do   # up to ~2 h waiting for room
  exec 9>"$RESERVE_LOCK"
  flock 9
  FREE_LOCAL=$(df -PBG "$FASTQ_DIR" | awk 'NR==2{gsub("G","",$4); print $4}')
  if (( FREE_LOCAL - NEED_GB > MIN_FREE_GB )); then
    GOT_SPACE=1
    log "space OK: ${FREE_LOCAL} GB free on $(hostname):/scratch, this run needs ~${NEED_GB} GB"
  fi
  flock -u 9; exec 9>&-
  (( GOT_SPACE == 1 )) && break
  log "  waiting for /scratch room: ${FREE_LOCAL} GB free, need ~${NEED_GB} GB + ${MIN_FREE_GB} GB floor (round $wait_round)"
  sleep 60
done

if (( GOT_SPACE == 0 )); then
  log "ERROR $ACC — no room on $(hostname):/scratch after waiting. Nothing downloaded."
  exit 1
fi

# --- ENA URLs + checksums --------------------------------------------------------
# API always prepends a run_accession column, so ftp=col2, md5=col3.
# It also returns an empty body now and then: retry, never read that as 'no FASTQ'.
FTP="" MD5=""
for attempt in 1 2 3 4 5 6; do
  META=$(curl -fsS --max-time 120 \
    "${ENA_API}?accession=${ACC}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" \
    2>/dev/null | tail -n +2)
  FTP=$(echo "$META" | cut -f2)
  MD5=$(echo "$META" | cut -f3)
  [[ -n "$FTP" && -n "$MD5" ]] && break
  log "  ENA API attempt $attempt returned nothing usable; retrying in $((attempt*15))s"
  sleep $((attempt*15))
done

if [[ -z "$FTP" || -z "$MD5" ]]; then
  log "ERROR $ACC — ENA API gave no FASTQ record after 6 attempts. This run IS listed as having FASTQ, so treat as a transient API failure and re-submit later. NOT marking verified."
  exit 1
fi

IFS=';' read -ra URLS <<< "$FTP"
IFS=';' read -ra MD5S <<< "$MD5"

# --- download + verify, with full retries ----------------------------------------
URLFILE="$FASTQ_DIR/.${ACC}.urls"
SUCCESS=0

for attempt in $(seq 1 $MAX_ATTEMPTS); do
  : > "$URLFILE"
  for u in "${URLS[@]}"; do
    [[ "$u" == http* ]] || u="https://$u"
    echo "$u" >> "$URLFILE"
  done

  log "DOWNLOAD $ACC attempt $attempt/$MAX_ATTEMPTS — ${#URLS[@]} file(s), ${CONN_PER_FILE} conn/file"
  START=$(date +%s)
  "$ARIA2C" -x"$CONN_PER_FILE" -s"$CONN_PER_FILE" -j"$FILES_AT_ONCE" \
            -c --auto-file-renaming=false --allow-overwrite=true \
            --retry-wait=15 --max-tries=8 --timeout=120 \
            --console-log-level=warn --summary-interval=0 \
            -d "$FASTQ_DIR" -i "$URLFILE" >>"$LOG" 2>&1
  RC=$?
  ELAPSED=$(( $(date +%s) - START ))
  rm -f "$URLFILE"

  if [[ $RC -ne 0 ]]; then
    log "  aria2c exited $RC after ${ELAPSED}s"
    for i in "${!URLS[@]}"; do rm -f "$FASTQ_DIR/$(basename "${URLS[$i]}")" "$FASTQ_DIR/$(basename "${URLS[$i]}").aria2"; done
    # 2026-08-14: the site web proxy (webproxy1:3128) refused connections for a
    # few minutes. aria2c returned in ~1 s, so all three attempts were spent
    # instantly and 337 runs failed together. An outage outlives an instant
    # retry — wait, and wait longer each time.
    if (( attempt < MAX_ATTEMPTS )); then
      BACKOFF=$(( attempt * RETRY_BACKOFF_S ))
      log "  waiting ${BACKOFF}s before attempt $((attempt+1)) — a network or proxy fault needs time to clear"
      sleep "$BACKOFF"
    fi
    continue
  fi

  # --- MD5 vs ENA (mandatory: corrupt files appear at the correct byte size) -----
  OK=1; BYTES=0
  for i in "${!URLS[@]}"; do
    fname=$(basename "${URLS[$i]}")
    if [[ ! -f "$FASTQ_DIR/$fname" ]]; then log "  MISSING $fname"; OK=0; continue; fi
    BYTES=$(( BYTES + $(stat -c %s "$FASTQ_DIR/$fname") ))
    got=$(md5sum "$FASTQ_DIR/$fname" | awk '{print $1}')
    if [[ "$got" == "${MD5S[$i]}" ]]; then
      log "  MD5 OK   $fname"
    else
      log "  MD5 FAIL $fname (want=${MD5S[$i]} got=$got)"; OK=0
    fi
  done

  if [[ $OK -ne 1 ]]; then
    log "  attempt $attempt failed verification — discarding and retrying"
    for i in "${!URLS[@]}"; do rm -f "$FASTQ_DIR/$(basename "${URLS[$i]}")" "$FASTQ_DIR/$(basename "${URLS[$i]}").aria2"; done
    if (( attempt < MAX_ATTEMPTS )); then sleep $(( attempt * RETRY_BACKOFF_S )); fi
    continue
  fi

  MBPS=$(( ELAPSED > 0 ? BYTES / ELAPSED / 1048576 : 0 ))
  log "  downloaded $(( BYTES / 1048576 )) MB in ${ELAPSED}s (~${MBPS} MB/s)"
  SUCCESS=1
  break
done

if [[ $SUCCESS -eq 1 ]]; then
  touch "$FASTQ_DIR/${ACC}.md5.ok"
  log "DONE $ACC — MD5 verified against ENA"
  exit 0
else
  log "ERROR $ACC — failed after $MAX_ATTEMPTS attempts. No .md5.ok written; re-submit the array to retry this run."
  for i in "${!URLS[@]}"; do rm -f "$FASTQ_DIR/$(basename "${URLS[$i]}")"; done
  exit 1
fi
