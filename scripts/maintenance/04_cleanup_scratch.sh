#!/usr/bin/env bash
# =============================================================================
# 04_cleanup_scratch.sh
# -----------------------------------------------------------------------------
# PURPOSE:
#   Find and delete leftover staging directories belonging to this user on every
#   compute node's /scratch. The download array cleans up after itself via an
#   EXIT/INT/TERM trap, but a hard kill (scancel -9, node crash, OOM) bypasses it.
#   Run this after any cancelled or crashed run so nothing of ours is left on
#   shared node disks — and so you can answer "have you cleaned up?" with data.
#
#   Every node's /scratch is NFS-visible from the login node as /mnt/<nodename>,
#   so this needs no Slurm allocation and takes seconds.
#
# INPUTS:  none
# OUTPUTS: stdout report; deletes /mnt/slurm-*/$USER/cockle_dl.*
#
# CONDA ENV: none needed
#
# RUN:
#   bash 04_cleanup_scratch.sh              # report + delete leftovers
#   bash 04_cleanup_scratch.sh --dry-run    # report only, delete nothing
#
# VERSION: 2.0  (2026-08-12)  use /mnt/<node> login-node mounts instead of srun
# =============================================================================
set -uo pipefail

DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

echo "=============================================="
echo " scratch cleanup for $USER"
[[ $DRY -eq 1 ]] && echo " MODE: dry run — nothing will be deleted"
echo "=============================================="

TOTAL_FOUND=0
for M in /mnt/slurm-*; do
  NODE=$(basename "$M")
  [[ -d "$M" ]] || continue

  # Skip mounts that are not actually readable (node down / not exported).
  ls "$M" >/dev/null 2>&1 || { echo "-- $NODE : unreachable, skipped"; continue; }

  DIRS=$(find "$M/$USER" -maxdepth 1 -name 'cockle_dl.*' 2>/dev/null)
  if [[ -z "$DIRS" ]]; then
    continue
  fi

  echo "-- $NODE"
  while IFS= read -r D; do
    [[ -n "$D" ]] || continue
    SZ=$(du -sh "$D" 2>/dev/null | cut -f1)
    echo "   LEFTOVER $D ($SZ)"
    TOTAL_FOUND=$((TOTAL_FOUND+1))
    if [[ $DRY -eq 0 ]]; then
      rm -rf "$D" && echo "   removed  $D"
    fi
  done <<< "$DIRS"
done

echo
if (( TOTAL_FOUND == 0 )); then
  echo " No leftover staging directories found — scratch is clean."
else
  [[ $DRY -eq 0 ]] && echo " Removed $TOTAL_FOUND leftover directory/ies." \
                   || echo " Found $TOTAL_FOUND leftover directory/ies (dry run — not deleted)."
fi

echo
echo " Current /scratch free space per node:"
for M in /mnt/slurm-*; do
  NODE=$(basename "$M")
  FREE=$(df -Ph "$M" 2>/dev/null | awk 'NR==2{print $4}')
  USED_BY_US=$(du -sh "$M/$USER" 2>/dev/null | cut -f1)
  [[ -n "$FREE" ]] && printf "   %-14s free=%-7s ours=%s\n" "$NODE" "$FREE" "${USED_BY_US:-0}"
done
echo "=============================================="
