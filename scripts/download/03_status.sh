#!/usr/bin/env bash
# =============================================================================
# 03_status.sh
# -----------------------------------------------------------------------------
# PURPOSE:
#   One-glance progress report for the download array: how many runs are
#   verified, how many are still outstanding, disk used, and which array indices
#   to re-submit if anything failed.
#
# INPUTS:  ../manifests/download_list.tsv , ../fastq/*.md5.ok
# OUTPUTS: stdout only (plus ../manifests/retry_indices.txt when runs are missing)
#
# CONDA ENV: none needed
# RUN:     bash 03_status.sh
# VERSION: 1.1  (2026-08-14 — repointed at slurm-406 after the slurm-404 migration)
# =============================================================================
set -uo pipefail

# Everything (data, manifests, logs) lives on one node's local /scratch.
# That same disk has two different paths depending on where you are standing:
# /scratch on the node itself, /mnt/<node> from the login node
# (slurm-dispatcher). Resolve whichever is visible so this script works from
# either. Override with COCKLE_DATA=... .
#
# NODE is whichever node 01_download_array.sbatch currently targets
# (--nodelist=). It has moved before (slurm-404 -> slurm-406, 2026-08-14,
# see ../README.md) and may move again — if this script starts erroring
# with "No such file or directory", check the sbatch's --nodelist first.
NODE=slurm-406
if [[ -n "${COCKLE_DATA:-}" ]]; then
  DATA_ROOT="$COCKLE_DATA"
elif [[ -d "/scratch/$USER/Alicia_nature_data" ]]; then
  DATA_ROOT="/scratch/$USER/Alicia_nature_data"     # running on the node
else
  DATA_ROOT="/mnt/$NODE/$USER/Alicia_nature_data"   # running on the login node
fi

if [[ ! -d "$DATA_ROOT" ]]; then
  echo "ERROR: $DATA_ROOT does not exist." >&2
  echo "  NODE is hardcoded to '$NODE' above — if the download moved to a" >&2
  echo "  different node, update NODE here (or set COCKLE_DATA=/path/to/Alicia_nature_data)." >&2
  exit 1
fi

ROOT="$DATA_ROOT"
FASTQ_DIR="$DATA_ROOT/fastq"
MANIFEST="$ROOT/manifests/download_list.tsv"
RETRY="$ROOT/manifests/retry_indices.txt"

TOTAL=$(wc -l < "$MANIFEST")
DONE=$(find "$FASTQ_DIR" -maxdepth 1 -name '*.md5.ok' 2>/dev/null | wc -l)
LEFT=$(( TOTAL - DONE ))

echo "=============================================="
echo " PRJEB58149 raw download status"
echo "=============================================="
echo " data location : $DATA_ROOT"
echo " runs verified : $DONE / $TOTAL"
echo " outstanding   : $LEFT"
[[ $TOTAL -gt 0 ]] && echo " progress      : $(( DONE * 100 / TOTAL ))%"
echo
echo " disk used     : $(du -sh "$FASTQ_DIR" 2>/dev/null | cut -f1)"
echo " share free    : $(df -Ph "$FASTQ_DIR" | awk 'NR==2{print $4}')"
echo
echo " GB downloaded (verified runs, from manifest):"
awk -F'\t' -v d="$FASTQ_DIR" '
  { if (system("test -f " d "/" $2 ".md5.ok") == 0) g += $6; t += $6 }
  END { printf "   %.1f / %.1f GB  (%.1f%%)\n", g, t, (t>0 ? g*100/t : 0) }
' "$MANIFEST"

if (( LEFT > 0 )); then
  awk -F'\t' -v d="$FASTQ_DIR" '
    { if (system("test -f " d "/" $2 ".md5.ok") != 0) print $1 }
  ' "$MANIFEST" > "$RETRY"
  echo
  echo " outstanding array indices written to:"
  echo "   $RETRY"
  echo
  echo " Re-submitting the array as-is is safe — verified runs exit immediately:"
  echo "   sbatch $(cd "$(dirname "${BASH_SOURCE[0]}")" \&\& pwd)/01_download_array.sbatch"
fi

echo
echo " running/pending array tasks:"
squeue -u "$USER" -n cockle_dl -h -o '   %i %T %M %N' 2>/dev/null | head -20
echo "=============================================="
