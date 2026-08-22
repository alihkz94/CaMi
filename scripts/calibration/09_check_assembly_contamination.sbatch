#!/usr/bin/env bash
# =============================================================================
# 09_check_assembly_contamination.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Make sure the cockle reference is clean before we use it to remove host reads.
#
#   The host filter deletes every read that aligns convincingly to the cockle
#   assembly. If the assembly itself contains bacterial contigs — which happens
#   often with molluscs, because the animal is full of bacteria when it is
#   sequenced — then the filter deletes real microbial reads, and the deletion is
#   invisible in the results.
#
#   This job classifies the assembly contigs themselves. A clean assembly gives
#   almost only unclassified contigs. Any bacterial contig is listed by name and
#   length, so we can judge the size of the problem.
#
# OUT: reference/contamination_check/cockle.report + bacterial_contigs.txt
#
# RUN: sbatch 09_check_assembly_contamination.sbatch
# VERSION: 1.0 (2026-08-15)
# =============================================================================
#SBATCH --job-name=asmcheck
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/asmcheck_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/asmcheck_%j.out

set -euo pipefail
export PATH=/home/ahakimzadeh/miniforge3/envs/cancer/bin:$PATH
PROJECT="/scratch/ahakimzadeh/Alicia_nature_data"
export TMPDIR="$PROJECT/tmp"
case "$TMPDIR" in /scratch/*) ;; *) echo "ERROR: TMPDIR not on /scratch" >&2; exit 1;; esac
mkdir -p "$TMPDIR"

REF="$PROJECT/reference/cockle_V1.2.fa.gz"   # kraken2 reads gzip directly
OUT="$PROJECT/reference/contamination_check"
mkdir -p "$OUT"; cd "$OUT"

[[ -s "$REF" ]] || { echo "ERROR: reference missing: $REF" >&2; exit 1; }

echo "[$(date '+%F %T')] classifying assembly contigs"
kraken2 --db "$PROJECT/kraken2_db/standard_16gb_20250402" \
        --threads "${SLURM_CPUS_PER_TASK:-8}" --memory-mapping --confidence 0 \
        --report cockle.report --output cockle.out "$REF" 2> cockle.log

echo
echo "=== contig-level summary ==="
awk -F'\t' '$4=="U"{printf "unclassified\t%d contigs\t%.2f%%\n",$2,$1}
            $4=="D"{gsub(/^ +/,"",$6); printf "domain %s\t%d contigs\t%.2f%%\n",$6,$2,$1}' cockle.report

echo
echo "=== bacterial / archaeal contigs, longest first ==="
# column 3 of kraken2 --output is the sequence id, column 4 the length
awk -F'\t' '$1=="C"{print $2"\t"$3"\t"$4}' cockle.out > classified_contigs.tsv
if [[ -s classified_contigs.tsv ]]; then
    sort -k3,3nr classified_contigs.tsv | head -25 | \
        awk -F'\t' '{printf "%-24s taxid=%-10s %s bp\n", $1, $2, $3}'
    echo
    echo "total classified contigs: $(wc -l < classified_contigs.tsv)"
    echo "total classified bases:   $(awk -F'\t' '{s+=$3} END{printf "%d", s}' classified_contigs.tsv)"
    echo "assembly bases:           $(pigz -dc "$REF" | awk '!/^>/{s+=length($0)} END{printf "%d", s}')"
else
    echo "none — the assembly is clean"
fi
echo
echo "[$(date '+%F %T')] DONE"
