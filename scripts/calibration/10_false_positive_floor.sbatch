#!/usr/bin/env bash
# =============================================================================
# 10_false_positive_floor.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Measure how often PURE HOST DNA is called bacterial. This sets the false
#   positive floor of the classifier, and it decides the confidence threshold.
#
#   Why we must know this: the contig-level check (script 09) showed that 99.998%
#   of the cockle assembly classifies at --confidence 0, and it classifies as
#   Homo sapiens — because the standard-16 database holds exactly one eukaryote,
#   human. Cockle is not in the database, so cockle sequence cannot be identified
#   correctly. It can only be unclassified or WRONG.
#
#   The pilot reported 0.0532% bacteria at --confidence 0. If pure cockle DNA
#   also gives about 0.05% bacteria, then the pilot measured noise, and every
#   sample in this study will "find" the same fake community.
#
#   METHOD: cut the assembly into read-length fragments, take a random sample,
#   and classify them at several thresholds. The fragments contain no microbe by
#   construction, so EVERY bacterial call is a false positive.
#
#   Note the fragments have no sequencing errors, so they match the database more
#   easily than real reads do. The measured rate is therefore an upper bound.
#
# OUT: calibration/false_positive_floor/summary.tsv
#
# RUN: sbatch 10_false_positive_floor.sbatch [NFRAG]
# VERSION: 1.0 (2026-08-15)
# =============================================================================
#SBATCH --job-name=fpfloor
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/fpfloor_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/fpfloor_%j.out

set -euo pipefail
export PATH=/home/ahakimzadeh/miniforge3/envs/cancer/bin:$PATH
PROJECT="/scratch/ahakimzadeh/Alicia_nature_data"
export TMPDIR="$PROJECT/tmp"
case "$TMPDIR" in /scratch/*) ;; *) echo "ERROR: TMPDIR not on /scratch" >&2; exit 1;; esac
mkdir -p "$TMPDIR"

NFRAG="${1:-2000000}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
REF="$PROJECT/reference/cockle_V1.2.fa.gz"
STD16="$PROJECT/kraken2_db/standard_16gb_20250402"
PLUSPF="/slurm-databases/Kraken2/PlusPF_20250402"

OUT="$PROJECT/calibration/false_positive_floor"
mkdir -p "$OUT"; cd "$OUT"

if [[ ! -s frags.fa.gz ]]; then
    echo "[$(date '+%F %T')] cutting the assembly into 150 bp fragments"
    # -s 150 -W 150 gives non-overlapping windows, so no fragment is counted twice.
    seqkit sliding -j "$THREADS" -s 150 -W 150 "$REF" 2> sliding.log \
      | seqkit sample -j "$THREADS" -s 11 -n "$NFRAG" -o frags.fa.gz 2> sample.log
fi
TOT=$(seqkit stats -T frags.fa.gz | awk 'NR==2{print $4}')
echo "[$(date '+%F %T')] $TOT host fragments — every bacterial call below is a false positive"

run_k2() {                      # $1=label $2=db $3=conf
    local label="$1" db="$2" conf="$3"
    [[ -d "$db" && -f "$db/hash.k2d" ]] || { echo "SKIP $label — no database at $db"; return 0; }
    [[ -s "${label}.report" ]] && { echo "SKIP $label — already done"; return 0; }
    echo "[$(date '+%F %T')] kraken2 $label conf=$conf"
    kraken2 --db "$db" --threads "$THREADS" --memory-mapping --confidence "$conf" \
            --report "${label}.report" --output /dev/null frags.fa.gz 2> "${label}.log"
}

run_k2 std16_conf0     "$STD16"  0
run_k2 std16_conf0.05  "$STD16"  0.05
run_k2 std16_conf0.1   "$STD16"  0.1
run_k2 pluspf_conf0    "$PLUSPF" 0
run_k2 pluspf_conf0.05 "$PLUSPF" 0.05
run_k2 pluspf_conf0.1  "$PLUSPF" 0.1

{
  printf "method\tbacteria_frags\tarchaea_frags\tviruses_frags\thuman_frags\tunclassified\tFALSE_POS_pct_bact\n"
  for r in *.report; do
      [[ -e "$r" ]] || continue
      awk -F'\t' -v m="${r%.report}" -v tot="$TOT" '
        $4=="U"{unc=$2}
        $4=="D" && $6 ~ /Bacteria/{bac=$2}
        $4=="D" && $6 ~ /Archaea/ {arc=$2}
        $4=="D" && $6 ~ /Viruses/ {vir=$2}
        $5==9606{hum=$2}
        END{printf "%s\t%d\t%d\t%d\t%d\t%d\t%.4f\n", m,bac,arc,vir,hum,unc,(tot?100*bac/tot:0)}
      ' "$r"
  done
} > summary.tsv

echo
echo "=== FALSE POSITIVE FLOOR — pure cockle DNA, no microbes present ==="
column -t -s $'\t' summary.tsv
echo
echo "Compare with the pilot's reported 0.0532% bacteria on real reads."
