#!/usr/bin/env bash
# =============================================================================
# 08_test_host_rule.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Decide the host-removal rule on evidence, not on opinion.
#
#   Three rules were considered for "is this read host?":
#     OLD  (inherited)  samtools view -f 12 — keep a pair only if BOTH mates are
#                       completely unmapped. Throws away any pair with a single
#                       spurious seed hit. Measured: 4,138 pairs, 3 bacteria.
#     MAPQ (v1.0)       host if mapped AND MAPQ>=20 AND low divergence.
#                       Measured: 65,436 pairs, 41 bacteria.
#     COV  (v2.0)       host if mapped AND aligned over >=50% of the read at
#                       <=10% divergence, WHATEVER the MAPQ.
#
#   Why COV is likely correct: MAPQ<20 means the aligner cannot choose between
#   equally good copies — a repeat. It does NOT mean the read is foreign. The
#   cross-tab on ERR10680552 (794,668 mates) shows the size of the problem:
#
#     unmapped                                12,227   1.54%
#     MAPQ>=20, div<=.10   (host, all rules) 632,174  79.55%
#     MAPQ>=20, div>.10                       17,428   2.19%
#     MAPQ<20,  cov>=.5, div<=.10            120,559  15.17%  <- disputed
#     MAPQ<20,  cov>=.5, div>.10               3,623   0.46%
#     MAPQ<20,  cov<.5   (short local)         8,657   1.09%
#
#   The disputed 120,559 mates align convincingly to the cockle genome but keep
#   MAPQ<20. The MAPQ rule calls them microbial. This job classifies them
#   directly. If they come back as host or unclassified, the MAPQ rule leaks
#   host sequence and COV is right.
#
# IN : calibration/<ACC>/bucket3/  (made by the caller)
# OUT: calibration/<ACC>/bucket3/b3.report
#
# RUN: sbatch 08_test_host_rule.sbatch [ACC]
# VERSION: 1.0 (2026-08-15)
# =============================================================================
#SBATCH --job-name=hostrule
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/hostrule_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/hostrule_%j.out

set -euo pipefail
export PATH=/home/ahakimzadeh/miniforge3/envs/cancer/bin:$PATH
PROJECT="/scratch/ahakimzadeh/Alicia_nature_data"
export TMPDIR="$PROJECT/tmp"
case "$TMPDIR" in /scratch/*) ;; *) echo "ERROR: TMPDIR not on /scratch" >&2; exit 1;; esac
mkdir -p "$TMPDIR"

ACC="${1:-ERR10680552}"
cd "$PROJECT/calibration/$ACC/bucket3"

kraken2 --db "$PROJECT/kraken2_db/standard_16gb_20250402" \
        --threads "${SLURM_CPUS_PER_TASK:-8}" --paired --memory-mapping --confidence 0 \
        --report b3.report --output /dev/null b3_1.fq.gz b3_2.fq.gz 2> b3.log

echo "=== identity of the disputed low-MAPQ, high-coverage reads ==="
awk -F'\t' '
  $4=="U"           {printf "unclassified\t%d\t%.2f%%\n", $2, $1}
  $4=="R"           {printf "classified(root)\t%d\t%.2f%%\n", $2, $1}
  $4=="D"           {gsub(/^ +/,"",$6); printf "domain %s\t%d\t%.2f%%\n", $6, $2, $1}
  $5==9606          {printf "  Homo sapiens\t%d\t%.2f%%\n", $2, $1}
  $5==6544||$5==6447{gsub(/^ +/,"",$6); printf "  %s\t%d\t%.2f%%\n", $6, $2, $1}
' b3.report
echo
echo "=== top 15 named taxa ==="
awk -F'\t' '$4 ~ /^[SGF]$/ {gsub(/^ +/,"",$6); printf "%8d  %-6s %s\n", $2, $4, $6}' b3.report | sort -k1,1nr | head -15
