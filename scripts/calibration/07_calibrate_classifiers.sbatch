#!/usr/bin/env bash
# =============================================================================
# 07_calibrate_classifiers.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Decide which classifier and which settings this project uses, on evidence.
#
#   The problem: the pilot measured 0.0532% bacteria in ERR10680552 using
#   k2_standard_16gb_20250402 with --confidence 0.1. The same run, the same
#   reads, through the shared PlusPF database at --confidence 0.1 gives 26x
#   fewer bacteria. Human reads match the pilot closely (0.380% vs 0.351%), so
#   trimming and subsampling are not the cause. It is the database/threshold pair.
#
#   This script runs every candidate on ONE identical read set and prints a
#   single table, so the choice is made on numbers rather than opinion:
#
#     Kraken2 / standard_16gb  conf 0.1   <- reproduces the pilot exactly
#     Kraken2 / standard_16gb  conf 0
#     Kraken2 / PlusPF         conf 0.1   <- what the pipeline does now
#     Kraken2 / PlusPF         conf 0.05
#     Kraken2 / PlusPF         conf 0
#     Kaiju    / nr_euk                   <- protein level, orthogonal method
#
#   Kaiju is not redundant: Kraken2 matches exact 31-base nucleotide k-mers, so
#   it misses organisms that are not close to a reference. Kaiju translates the
#   reads and searches protein space, which is far more conserved. If Kaiju also
#   reports few bacteria, the low count is real biology, not a Kraken2 artefact.
#
#   Reads are subsampled fresh from fastq/ with a fixed seed, so the input is
#   reproducible and does not depend on a Nextflow work directory.
#
#   NOTE: the pilot's denominator is TOTAL reads, before host removal. This
#   script therefore classifies TRIMMED reads, not host-removed reads, so the
#   percentages are directly comparable to PILOT_FINDINGS.md.
#
# IN : $1 = run accession (default ERR10680552, a pilot sample)
#      $2 = read pairs to sample (default 400000)
# OUT: calibration/<ACC>/summary.tsv  + one report per method
#
# RUN: sbatch 07_calibrate_classifiers.sbatch [ACC] [NPAIRS]
# VERSION: 1.0 (2026-08-14)
# =============================================================================
#SBATCH --job-name=calib
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=260G
#SBATCH --time=08:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/calib_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/calib_%j.out

set -euo pipefail

export PATH=/home/ahakimzadeh/miniforge3/envs/cancer/bin:$PATH
PROJECT="/scratch/ahakimzadeh/Alicia_nature_data"
export TMPDIR="$PROJECT/tmp"
case "$TMPDIR" in
    /scratch/*) ;;
    *) echo "ERROR: TMPDIR must be under /scratch. Got: $TMPDIR" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"

ACC="${1:-ERR10680552}"
NPAIRS="${2:-400000}"
SEED=11
THREADS="${SLURM_CPUS_PER_TASK:-16}"

PLUSPF="/slurm-databases/Kraken2/PlusPF_20250402"
STD16="$PROJECT/kraken2_db/standard_16gb_20250402"
KAIJU_DB="/slurm-databases/Kaiju/kaiju_db_nr_euk.fmi"
KAIJU_NODES="/slurm-databases/Kaiju/nodes.dmp"
KAIJU_NAMES="/slurm-databases/Kaiju/names.dmp"

OUT="$PROJECT/calibration/$ACC"
mkdir -p "$OUT"
cd "$OUT"

log() { echo "[$(date '+%F %T')] $*"; }

log "calibrating on $ACC, $NPAIRS pairs, $THREADS threads"

# --- 1. reproducible read set -------------------------------------------------
R1="$PROJECT/fastq/${ACC}_1.fastq.gz"
R2="$PROJECT/fastq/${ACC}_2.fastq.gz"
[[ -f "$PROJECT/fastq/${ACC}.md5.ok" ]] || { echo "ERROR: $ACC is not MD5-verified" >&2; exit 1; }

if [[ ! -s trim_1.fastq.gz ]]; then
    log "counting reads"
    RAW=$(( $(pigz -dc -p "$THREADS" "$R1" | wc -l) / 4 ))
    log "  $RAW pairs in the run"
    PROP=$(awk -v n="$NPAIRS" -v t="$RAW" 'BEGIN{printf "%.6f", (n<t? n/t : 1)}')
    log "subsampling proportion $PROP (seed $SEED)"
    seqkit sample -j "$THREADS" -s "$SEED" -p "$PROP" -o sub_1.fastq.gz "$R1"
    seqkit sample -j "$THREADS" -s "$SEED" -p "$PROP" -o sub_2.fastq.gz "$R2"
    log "trimming (fastp), same settings as the pipeline"
    fastp --thread "$THREADS" --dont_eval_duplication \
          -i sub_1.fastq.gz -I sub_2.fastq.gz \
          -o trim_1.fastq.gz -O trim_2.fastq.gz \
          -j fastp.json -h fastp.html 2> fastp.log
    rm -f sub_1.fastq.gz sub_2.fastq.gz
fi

PAIRS=$(( $(pigz -dc -p "$THREADS" trim_1.fastq.gz | wc -l) / 4 ))
log "classifying $PAIRS trimmed pairs"

# --- 2. Kraken2 arms ----------------------------------------------------------
run_kraken() {           # $1=label  $2=db  $3=confidence
    local label="$1" db="$2" conf="$3"
    [[ -d "$db" && -f "$db/hash.k2d" ]] || { log "SKIP $label — database not present at $db"; return 0; }
    [[ -s "${label}.report" ]] && { log "SKIP $label — already done"; return 0; }
    log "kraken2 $label (conf=$conf)"
    kraken2 --db "$db" --threads "$THREADS" --paired --memory-mapping \
            --confidence "$conf" --report "${label}.report" --output /dev/null \
            trim_1.fastq.gz trim_2.fastq.gz 2> "${label}.log"
}

run_kraken "k2_std16_conf0.1"  "$STD16"  0.1
run_kraken "k2_std16_conf0"    "$STD16"  0
run_kraken "k2_pluspf_conf0.1" "$PLUSPF" 0.1
run_kraken "k2_pluspf_conf0.05" "$PLUSPF" 0.05
run_kraken "k2_pluspf_conf0"   "$PLUSPF" 0

# --- 3. Kaiju arm (protein level) --------------------------------------------
if [[ -f "$KAIJU_DB" && ! -s kaiju.summary.tsv ]]; then
    log "kaiju (nr_euk, protein level) — the index is ~174 GB, first load is slow"
    kaiju -t "$KAIJU_NODES" -f "$KAIJU_DB" \
          -i trim_1.fastq.gz -j trim_2.fastq.gz \
          -z "$THREADS" -o kaiju.out 2> kaiju.log
    # kaiju2table -r accepts only phylum..species (no domain/superkingdom), so ask
    # for phylum rows and print the superkingdom in the path, then sum per kingdom.
    kaiju2table -t "$KAIJU_NODES" -n "$KAIJU_NAMES" -r phylum \
                -l superkingdom,phylum \
                -o kaiju.phylum.tsv kaiju.out 2>> kaiju.log
    awk -F'\t' 'NR>1 && $NF!="" {
            split($NF, a, ";"); k=a[1]; gsub(/^ +| +$/, "", k);
            if (k != "") reads[k] += $3
        }
        END { for (k in reads) printf "%s\t%d\n", k, reads[k] }' \
        kaiju.phylum.tsv | sort -k2,2nr > kaiju.summary.tsv
    gzip -f kaiju.out
else
    log "SKIP kaiju — index missing or already done"
fi

# --- 4. one comparable table --------------------------------------------------
log "building summary"
{
    printf "method\tclassified_pairs\tbacteria_pairs\tarchaea_pairs\tviruses_pairs\thuman_pairs\tpct_bacteria\n"
    for r in *.report; do
        [[ -e "$r" ]] || continue
        m="${r%.report}"
        awk -F'\t' -v m="$m" -v tot="$PAIRS" '
            $4=="U" {unc=$2}
            $4=="D" && $6 ~ /Bacteria/  {bac=$2}
            $4=="D" && $6 ~ /Archaea/   {arc=$2}
            $4=="D" && $6 ~ /Viruses/   {vir=$2}
            $5==9606                    {hum=$2}
            END {printf "%s\t%d\t%d\t%d\t%d\t%d\t%.4f\n", m, tot-unc, bac, arc, vir, hum, (tot? 100*bac/tot : 0)}
        ' "$r"
    done
    if [[ -s kaiju.summary.tsv ]]; then
        awk -F'\t' -v tot="$PAIRS" '
            $1=="Bacteria" {bac=$2} $1=="Archaea" {arc=$2} $1=="Viruses" {vir=$2} $1=="Eukaryota" {euk=$2}
            END {printf "kaiju_nr_euk\t%d\t%d\t%d\t%d\t%s\t%.4f\n", bac+arc+vir+euk, bac, arc, vir, "NA", (tot? 100*bac/tot : 0)}
        ' kaiju.summary.tsv
    fi
} > summary.tsv

log "DONE — $OUT/summary.tsv"
echo
echo "reads classified: $PAIRS pairs from $ACC"
echo "pilot baseline for ERR10680552: 0.0532% bacteria, 0.3506% human (total-read denominator)"
echo
column -t -s $'\t' summary.tsv
