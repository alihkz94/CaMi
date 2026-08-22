#!/usr/bin/env bash
# =============================================================================
# 11_compare_host_aligners.sbatch
# -----------------------------------------------------------------------------
# PURPOSE:
#   Decide, on evidence, whether host removal can move from BWA-MEM2 to minimap2.
#
#   WHY IT MATTERS: host removal is ~94% of all compute in this pipeline.
#   Measured on a real production task (FACE17_1575): 213,880,598 reads,
#   84,557 CPU-seconds, 11.9 of 12 threads busy. It is fully CPU-bound, so the
#   only way to go faster is a cheaper algorithm — not more cores, and not
#   chunking (the node is already saturated by 7 concurrent samples).
#
#   WHY IT IS NOT A FREE SWAP: minimap2 is less sensitive than BWA-MEM2 on
#   divergent reads. Less sensitive means FEWER host reads are removed, so MORE
#   cockle DNA survives into the results — exactly the failure this pipeline was
#   just fixed for (see docs/methods.md, "Host removal").
#   Speed is worthless if it reintroduces that.
#
#   THE TEST: same sample, same trimmed reads, same keep_nonhost.sh filter with
#   the same thresholds, same human-removal and dedup steps, same Kraken 2
#   settings. ONLY the host aligner differs.
#
#   Arm A (BWA-MEM2) is already done by the production run; its output is read
#   from results/04_dedup/. Arm B (minimap2) is produced here.
#
#   DECISION RULE — minimap2 wins only if it is both faster AND equivalent:
#     - surviving pairs within a few percent of Arm A
#     - bacteria count comparable
#     - Homo count NOT materially higher (higher Homo = more host leaked through)
#   If surviving pairs or Homo rise sharply, minimap2 is removing less host and
#   must be rejected however fast it is.
#
# IN : $1 = sample code already completed by the production run (default MOCE18_1282)
# OUT: aligner_comparison/<sample>/summary.txt
#
# RUN: sbatch 11_compare_host_aligners.sbatch [SAMPLE]
# VERSION: 1.0 (2026-08-16)
# =============================================================================
#SBATCH --job-name=aligncmp
#SBATCH --partition=slurm
#SBATCH --nodelist=slurm-406
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/scratch/ahakimzadeh/Alicia_nature_data
#SBATCH --output=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/aligncmp_%j.out
#SBATCH --error=/scratch/ahakimzadeh/Alicia_nature_data/slurm_logs/aligncmp_%j.out

set -euo pipefail

PROJECT="/scratch/ahakimzadeh/Alicia_nature_data"
# CODE is derived from this script's own location, not hardcoded. The repository
# was renamed to CaMi and moved out of the data root; a script that assumes either
# would break silently the next time something moves.
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PATH="/home/ahakimzadeh/miniforge3/envs/cancer/bin:$CODE/bin:$PATH"

export TMPDIR="$PROJECT/tmp"
case "$TMPDIR" in /scratch/*) ;; *) echo "ERROR: TMPDIR not on /scratch" >&2; exit 1;; esac
mkdir -p "$TMPDIR"

SAMPLE="${1:-MOCE18_1282}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
ALN=$(( THREADS - 4 ))          # identical split to the production module

# Filter thresholds — must match nextflow.config exactly, or this compares two
# things at once and proves nothing.
export MIN_COV=0.50
export MAX_DIV=0.10
export MIN_MAPQ=0

STD16="$PROJECT/kraken2_db/standard_16gb_20250402"
HUMAN_IDX="$PROJECT/reference/GRCh38.primary.mmi"
COCKLE_FA="$PROJECT/reference/cockle_V1.2.fa.gz"
COCKLE_MMI="$PROJECT/reference/cockle_V1.2.sr.mmi"

OUT="$PROJECT/aligner_comparison/$SAMPLE"
mkdir -p "$OUT"; cd "$OUT"

log(){ echo "[$(date '+%F %T')] $*"; }

# --- locate the EXACT trimmed reads the production run used -------------------
T1=$(find "$PROJECT/work" -name "${SAMPLE}.trim_1.fastq.gz" -size +1M 2>/dev/null | head -1)
T2="${T1%_1.fastq.gz}_2.fastq.gz"
[[ -s "$T1" && -s "$T2" ]] || { log "ERROR: trimmed reads for $SAMPLE not found in work/"; exit 1; }
log "input reads : $T1"

# --- Arm A ground truth: the production BWA-MEM2 result -----------------------
A1="$PROJECT/results/04_dedup/${SAMPLE}.nonhost_1.fastq.gz"
[[ -s "$A1" ]] || { log "ERROR: $SAMPLE has no completed BWA-MEM2 output yet"; exit 1; }
log "arm A (bwa-mem2) already complete: $A1"

# --- one-time: minimap2 short-read index of the cockle genome -----------------
if [[ ! -s "$COCKLE_MMI" ]]; then
    log "building minimap2 short-read index of the cockle genome (one time)"
    _t0=$(date +%s)
    minimap2 -x sr -t "$THREADS" -d "$COCKLE_MMI" "$COCKLE_FA" 2>&1 | tail -3
    log "  index build took $(( $(date +%s) - _t0 )) s"
fi

# --- Arm B: minimap2 host removal, then the identical downstream --------------
if [[ ! -s mm2.nonhost_1.fastq.gz ]]; then
    log "arm B step 1/3: minimap2 host removal"
    _t0=$(date +%s)
    minimap2 -ax sr -t "$ALN" "$COCKLE_MMI" "$T1" "$T2" 2> mm2.host.log \
      | samtools view -@ 2 -h -F 0x900 - \
      | MIN_COV=$MIN_COV MAX_DIV=$MAX_DIV MIN_MAPQ=$MIN_MAPQ keep_nonhost.sh \
      | samtools collate -@ 2 -O -u - \
      | samtools fastq -@ 2 -1 mm2.nohost_1.fastq.gz -2 mm2.nohost_2.fastq.gz \
                       -0 /dev/null -s /dev/null -n -c 6
    echo "MINIMAP2_HOST_REMOVAL $(( $(date +%s) - _t0 )) wall_s" > mm2.host.time
    cat mm2.host.time

    log "arm B step 2/3: human removal (identical to production)"
    minimap2 -ax sr -t "$ALN" "$HUMAN_IDX" mm2.nohost_1.fastq.gz mm2.nohost_2.fastq.gz 2> mm2.human.log \
      | samtools view -@ 2 -h -F 0x900 - \
      | MIN_COV=$MIN_COV MAX_DIV=$MAX_DIV MIN_MAPQ=$MIN_MAPQ keep_nonhost.sh \
      | samtools collate -@ 2 -O -u - \
      | samtools fastq -@ 2 -1 mm2.predup_1.fastq.gz -2 mm2.predup_2.fastq.gz \
                       -0 /dev/null -s /dev/null -n -c 6

    log "arm B step 3/3: dedup (identical to production)"
    fastp --thread "$THREADS" --dedup --dup_calc_accuracy 3 \
        --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
        -i mm2.predup_1.fastq.gz -I mm2.predup_2.fastq.gz \
        -o mm2.nonhost_1.fastq.gz -O mm2.nonhost_2.fastq.gz \
        -j mm2.dedup.json -h /dev/null 2> mm2.dedup.log
    rm -f mm2.nohost_*.fastq.gz mm2.predup_*.fastq.gz
fi

# --- classify BOTH arms with identical settings -------------------------------
run_k2(){   # $1=label $2=r1 $3=r2
    [[ -s "${1}.kraken2.report" ]] && { log "SKIP kraken2 $1 (done)"; return 0; }
    log "kraken2 on $1"
    kraken2 --db "$STD16" --threads "$THREADS" --paired --memory-mapping --confidence 0 \
            --report "${1}.kraken2.report" --output /dev/null "$2" "$3" 2> "${1}.kraken2.log"
}
run_k2 bwa    "$A1" "${A1%_1.fastq.gz}_2.fastq.gz"
run_k2 mm2    mm2.nonhost_1.fastq.gz mm2.nonhost_2.fastq.gz

# --- report -------------------------------------------------------------------
pairs(){ echo $(( $(pigz -dc "$1" | wc -l) / 4 )); }
tax(){ awk -F'\t' -v want="$2" '
        $4=="D" && $6 ~ want {gsub(/^ +/,"",$6); print $2; found=1}
        END{ if(!found) print 0 }' "$1" | head -1; }
homo(){ awk -F'\t' '$5==9606{print $2; f=1} END{if(!f)print 0}' "$1" | head -1; }
uncl(){ awk -F'\t' '$4=="U"{print $2; f=1} END{if(!f)print 0}' "$1" | head -1; }

BP=$(pairs "$A1"); MP=$(pairs mm2.nonhost_1.fastq.gz)
BB=$(tax bwa.kraken2.report Bacteria);  MB=$(tax mm2.kraken2.report Bacteria)
BA=$(tax bwa.kraken2.report Archaea);   MA=$(tax mm2.kraken2.report Archaea)
BH=$(homo bwa.kraken2.report);          MH=$(homo mm2.kraken2.report)
BU=$(uncl bwa.kraken2.report);          MU=$(uncl mm2.kraken2.report)

# production BWA-MEM2 timing for this sample, summed from its own batch log
BWA_LOG="$PROJECT/results/02_host_removed/${SAMPLE}.bwa-mem2.cockle.log"
BWA_REAL=$(awk 'match($0,/in ([0-9.]+) CPU sec, ([0-9.]+) real sec/,m){r+=m[2]} END{printf "%.0f", r}' "$BWA_LOG" 2>/dev/null)
MM2_REAL=$(awk '/MINIMAP2_HOST_REMOVAL/{print $2}' mm2.host.time 2>/dev/null)

{
echo "=========================================================================="
echo " HOST ALIGNER COMPARISON — $SAMPLE"
echo " identical input, identical filter (cov>=$MIN_COV div<=$MAX_DIV mapq>=$MIN_MAPQ),"
echo " identical human removal, dedup and Kraken 2 settings. Only the aligner differs."
echo "=========================================================================="
printf "%-26s %16s %16s\n" "metric" "BWA-MEM2" "minimap2"
printf "%-26s %16s %16s\n" "--------------------------" "----------------" "----------------"
printf "%-26s %16s %16s\n" "host-removal wall sec" "${BWA_REAL:-NA}" "${MM2_REAL:-NA}"
printf "%-26s %16d %16d\n" "surviving pairs"        "$BP" "$MP"
printf "%-26s %16d %16d\n" "bacteria"               "$BB" "$MB"
printf "%-26s %16d %16d\n" "archaea"                "$BA" "$MA"
printf "%-26s %16d %16d\n" "Homo sapiens"           "$BH" "$MH"
printf "%-26s %16d %16d\n" "unclassified"           "$BU" "$MU"
echo
awk -v bp="$BP" -v mp="$MP" -v bh="$BH" -v mh="$MH" -v br="${BWA_REAL:-0}" -v mr="${MM2_REAL:-0}" 'BEGIN{
  printf "  surviving pairs : minimap2 keeps %+.1f%% vs BWA-MEM2\n", (bp?100*(mp-bp)/bp:0);
  printf "  Homo            : minimap2 keeps %+.1f%% vs BWA-MEM2\n", (bh?100*(mh-bh)/bh:0);
  if(mr>0 && br>0) printf "  speed           : minimap2 is %.1fx faster on host removal\n", br/mr;
  print "";
  print "  DECISION RULE: adopt minimap2 only if it is faster AND surviving pairs";
  print "  and Homo are not materially higher. A large rise in either means it is";
  print "  removing less host — the regression this pipeline was just fixed for.";
}'
echo "=========================================================================="
} | tee summary.txt

log "DONE -> $OUT/summary.txt"
