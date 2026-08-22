#!/usr/bin/env bash
# =============================================================================
# fetch_unmapped.sh
# -----------------------------------------------------------------------------
# PURPOSE:
#   Recover the both-mates-unmapped (non-cockle) reads from a coordinate-sorted
#   submitted BAM at ENA WITHOUT downloading the whole (tens-of-GB) file.
#
#   The SCUBA CANCERS BAMs are sorted against the C. edule V1.2 reference, so the
#   fully-unmapped ('*') reads sit in one contiguous tail at the end of the file.
#   The .bai index's largest alignment-END virtual offset marks exactly where the
#   mapped data ends and the unmapped tail begins. A BGZF virtual offset packs
#   (coffset<<16 | uoffset): coffset = the BGZF block's byte position (where we
#   start the ranged download); uoffset = how many bytes into THAT block's
#   decompressed data the first unmapped record starts (block boundaries and BAM
#   record boundaries do not coincide, so this skip is essential).
#
#   Method:  (1) grab the small .bai; (2) compute coffset/uoffset; (3) range-
#   download only [coffset:EOF] with parallel curl workers (ENA throttles per
#   connection, so we fan out ~16 ranges); (4) reconstruct an uncompressed BAM =
#   decompressed(header) + decompressed(tail)[uoffset:], pipe straight to
#   `samtools fastq`. Downloads ~0.7-3 GB instead of 55-65 GB per sample.
#
#   The result is the study's own "reads that did not map to the cockle genome"
#   set — i.e. host removal already performed by the submitters with bwa-mem
#   (Bruzos et al. Methods S6). We still re-run our bwa-mem2 host removal + human
#   removal downstream so the golden set is processed identically to batch 1.
#
# INPUTS:  $1 = ENA run accession (ERR…)   $2 = sample title used in the BAM name
# ENV:     CONN (parallel range workers, default 16)  OUTDIR (default golden_set/unmapped)
# OUTPUT:  <OUTDIR>/<ACC>_1.fastq.gz, <ACC>_2.fastq.gz  (+ .fetch.log, .idxstats)
# NEEDS:   curl, aria2c, samtools, bgzip, python3, awk   (conda env cockleseq)
# =============================================================================
set -euo pipefail
ACC="$1"; TITLE="$2"
CONN="${CONN:-12}"   # parallel range workers; EBI limits concurrent conns per IP
HERE="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${OUTDIR:-$HERE/unmapped}"; mkdir -p "$OUTDIR"
# TMPDIR must be set, and it must be on /scratch. The original line defaulted to /tmp.
# /tmp is the root filesystem: filling it stopped slurm-001 on 2026-08-12 (see README.md,
# trap 1). This script stages multi-GB BAM ranges, so there is no safe fallback. Fail loudly.
: "${TMPDIR:?TMPDIR must be set to a path on /scratch (never /tmp)}"
case "$TMPDIR" in
    /scratch/*) ;;
    *) echo "ERROR: TMPDIR must be under /scratch. Got: $TMPDIR" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"
WORK="$(mktemp -d "${TMPDIR}/unmap.${ACC}.XXXX")"
# keep WORK on failure for debugging; remove only on clean exit
cleanup(){ local rc=$?; if [[ $rc -eq 0 ]]; then rm -rf "$WORK"; else echo "[FAIL rc=$rc] kept $WORK" | tee -a "$LOG"; fi; }
trap cleanup EXIT
LOG="$OUTDIR/${ACC}.fetch.log"; : > "$LOG"
log(){ echo "[$(date +%T)] $*" | tee -a "$LOG"; }

# ENA submitted-file base path for this run (…/ERR106/ERR10681536/TITLE.bam)
# NB: ENA filereport prepends run_accession, so submitted_ftp=f2, submitted_bytes=f3.
# Both fields are ';'-joined (bam;bai) in the SAME order -> take element 1 of each.
sub=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=submitted_ftp,submitted_bytes&format=tsv" | sed -n 2p)
BAMPATH=$(echo "$sub" | cut -f2 | cut -d';' -f1)
SIZE=$(echo "$sub" | cut -f3 | cut -d';' -f1)
[[ -n "$BAMPATH" && -n "$SIZE" ]] || { log "ERROR: no submitted BAM for $ACC"; exit 1; }
BAMURL="https://${BAMPATH}"; BAIURL="${BAMURL}.bai"
log "BAM  $BAMURL  ($SIZE bytes)"

# 1) fetch .bai and compute the unmapped-tail virtual offset (block + within-block)
aria2c -x8 -s8 -j1 --file-allocation=none --console-log-level=error \
       -d "$WORK" -o idx.bai "$BAIURL" >>"$LOG" 2>&1
read VOFF COFF UOFF < <(python3 "$HERE/bai_offsets.py" "$WORK/idx.bai")
TAIL=$(( SIZE - COFF ))
log "unmapped tail: coffset=$COFF uoffset=$UOFF  tail=$TAIL bytes ($(awk -v b=$TAIL 'BEGIN{printf "%.2f", b/1e6}') MB)"
[[ "$COFF" -gt 0 && "$TAIL" -gt 0 ]] || { log "ERROR: bad offsets"; exit 1; }

# true read totals straight from the index (fast; reads only the header + index,
# not the data) — gives the real host fraction, AND the number this script cross-
# checks its own extraction against below. Required, not best-effort: a range
# fetch can return the wrong bytes at the RIGHT size (README trap: 2026-08-14, a
# proxy fault delivered a full-size but wrong-content response). idxstats reads
# the remote file independently of the range fetch below, so agreement between
# the two is real evidence, not a repeated assumption.
samtools idxstats "$BAMURL" > "$OUTDIR/${ACC}.idxstats" 2>>"$LOG"
[[ -s "$OUTDIR/${ACC}.idxstats" ]] || { log "ERROR: idxstats failed or returned nothing for $ACC"; exit 1; }
awk 'BEGIN{OFS="\t"} {m+=$3; if($1=="*")star=$4; else pm+=$4}
     END{tot=m+pm+star; printf "  totals: reads=%d mapped=%d both_unmapped=%d host=%.2f%%\n",
         tot, m, star, tot?100*m/tot:0}' "$OUTDIR/${ACC}.idxstats" | tee -a "$LOG"

# resumable single-range fetch: keeps appending the still-missing bytes of
# [a:b] to $out until complete. EBI frequently drops long connections mid-stream,
# so resuming (rather than curl --retry, which restarts the whole range) means a
# drop only costs the un-fetched remainder. Bails after $MAXTRY stalled attempts.
fetch_one(){
    local a="$1" b="$2" out="$3"          # assign first…
    local want=$(( b - a + 1 )) try=0 have laststall=0   # …then compute (set -u safe)
    while :; do
        have=0; [[ -f "$out" ]] && have=$(stat -c%s "$out")
        (( have >= want )) && break
        curl -sS -f --connect-timeout 30 --max-time 900 \
             -r "$(( a + have ))-${b}" "$BAMURL" >> "$out" 2>>"$LOG" || true
        local now; now=$(stat -c%s "$out" 2>/dev/null || echo 0)
        if (( now <= have )); then laststall=$(( laststall + 1 )); else laststall=0; fi
        try=$(( try + 1 ))
        (( laststall >= 8 )) && { log "  range $a-$b STALLED at $now/$want"; return 1; }
        (( try > 60 )) && { log "  range $a-$b too many tries ($now/$want)"; return 1; }
    done
    return 0
}

# parallel range downloader: split [start,end) into $CONN contiguous ranges,
# fetch each resumably in parallel, verify every worker succeeded and the
# reassembled size is exact.
fetch_range(){
    local start="$1" end="$2" out="$3" total=$(( $2 - $1 )) n="$CONN"
    local chunk=$(( (total + n - 1) / n )) i parts=() pids=() idx=()
    rm -f "${out}".part* "$out"
    for ((i=0;i<n;i++)); do
        local a=$(( start + i*chunk )) b=$(( start + (i+1)*chunk - 1 ))
        (( a > end-1 )) && break
        (( b > end-1 )) && b=$(( end-1 ))
        fetch_one "$a" "$b" "${out}.part${i}" &
        pids+=($!); idx+=($i)
    done
    local k rc=0
    for ((k=0;k<${#pids[@]};k++)); do
        if ! wait "${pids[$k]}"; then rc=1; log "  range worker ${idx[$k]} FAILED"; fi
    done
    (( rc == 0 )) || return 1
    for ((k=0;k<${#idx[@]};k++)); do parts+=("${out}.part${idx[$k]}"); done
    cat "${parts[@]}" > "$out"; rm -f "${parts[@]}"
    local got=$(stat -c%s "$out")
    log "  fetched $(basename "$out"): $got / $total bytes"
    (( got == total )) || { log "  SIZE MISMATCH ($got != $total)"; return 1; }
}

# 2) download the header region (first 1 MB) and the unmapped tail [COFF:SIZE),
#    both via the resumable parallel fetcher (no flaky remote htslib reads).
log "downloading header region + unmapped tail ($TAIL B) with $CONN workers…"
fetch_one 0 1048575 "$WORK/head_region.bam" || { log "ERROR: header region fetch failed"; exit 1; }
fetch_range "$COFF" "$SIZE" "$WORK/tail.bgzf"

# 3+4) reconstruct an uncompressed BAM in a pipe and stream to paired FASTQ:
#   decompressed(header)               -> BAM magic + @SQ dictionary (from the
#                                         local truncated first-slice; -H reads
#                                         only the header blocks, EOF warning ok)
#   decompressed(tail)[uoffset:]       -> the '*' unmapped alignment records
# samtools autodetects the uncompressed BAM on stdin; -f 12 -F 0x900 keeps only
# both-mates-unmapped primary reads (belt-and-suspenders; the tail is all '*').
log "reconstructing unmapped BAM and extracting paired FASTQ…"
{ samtools view -b -H "$WORK/head_region.bam" 2>>"$LOG" | bgzip -d ;
  bgzip -cd "$WORK/tail.bgzf" | tail -c +$(( UOFF + 1 )) ; } 2>>"$LOG" \
  | samtools view -@ 4 -u -f 12 -F 0x900 - 2>>"$LOG" \
  | samtools collate -@ 4 -O -u - 2>>"$LOG" \
  | samtools fastq -@ 4 -1 "$OUTDIR/${ACC}_1.fastq.gz" -2 "$OUTDIR/${ACC}_2.fastq.gz" \
                   -0 /dev/null -s /dev/null -n -c 6 2>>"$LOG"
p1=$(( $(zcat "$OUTDIR/${ACC}_1.fastq.gz" | wc -l) / 4 ))
p2=$(( $(zcat "$OUTDIR/${ACC}_2.fastq.gz" | wc -l) / 4 ))
[[ "$p1" -eq "$p2" ]] || { log "ERROR: mate count mismatch R1=$p1 R2=$p2"; exit 1; }

# HARD cross-check against the index's own '*' (unmapped) count. Every pair we
# extracted has BOTH mates unmapped (samtools -f 12 above), so it must be a
# SUBSET of the '*' bucket: 2*p1 can never legitimately exceed star. If it does,
# the reconstructed tail included data it should not have — the exact failure
# mode a wrong byte range or a tampered proxy response produces, and the one a
# size check alone cannot see (the file can be the right size and still wrong).
star=$(awk '$1=="*"{print $4}' "$OUTDIR/${ACC}.idxstats" 2>/dev/null)
[[ -n "$star" ]] || { log "ERROR: no '*' row in idxstats for $ACC — cannot verify the extraction"; exit 1; }
if (( 2 * p1 > star )); then
    log "ERROR: extracted $p1 pairs (=$(( 2*p1 )) reads) but the index reports only $star unmapped reads."
    log "  The range fetch returned data it should not have. Treating this as corrupt. Not writing .md5-equivalent success."
    exit 1
fi
# A shortfall is expected sometimes (a read can be unmapped while its mate is
# mapped; -f 12 correctly excludes those). Flag a LARGE shortfall for a human to
# look at, but do not fail the run on it alone — it is not proof of corruption.
if (( star > 0 && 4 * p1 < star )); then
    log "WARN: only $p1 pairs recovered against $star unmapped reads in the index (<50%). Inspect before trusting this sample."
fi
log "DONE $ACC -> ${ACC}_[12].fastq.gz  ($p1 read pairs; idx '*'=$star reads, upper bound $(( star/2 )) pairs)"
