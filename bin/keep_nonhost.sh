#!/usr/bin/env bash
# =============================================================================
# keep_nonhost.sh — decide which reads are NOT host, from a SAM stream
# -----------------------------------------------------------------------------
# Reads SAM on stdin, writes SAM on stdout, and keeps only the reads that are not
# confidently host. `samtools fastq` downstream sends unpaired survivors to -s,
# so a fragment survives only when NEITHER mate is host. That is correct: two
# mates come from one DNA fragment, so one host mate condemns the pair.
#
# THE RULE
#   A read is host when it aligns CONVINCINGLY:
#       mapped  AND  coverage >= MIN_COV  AND  divergence <= MAX_DIV
#   where
#       coverage   = aligned query bases / total query bases   (from the CIGAR)
#       divergence = NM / aligned query bases
#   Everything else is kept: unmapped reads, short local hits, and divergent hits.
#
#   MAPQ IS NOT USED BY DEFAULT (MIN_MAPQ=0). This is the important part, and it
#   was learned the hard way — see the history below.
#
# HISTORY — two wrong rules came before this one
#
#   1. `samtools view -f 12` (inherited). Keeps a pair only when BOTH mates are
#      completely unmapped. bwa-mem2 reports an alignment for almost any read
#      with a single seed match, so one spurious hit threw the pair away.
#      Measured on ERR10680552: 4,138 pairs survived, 3 bacterial pairs.
#
#   2. `MAPQ >= 20 AND low divergence` (this script, v1.0). Better, but wrong in
#      the opposite direction. MAPQ < 20 means the aligner cannot choose between
#      several equally good copies — a REPEAT. It does not mean the read is
#      foreign. Molluscan genomes are full of repeats, so this rule called a
#      large block of ordinary cockle sequence "microbial".
#      Measured: 65,436 pairs survived, 41 bacterial pairs.
#
#   The cross-tab that settled it (ERR10680552, 794,668 mates vs the cockle
#   assembly):
#
#     unmapped                                 12,227    1.54%
#     MAPQ>=20, div<=.10   (host, every rule) 632,174   79.55%
#     MAPQ>=20, div>.10                        17,428    2.19%
#     MAPQ<20,  cov>=.5, div<=.10             120,559   15.17%   <- disputed
#     MAPQ<20,  cov>=.5, div>.10                3,623    0.46%
#     MAPQ<20,  cov<.5     (short local)        8,657    1.09%
#
#   The disputed 120,559 mates align over most of their length at high identity
#   but keep MAPQ<20. Rule 2 called every one of them microbial. Classifying them
#   directly with Kraken2 gave the answer:
#
#     unclassified   51,814 pairs   99.33%   <- cockle: absent from the database
#     Homo sapiens      320 pairs    0.61%   <- human repeats
#     bacteria           26 pairs    0.05%   <- 2-3 reads each, unrelated families
#
#   So the disputed block is host sequence, and rule 2 leaked it into the
#   results. It is what made `Homo` the dominant genus in the QC report. The 26
#   bacterial pairs are scattered singletons across unrelated families, which is
#   the signature of a spurious LCA assignment, not of a real community.
#
# ON THE THRESHOLDS
#   C. edule is highly polymorphic and the reference is one individual, so real
#   host reads carry real mismatches. MAX_DIV 0.10 allows for that. MIN_COV 0.50
#   keeps reads that touch the assembly only over a short stretch, which is what
#   a genuine foreign read with one conserved domain looks like.
#
#   The assembly was checked for bacterial contigs before this rule was adopted
#   (scripts/calibration/09_check_assembly_contamination.sbatch). A dirty assembly would make
#   this rule delete real microbes, and the deletion would be invisible.
#
# ENV: MIN_COV (default 0.50)  MAX_DIV (default 0.10)  MIN_MAPQ (default 0, off)
# USE: samtools view -h in.bam | keep_nonhost.sh | samtools collate -O -u - | samtools fastq ...
# VERSION: 2.0 (2026-08-15)
# =============================================================================
set -euo pipefail

MIN_COV="${MIN_COV:-0.50}"
MAX_DIV="${MAX_DIV:-0.10}"
MIN_MAPQ="${MIN_MAPQ:-0}"

exec awk -v min_cov="$MIN_COV" -v max_div="$MAX_DIV" -v min_mapq="$MIN_MAPQ" '
    BEGIN { OFS = "\t" }
    /^@/ { print; next }
    {
        host = 0
        if ( ! and($2, 4) ) {                       # 0x4 = unmapped
            nm = -1
            for (i = 12; i <= NF; i++)
                if ($i ~ /^NM:i:/) { split($i, a, ":"); nm = a[3]; break }

            # Walk the CIGAR: M/I/=/X consume the query and are aligned;
            # S is soft-clipped, so it is query but not aligned. H is not in SEQ.
            aln = 0; qlen = 0; num = ""
            len6 = length($6)
            for (j = 1; j <= len6; j++) {
                c = substr($6, j, 1)
                if (c ~ /[0-9]/) { num = num c; continue }
                v = num + 0; num = ""
                if (c == "M" || c == "I" || c == "=" || c == "X") { aln += v; qlen += v }
                else if (c == "S")                                { qlen += v }
            }
            if (qlen == 0) qlen = length($10)       # CIGAR "*" — fall back to SEQ

            cov = (qlen > 0) ? aln / qlen : 0
            # No NM tag means divergence cannot be judged. A read that aligns over
            # most of its length is treated as host, which is the safe direction:
            # a leaked host read is classified and pollutes the table.
            div = (aln > 0 && nm >= 0) ? nm / aln : 0

            if (cov >= min_cov && div <= max_div && $5 >= min_mapq) host = 1
        }
        if ( ! host ) print
    }
'
