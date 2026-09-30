/*
 * CONFIRM_BLAST — megablast of one query chunk against nt
 * -----------------------------------------------------------------------------
 * The database is passed as a PATH STRING and never staged (nt is ~700 GB).
 * taxdb.btd/.bti must sit beside it for the staxids column.
 *
 * Speed: each search scans nt once, then costs ~0.1 s per read. The scan takes
 * ~24 min from a cold file system and ~2 min with nt in the page cache. On Slurm
 * the cache is charged to the task that first reads it, so a task smaller than
 * nt evicts it while scanning. docs/methods.md, section 9, gives the fix.
 */
process CONFIRM_BLAST {
    tag "${chunk.baseName}"

    input:
    path chunk

    output:
    path "${chunk.baseName}.blast.tsv.gz", emit: hits

    script:
    """
    export BLASTDB=\$(dirname ${params.blast_db})
    blastn -task megablast -db ${params.blast_db} -query ${chunk} -num_threads ${task.cpus} \\
           -evalue 1e-5 -max_target_seqs 10 \\
           -outfmt "6 qseqid sacc staxids pident length qlen evalue bitscore" \\
           -out hits.tsv

    # Every line has 8 columns and every hit belongs to a query of this chunk. A
    # query without a hit writes no line; confirm_combine.py counts it as nohit.
    awk -F'\\t' 'NF!=8{bad=1} END{exit bad}' hits.tsv
    grep '^>' ${chunk} | cut -c2- | sort > q.ids
    cut -f1 hits.tsv | sort -u > hit.ids
    if [ -n "\$(comm -13 q.ids hit.ids)" ]; then echo "hits for queries not in ${chunk}" >&2; exit 1; fi

    gzip -c hits.tsv > ${chunk.baseName}.blast.tsv.gz
    rm -f hits.tsv
    """
}
