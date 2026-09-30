/*
 * CONFIRM_BLAST_PREP — one query set for the whole cohort, cut into chunks
 * -----------------------------------------------------------------------------
 * A BLAST search against nt costs one full scan of the database plus a small
 * amount per read, so the drawn reads of all samples are pooled and searched in
 * chunks of --confirm_chunk reads rather than one search per sample.
 */
process CONFIRM_BLAST_PREP {

    input:
    path draws

    output:
    path "queries.tsv", emit: queries
    path "chunk_*.fa",  emit: chunks, optional: true

    script:
    """
    confirm_blast_prep.py --chunk ${params.confirm_chunk} --seed ${params.seed}
    """
}
