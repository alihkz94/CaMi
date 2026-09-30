/*
 * CONFIRM_COMBINE — the three legs, one row per (sample, taxon) detection
 * -----------------------------------------------------------------------------
 * Refuses unless every sample it is given has all its pieces, every BLAST chunk
 * is present, and Kaiju and Kraken 2 classified the same number of pairs.
 * See bin/confirm_combine.py and docs/methods.md, section 9.
 */
process CONFIRM_COMBINE {

    input:
    path per_sample
    path queries
    path hits
    path samples
    path taxdump

    output:
    path "confirmation_long.tsv", emit: long
    path "summary.tsv",           emit: summary
    path "blast_reads.tsv.gz"
    path "params.json"
    path "per_sample"

    script:
    def merge = params.confirm_merge_family ? "--merge-family ${params.confirm_merge_family}" : ''
    """
    confirm_combine.py --samples ${samples} --taxdump ${taxdump} --jobs ${task.cpus} ${merge} \\
        --min-agree ${params.confirm_min_agree} --min-frac ${params.confirm_min_share} \\
        --genus-pident ${params.confirm_genus_pident} --species-pident ${params.confirm_species_pident} \\
        --min-qcov ${params.confirm_min_qcov}
    """
}
