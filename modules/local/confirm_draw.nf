/*
 * CONFIRM_DRAW — a sample's detections, their reads, and the reads to BLAST
 * -----------------------------------------------------------------------------
 * The detection rule is applied to the Kraken 2 report and to Kaiju's own
 * per-read calls; every genus and species either method detects is one unit to
 * confirm. The per-read Kraken 2 calls inside each Kraken 2 unit's clade must add
 * up to the report's clade count, or the task fails. Up to --confirm_draw reads
 * per unit are drawn with a fixed seed for BLAST. See bin/confirm_draw.py.
 */
process CONFIRM_DRAW {
    tag "$sample"

    input:
    tuple val(sample), path(n1), path(report), path(calls), path(kaiju)
    path taxdump

    output:
    tuple val(sample), path("${sample}.units.tsv"), path("${sample}.k2reads.tsv.gz"),
          path("${sample}.draw.tsv"), path("${sample}.draw.fa"),
          path("${sample}.kunits.tsv"), path("${sample}.kk2reads.tsv.gz"), emit: draw
    path "${sample}.k2check.tsv", emit: check

    script:
    def merge = params.confirm_merge_family ? "--merge-family ${params.confirm_merge_family}" : ''
    """
    confirm_draw.py --sample ${sample} --report ${report} --calls ${calls} --mate1 ${n1} \\
        --kaiju ${kaiju} --taxdump ${taxdump} \\
        --draw ${params.confirm_draw} --seed ${params.seed} \\
        --min-reads ${params.confirm_min_reads} --min-frac ${params.confirm_min_frac} \\
        --exclude '${params.confirm_exclude}' ${merge}
    """
}
