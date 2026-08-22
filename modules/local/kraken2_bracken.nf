/*
 * KRAKEN2_BRACKEN — taxonomy from nucleotide k-mers, then abundance
 * -----------------------------------------------------------------------------
 * The database is passed as a PATH STRING, never as a staged input. Staging it
 * would copy the whole database into every task directory. --memory-mapping lets
 * all concurrent tasks share one copy in the page cache.
 *
 * ON --confidence: the value MUST come from params, and the default is 0.
 * The inherited pipeline used 0.1. Measured on 2,000,000 fragments of pure cockle
 * DNA — where every bacterial call is false by construction — and on real reads:
 *
 *   setting          real reads   host-only floor   signal/noise
 *   std16 conf 0       0.0529%        0.0101%          5.2x
 *   std16 conf 0.1     0.0003%        0.0009%          0.33x   <- below the floor
 *
 * At confidence 0.1 the few bacteria reported cannot be told apart from
 * misclassified host DNA. See docs/methods.md.
 */
process KRAKEN2_BRACKEN {
    tag "$sample"

    input:
    tuple val(sample), path(n1), path(n2)

    output:
    tuple val(sample), path("${sample}.kraken2.report"),
          path("${sample}.bracken.genus.tsv"), path("${sample}.bracken.species.tsv"), emit: profiled
    path "${sample}.kraken2.log",         emit: log
    path "${sample}.bracken.genus.log",   emit: glog
    path "${sample}.bracken.species.log", emit: slog

    script:
    """
    kraken2 --db ${params.kraken_db} --threads ${task.cpus} --paired --memory-mapping \\
        --confidence ${params.kraken_confidence} \\
        --report ${sample}.kraken2.report --output /dev/null ${n1} ${n2} 2> ${sample}.kraken2.log

    # Bracken can fail on a near-empty report. Keep the run moving and let the QC
    # gate flag the sample, rather than stopping everything.
    bracken -d ${params.kraken_db} -i ${sample}.kraken2.report -o ${sample}.bracken.genus.tsv \\
        -r ${params.readlen} -l G -t ${params.bracken_threshold} 2> ${sample}.bracken.genus.log \\
        || : > ${sample}.bracken.genus.tsv
    bracken -d ${params.kraken_db} -i ${sample}.kraken2.report -o ${sample}.bracken.species.tsv \\
        -r ${params.readlen} -l S -t ${params.bracken_threshold} 2> ${sample}.bracken.species.log \\
        || : > ${sample}.bracken.species.tsv
    """
}
