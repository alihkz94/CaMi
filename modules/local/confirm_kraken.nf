/*
 * CONFIRM_KRAKEN — the per-read Kraken 2 calls behind a sample's report
 * -----------------------------------------------------------------------------
 * KRAKEN2_BRACKEN keeps only the report. This reruns Kraken 2 with the same
 * database and flags, adds --output, and refuses unless the new report is
 * byte-identical to the one the profile was built from. The per-read calls are
 * then known to be the ones that produced the reported counts.
 *
 * Kraken 2 memory-maps its database, so the rerun costs seconds to minutes and
 * shares the page cache with KRAKEN2_BRACKEN.
 */
process CONFIRM_KRAKEN {
    tag "$sample"

    input:
    tuple val(sample), path(n1), path(n2), path(report)

    output:
    tuple val(sample), path("${sample}.kraken2.calls.tsv.gz"), path("${sample}.pairs"), emit: calls

    script:
    """
    kraken2 --db ${params.kraken_db} --threads ${task.cpus} --paired --memory-mapping \\
        --confidence ${params.kraken_confidence} \\
        --report rerun.report --output rerun.out ${n1} ${n2} 2> rerun.log

    if ! cmp -s rerun.report ${report}; then
        echo "${sample}: the rerun report differs from ${report}. The database or the" >&2
        echo "flags changed since profiling; rerun PROFILE before confirming." >&2
        diff ${report} rerun.report | head -20 >&2
        exit 1
    fi

    # Classified pairs only (read id, taxid). The pair count covers every line, so
    # an unclassified pair is known to be unclassified, not missing.
    awk -F'\\t' '\$1=="C"{print \$2"\\t"\$3}' rerun.out | gzip -1 > ${sample}.kraken2.calls.tsv.gz
    wc -l < rerun.out | tr -d ' ' > ${sample}.pairs
    rm -f rerun.out
    """
}
