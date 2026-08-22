/*
 * DEDUP — remove duplicate read pairs
 * -----------------------------------------------------------------------------
 * This is the last step that changes the reads. Its output IS the microbial read
 * set, and it is always published: after host and human removal the set is a few
 * percent of the raw data, so the storage cost is small, and it is the input for
 * any later assembly, Kaiju run, or targeted BLAST search.
 *
 * --dup_calc_accuracy 3, not 6: level 6 runs out of memory (README trap 6).
 */
process DEDUP {
    tag "$sample"

    input:
    tuple val(sample), path(p1), path(p2)

    output:
    tuple val(sample), path("${sample}.nonhost_1.fastq.gz"), path("${sample}.nonhost_2.fastq.gz"), emit: dedup
    path "${sample}.dedup.json", emit: qc
    path "${sample}.dedup.log",  emit: log

    script:
    """
    fastp --thread ${task.cpus} --dedup --dup_calc_accuracy 3 \\
        --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \\
        -i ${p1} -I ${p2} -o ${sample}.nonhost_1.fastq.gz -O ${sample}.nonhost_2.fastq.gz \\
        -j ${sample}.dedup.json -h /dev/null 2> ${sample}.dedup.log
    """
}
