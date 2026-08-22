/*
 * HUMAN_REMOVAL — remove human contamination
 * -----------------------------------------------------------------------------
 * minimap2 is enough here. Human contamination comes from the laboratory, so it
 * matches GRCh38 closely, and the short-read preset finds it. A bwa-mem2 index of
 * the human genome is also far larger than the minimap2 index.
 *
 * The same identity rule as HOST_REMOVAL is used. See bin/keep_nonhost.sh.
 */
process HUMAN_REMOVAL {
    tag "$sample"

    input:
    tuple val(sample), path(h1), path(h2)
    path hidx

    output:
    tuple val(sample), path("${sample}.nonhost_predup_1.fastq.gz"), path("${sample}.nonhost_predup_2.fastq.gz"), emit: nonhost
    path "${sample}.minimap2.human.log", emit: log

    script:
    def alnCpus = Math.max(1, task.cpus - 4)
    """
    minimap2 -ax sr -t ${alnCpus} ${hidx} ${h1} ${h2} 2> ${sample}.minimap2.human.log \\
      | samtools view -@ 2 -h -F 0x900 - \\
      | MIN_COV=${params.host_min_cov} MAX_DIV=${params.host_max_div} MIN_MAPQ=${params.host_min_mapq} keep_nonhost.sh \\
      | samtools collate -@ 2 -O -u - \\
      | samtools fastq -@ 2 -1 ${sample}.nonhost_predup_1.fastq.gz -2 ${sample}.nonhost_predup_2.fastq.gz -0 /dev/null -s /dev/null -n -c 6
    """
}
