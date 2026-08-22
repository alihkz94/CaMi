/*
 * HOST_REMOVAL — remove Cerastoderma edule reads
 * -----------------------------------------------------------------------------
 * bwa-mem2 is used, not minimap2: it is more sensitive, and it matches the
 * `bwa mem` that the study submitters used to make the original BAM files.
 *
 * -F 0x900 drops secondary and supplementary records. `keep_nonhost.sh` then
 * applies the identity rule. Read that script before you change any threshold:
 * the rule was wrong twice, and each wrong version changed the bacterial count
 * by more than a factor of 10.
 */
process HOST_REMOVAL {
    tag "$sample"

    input:
    tuple val(sample), path(t1), path(t2)
    path idx

    output:
    tuple val(sample), path("${sample}.nohost_1.fastq.gz"), path("${sample}.nohost_2.fastq.gz"), emit: nohost
    path "${sample}.bwa-mem2.cockle.log", emit: log

    script:
    // Split the CPU budget instead of over-subscribing it. The original asked for
    // task.cpus on the aligner AND 4 on each of three samtools stages, which is
    // 28 threads inside a 16-core cgroup. The samtools stages only filter and
    // recompress, so a few threads each is enough.
    def alnCpus = Math.max(1, task.cpus - 4)
    """
    bwa-mem2 mem -t ${alnCpus} -K 100000000 ${params.host_prefix} ${t1} ${t2} 2> ${sample}.bwa-mem2.cockle.log \\
      | samtools view -@ 2 -h -F 0x900 - \\
      | MIN_COV=${params.host_min_cov} MAX_DIV=${params.host_max_div} MIN_MAPQ=${params.host_min_mapq} keep_nonhost.sh \\
      | samtools collate -@ 2 -O -u - \\
      | samtools fastq -@ 2 -1 ${sample}.nohost_1.fastq.gz -2 ${sample}.nohost_2.fastq.gz -0 /dev/null -s /dev/null -n -c 6
    """
}
