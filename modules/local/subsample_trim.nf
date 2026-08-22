/*
 * SUBSAMPLE_TRIM — optional subsampling, then adapter and quality trimming
 * -----------------------------------------------------------------------------
 * Use --subsample N to work on a fixed number of read pairs. Use 0, the default,
 * to keep every read. A subsample is for tests only: at 400,000 pairs a sample
 * shows no microbial signal above the noise.
 */
process SUBSAMPLE_TRIM {
    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}.trim_1.fastq.gz"), path("${sample}.trim_2.fastq.gz"), emit: trimmed
    path "${sample}.fastp.json", emit: qc
    path "${sample}.fastp.html", emit: html
    path "${sample}.fastp.log",  emit: log

    script:
    // The original counted reads with `zcat | wc -l` on every sample, before
    // testing whether subsampling was even requested. At subsample=0 (the default
    // for a full run) that decompresses every FASTQ end-to-end for a number that
    // is then discarded — about 10 TB of pointless I/O across the cohort.
    // Count only when a cap is actually set, and use all the cores to do it.
    """
    in1=${r1}; in2=${r2}
    if [ ${params.subsample} -gt 0 ]; then
        raw=\$(( \$(pigz -dc -p ${task.cpus} ${r1} | wc -l) / 4 ))
        if [ \$raw -gt ${params.subsample} ]; then
            prop=\$(awk -v n=${params.subsample} -v t=\$raw 'BEGIN{printf "%.6f", n/t}')
            seqkit sample -j ${task.cpus} -s ${params.seed} -p \$prop -o sub_1.fastq.gz ${r1}
            seqkit sample -j ${task.cpus} -s ${params.seed} -p \$prop -o sub_2.fastq.gz ${r2}
            in1=sub_1.fastq.gz; in2=sub_2.fastq.gz
        fi
    fi
    fastp --thread ${task.cpus} --dont_eval_duplication \\
        -i \$in1 -I \$in2 -o ${sample}.trim_1.fastq.gz -O ${sample}.trim_2.fastq.gz \\
        -j ${sample}.fastp.json -h ${sample}.fastp.html 2> ${sample}.fastp.log
    """
}
