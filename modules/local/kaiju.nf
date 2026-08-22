/*
 * KAIJU — taxonomy from translated protein sequence
 * -----------------------------------------------------------------------------
 * WHY THIS IS NOT REDUNDANT WITH KRAKEN 2
 *   Kraken 2 matches exact 31-base nucleotide k-mers. A marine bacterium that has
 *   no close relative in the database matches nothing, and the read is reported as
 *   unclassified. Kaiju translates the read and searches protein space, which is
 *   far more conserved, so it finds organisms that the nucleotide method misses.
 *
 *   On the same 397,334 trimmed pairs of ERR10680552:
 *       Kraken 2 / standard-16, confidence 0 :   210 bacterial pairs
 *       Kaiju    / nr_euk                    : 3,249 bacterial pairs
 *
 *   Kaiju finding MORE than the largest nucleotide database is the expected
 *   result for a marine sample, and it is the reason this step exists. Kraken 2
 *   and Bracken give the abundance backbone; Kaiju gives the sensitivity.
 *
 * RESOURCES
 *   The nr_euk index is about 187 GB and Kaiju loads it into memory — it does not
 *   memory-map it like Kraken 2. Each concurrent task therefore needs its own
 *   copy. slurm-406 has 4.1 TB, so params.kaiju_forks must stay small.
 */
process KAIJU {
    tag "$sample"

    input:
    tuple val(sample), path(n1), path(n2)

    output:
    tuple val(sample), path("${sample}.kaiju.summary.tsv"), emit: summary
    path "${sample}.kaiju.phylum.tsv", emit: phylum
    path "${sample}.kaiju.species.tsv", emit: species
    path "${sample}.kaiju.log", emit: log

    script:
    """
    kaiju -t ${params.kaiju_nodes} -f ${params.kaiju_db} \\
          -i ${n1} -j ${n2} -z ${task.cpus} -o ${sample}.kaiju.out 2> ${sample}.kaiju.log

    # kaiju2table -r accepts phylum..species only — it rejects domain and
    # superkingdom. Ask for phylum rows, print the superkingdom inside the path,
    # then add the rows per kingdom.
    # Same guard as Bracken: a sample with almost no surviving reads makes
    # kaiju2table fail, and without this the whole sample would be dropped from
    # the run. Leave an empty table instead and let the QC gate flag it.
    kaiju2table -t ${params.kaiju_nodes} -n ${params.kaiju_names} -r phylum \\
                -l superkingdom,phylum \\
                -o ${sample}.kaiju.phylum.tsv ${sample}.kaiju.out 2>> ${sample}.kaiju.log \\
        || : > ${sample}.kaiju.phylum.tsv
    kaiju2table -t ${params.kaiju_nodes} -n ${params.kaiju_names} -r species \\
                -o ${sample}.kaiju.species.tsv ${sample}.kaiju.out 2>> ${sample}.kaiju.log \\
        || : > ${sample}.kaiju.species.tsv

    awk -F'\\t' 'NR>1 && \$NF!="" {
            split(\$NF, a, ";"); k=a[1]; gsub(/^ +| +\$/, "", k);
            if (k != "") reads[k] += \$3
        }
        END { for (k in reads) printf "%s\\t%d\\n", k, reads[k] }' \\
        ${sample}.kaiju.phylum.tsv | sort -k2,2nr > ${sample}.kaiju.summary.tsv
    [ -s ${sample}.kaiju.summary.tsv ] || : > ${sample}.kaiju.summary.tsv

    # The per-read output is large and is not needed downstream. Compress it and
    # leave it in the work directory; it is not published.
    gzip -f ${sample}.kaiju.out
    """
}
