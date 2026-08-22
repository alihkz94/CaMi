/*
 * FETCH_UNMAPPED_BAM — get the 65 runs that have no FASTQ at ENA
 * -----------------------------------------------------------------------------
 * 65 of the 563 runs are submitted as BAM only, and 64 of those are tumour
 * samples — so they cannot be dropped without biasing the study.
 *
 * bin/fetch_unmapped.sh reads the BAM index over HTTP and downloads only the byte
 * ranges that hold unmapped reads. That is 0.7-3 GB for each run instead of the
 * full 55-65 GB.
 */
process FETCH_UNMAPPED_BAM {
    tag "$sample"

    input:
    tuple val(sample), val(acc)

    output:
    tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz"), emit: reads
    path "${sample}.fetch.log", emit: log
    path "${sample}.idxstats", optional: true, emit: idxstats

    script:
    // The site web proxy (webproxy1 -> webproxy-02) flaps. On 2026-08-14 it
    // refused every connection for minutes at a time. aria2c returns in about a
    // second, so an immediate retry just fails again. Retry in the TASK, with a
    // growing wait — the same fix applied to the download array.
    // Doing it here rather than in errorStrategy keeps the Nextflow head process
    // free instead of sleeping through the outage.
    """
    export OUTDIR=\$PWD
    ok=0
    for attempt in 1 2 3 4; do
        if fetch_unmapped.sh ${acc} ${sample}; then
            ok=1; break
        fi
        # fetch_unmapped.sh keeps its staging directory on failure, for debugging.
        # Each attempt makes a new one holding the whole BAM tail (~11 GB here), so
        # four attempts across 65 runs would strand terabytes on /scratch. Keep the
        # log, drop the payload.
        rm -f ${acc}_1.fastq.gz ${acc}_2.fastq.gz
        rm -rf "\${TMPDIR}"/unmap.${acc}.*
        if [ \$attempt -lt 4 ]; then
            wait_s=\$(( attempt * 300 ))
            echo "attempt \$attempt failed; waiting \${wait_s}s for the network/proxy to clear" >&2
            sleep \$wait_s
        fi
    done
    if [ \$ok -ne 1 ]; then
        echo "ERROR: ${sample} (${acc}) failed after 4 attempts" >&2
        exit 1
    fi

    # fetch_unmapped.sh names its output by accession; switch to the sample code
    # so this set matches the rest of the project.
    mv ${acc}_1.fastq.gz  ${sample}_1.fastq.gz
    mv ${acc}_2.fastq.gz  ${sample}_2.fastq.gz
    mv ${acc}.fetch.log   ${sample}.fetch.log
    [ -f ${acc}.idxstats ] && mv ${acc}.idxstats ${sample}.idxstats || true
    """
}
