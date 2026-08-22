/*
 * BUILD_HOST_INDEX — bwa-mem2 index of the cockle genome
 * -----------------------------------------------------------------------------
 * Run scripts/calibration/09_check_assembly_contamination.sbatch on any new
 * assembly BEFORE you use it. The host filter deletes every read that aligns to
 * this reference. If the assembly holds bacterial contigs, the filter deletes
 * real microbes, and nothing in the results shows that it happened.
 * GCA_947846245.1 was tested: no bacterial contigs.
 */
process BUILD_HOST_INDEX {
    tag "cockle"

    input:
    path fagz

    output:
    path "${params.host_prefix}.{0123,amb,ann,bwt.2bit.64,pac}"

    script:
    """
    # bwa-mem2 index needs an uncompressed FASTA. Delete it afterwards: mapping
    # only needs the .0123/.amb/.ann/.bwt.2bit.64/.pac files.
    pigz -dc -p ${task.cpus} ${fagz} > ${params.host_prefix}
    bwa-mem2 index ${params.host_prefix}
    rm -f ${params.host_prefix}
    """
}
