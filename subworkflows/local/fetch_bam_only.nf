/*
 * FETCH_BAM_ONLY — recover the 65 runs that ENA holds only as BAM
 * -----------------------------------------------------------------------------
 * 64 of these 65 runs are tumour samples. To drop them would bias the whole
 * disease comparison, so they are fetched rather than skipped.
 *
 * The runs are fetched ONE AT A TIME (maxForks 1 in conf/base.config). ENA
 * throttles per client, and this step normally runs while the main FASTQ
 * download array is still working.
 */

include { FETCH_UNMAPPED_BAM } from '../../modules/local/fetch_unmapped_bam'
include { parseBamRuns       } from '../../lib/helpers'

workflow FETCH_BAM_ONLY {

    take:
    strict          // true when the user asked for FETCH_BAM directly

    main:
    todo = parseBamRuns()

    if( strict && params.samples && !todo )
        error "--samples matched no BAM-only run in ${params.bam_manifest}.\n" +
              "  That file holds the 65 runs which have no FASTQ at ENA.\n" +
              "  Selection given: ${params.samples}"

    if( todo )
        FETCH_UNMAPPED_BAM( channel.fromList( todo ) )
    else
        log.info "FETCH_BAM: nothing to do."
}
