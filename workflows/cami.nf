/*
 * =============================================================================
 * cami.nf — the top level workflow
 * -----------------------------------------------------------------------------
 * Selects the step and connects the subworkflows.
 *
 * Nextflow 26 uses a strict parser, which does not accept -entry with this
 * layout. The step is selected with --step instead.
 * =============================================================================
 */

include { PREPARE_REFERENCE  } from '../subworkflows/local/prepare_reference'
include { FETCH_BAM_ONLY     } from '../subworkflows/local/fetch_bam_only'
include { PROFILE_MICROBIOME } from '../subworkflows/local/profile_microbiome'
include { DOWNSTREAM_STATS   } from '../subworkflows/local/downstream_stats'
include { checkLocation      } from '../lib/helpers'
include { existingHostIdx    } from '../lib/helpers'

workflow CAMI {

    main:
    checkLocation()

    if( params.step == 'REFERENCE' ) {
        PREPARE_REFERENCE()
    }
    else if( params.step == 'FETCH_BAM' ) {
        FETCH_BAM_ONLY( true )
    }
    else if( params.step == 'PROFILE' ) {
        PROFILE_MICROBIOME(
            channel.value( existingHostIdx() ),
            channel.value( file(params.human_idx, checkIfExists: true) )
        )
    }
    else if( params.step == 'STATS' ) {
        // Reads what PROFILE has already written. Safe to run while PROFILE is
        // still going — it uses the finished samples and quarantines the rest.
        DOWNSTREAM_STATS()
    }
    else if( params.step == 'all' ) {
        ref = PREPARE_REFERENCE()
        FETCH_BAM_ONLY( false )
        PROFILE_MICROBIOME( ref.host, ref.human )
        // STATS is deliberately NOT chained here. It is an analysis step whose
        // settings get revisited many times; running it automatically at the end
        // of a 5-day profiling run would only invite a rushed interpretation.
    }
    else {
        error "Unknown --step '${params.step}'. Use REFERENCE, FETCH_BAM, PROFILE, STATS or all."
    }
}
