/*
 * PREPARE_REFERENCE — download the genomes and build the indexes
 * -----------------------------------------------------------------------------
 * Runs once. About 40 minutes. The results are published to params.refdir and are
 * reused by every later run, so -resume skips this step completely.
 */

include { FETCH_GENOME     } from '../../modules/local/fetch_genome'
include { BUILD_HOST_INDEX } from '../../modules/local/build_host_index'
include { BUILD_MM2_INDEX  } from '../../modules/local/build_mm2_index'

workflow PREPARE_REFERENCE {

    main:
    genomes = FETCH_GENOME( channel.fromList([
        [ 'cockle_V1.2.fa.gz',    params.cockle_url ],
        [ 'GRCh38.primary.fa.gz', params.human_url  ]
    ]) )

    cockle_fa = genomes.filter { name, _f -> name.startsWith('cockle') }.map { _name, f -> f }
    human_fa  = genomes.filter { name, _f -> name.startsWith('GRCh38') }.map { _name, f -> f }

    host_idx = BUILD_HOST_INDEX( cockle_fa )

    // The human genome always, and the vibrio panel only if one was given.
    // params.vibrio_fa is null by default: it is a validation panel used by
    // scripts/calibration/, not by the profiling chain, and a fresh installation
    // has no copy of it. Passing a null through file() would fail the whole
    // REFERENCE step over an index nothing downstream reads.
    to_index = human_fa.map { f -> tuple('GRCh38.primary', f) }
    if( params.vibrio_fa )
        to_index = to_index.mix( channel.of( tuple('vibrio_panel', file(params.vibrio_fa, checkIfExists: true)) ) )

    mm2_idx   = BUILD_MM2_INDEX( to_index )
    human_idx = mm2_idx.filter { f -> f.name.startsWith('GRCh38') }

    emit:
    host  = host_idx.collect()
    human = human_idx.first()
}
