/*
 * =============================================================================
 * DOWNSTREAM_STATS — the statistics that consume the profiling output
 * -----------------------------------------------------------------------------
 * Run with:  nextflow run alihkz94/CaMi --step STATS ... -resume
 *
 * WHY THIS IS A SEPARATE STEP
 *   It reads what PROFILE has already written and adds nothing to the profiling
 *   run itself. Keeping it separate means it can be re-run, re-tuned, and
 *   re-interpreted as often as the analysis needs, without ever putting the
 *   563-sample profiling run at risk. It is also safe to run while PROFILE is
 *   still going: it simply uses whatever samples are finished, and the
 *   integrity gate quarantines anything half-written.
 *
 * WHY IT IS IN NEXTFLOW AT ALL
 *   Crash safety and provenance. Each stage is a cached task, so -resume picks
 *   up exactly where a stopped run left off, and every input hash, command and
 *   output is recorded in the same trace as the rest of the pipeline.
 *
 * THE CHAIN
 *   AGGREGATE_COUNTS  every classifier output -> one set of tidy tables
 *   BUILD_PAIRS       -> the within-individual tumour vs matched-host design
 *   PAIRED_TEST       -> exact paired test, once per (method, rank)
 *
 * NOT INCLUDED ON PURPOSE
 *   Alpha/beta diversity, PERMANOVA and the unpaired tumour-vs-healthy contrast
 *   are deliberately absent. Each needs a decision about the cohort that belongs
 *   with the person who collected it — above all whether healthy samples of a
 *   DIFFERENT tissue may enter a disease comparison at all. Answering that with
 *   a default would produce a number nobody should trust.
 * =============================================================================
 */

include { AGGREGATE_COUNTS } from '../../modules/local/aggregate_counts'
include { BUILD_PAIRS      } from '../../modules/local/build_pairs'
include { PAIRED_TEST      } from '../../modules/local/paired_test'
include { resolveDesign    } from '../../lib/helpers'
include { designColumns    } from '../../lib/helpers'

workflow DOWNSTREAM_STATS {

    main:
    // Every per-sample file the aggregation reads becomes a real input, so the
    // resume cache invalidates correctly when more samples finish.
    reports = channel.empty()
        .mix( channel.fromPath("${params.outdir}/04_dedup/*.dedup.json") )
        .mix( channel.fromPath("${params.outdir}/05_kraken2/*.kraken2.report") )
        .mix( channel.fromPath("${params.outdir}/06_bracken/*.bracken.*.tsv") )
        .mix( channel.fromPath("${params.outdir}/07_kaiju/*.kaiju.*.tsv") )
        .collect()

    // ifEmpty([]) is load-bearing, not defensive. collect() on an empty channel
    // emits NOTHING, so a cohort with no manifests directory — every cohort that
    // did not come from ENA — left AGGREGATE_COUNTS with an input that never
    // arrived, and --step STATS finished in seconds having run no task at all.
    manifests = channel.fromPath("${params.projectRoot}/manifests/*").collect().ifEmpty([])

    // The design, when a samplesheet supplied one. One channel element holding
    // the whole file: a header alone means "no samplesheet design", and
    // aggregate_counts.py then reads the ENA manifests exactly as before.
    //
    // TAB-separated, not comma. readSamplesheet() accepts either separator, so a
    // tab-separated sheet may hold a value with a comma in it; writing this file
    // with commas then produced a row with more fields than the header, and the
    // parse error named this file rather than the sheet that caused it.
    // readSamplesheet() rejects both characters in a value, so a tab here is safe.
    def designRows = resolveDesign()
    def designText = (['sample'] + designColumns()).join('\t') + '\n' +
                     (designRows ?: []).collect { r -> r.join('\t') }.join('\n') +
                     (designRows ? '\n' : '')
    design = channel.of( designText ).collectFile( name: 'design.tsv' )

    AGGREGATE_COUNTS( reports, manifests, design )
    BUILD_PAIRS( AGGREGATE_COUNTS.out.all.collect() )

    // One task per view of the data. Bracken gives the abundance backbone;
    // Kaiju is the independent protein-level check, and the two disagreeing is
    // informative rather than a fault.
    views = channel.fromList([
        ['bracken', 'genus'  ],
        ['bracken', 'species'],
        ['kaiju',   'phylum' ],
        ['kaiju',   'species'],
    ])

    PAIRED_TEST(
        views,
        AGGREGATE_COUNTS.out.all.collect(),
        BUILD_PAIRS.out.all.collect()
    )

    emit:
    tables = AGGREGATE_COUNTS.out.all
    pairs  = BUILD_PAIRS.out.pairs
    paired = PAIRED_TEST.out.all
}
