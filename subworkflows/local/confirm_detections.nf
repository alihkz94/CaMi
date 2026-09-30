/*
 * CONFIRM_DETECTIONS — every Kraken 2 detection checked by two independent methods
 * -----------------------------------------------------------------------------
 *   CONFIRM_KRAKEN      per-read Kraken 2 calls, report must be byte-identical
 *   CONFIRM_DRAW        detections (units), clade reads, <= 20 reads each for BLAST
 *   CONFIRM_BLAST_PREP  one pooled query set in chunks
 *   CONFIRM_BLAST       megablast against nt, one task per chunk
 *   CONFIRM_COMBINE     Kaiju and BLAST verdicts per unit -> triple/double/single
 *
 * Kaiju's per-read output comes from the KAIJU step itself. A sample that is
 * missing any piece is not silently absent: it is listed in
 * 12_confirmation/confirm_incomplete.txt.
 */

include { CONFIRM_KRAKEN     } from '../../modules/local/confirm_kraken'
include { CONFIRM_DRAW       } from '../../modules/local/confirm_draw'
include { CONFIRM_BLAST_PREP } from '../../modules/local/confirm_blast_prep'
include { CONFIRM_BLAST      } from '../../modules/local/confirm_blast'
include { CONFIRM_COMBINE    } from '../../modules/local/confirm_combine'

workflow CONFIRM_DETECTIONS {

    take:
    reads     // [sample, read1, read2]
    reports   // [sample, kraken2.report]
    kaiju     // [sample, kaiju.out.gz]

    main:
    CONFIRM_KRAKEN( reads.join(reports) )
    calls = CONFIRM_KRAKEN.out.calls

    CONFIRM_DRAW(
        reads.map { sample, r1, _r2 -> tuple(sample, r1) }
             .join(reports)
             .join(calls.map { sample, c, _pairs -> tuple(sample, c) })
    )

    CONFIRM_BLAST_PREP(
        CONFIRM_DRAW.out.draw.map { _s, _u, _k, draw, fa -> [draw, fa] }.flatten().collect()
    )
    CONFIRM_BLAST( CONFIRM_BLAST_PREP.out.chunks.flatten() )

    // [sample, report, units, k2reads, draw, fasta, calls, pairs, kaiju]
    complete = reports.join(CONFIRM_DRAW.out.draw).join(calls).join(kaiju)

    reports.map { row -> [row[0]] }
        .join( complete.map { row -> [row[0], true] }, remainder: true )
        .filter { row -> row[1] == null }
        .map { row -> "${row[0]}\n" }
        .collectFile( name: 'confirm_incomplete.txt', storeDir: "${params.outdir}/12_confirmation", sort: true )

    CONFIRM_COMBINE(
        complete.map { row -> [row[1], row[2], row[3], row[4], row[7], row[8]] }.flatten().collect(),
        CONFIRM_BLAST_PREP.out.queries,
        CONFIRM_BLAST.out.hits.collect().ifEmpty([]),
        complete.map { row -> row[0] }.collectFile( name: 'confirm_samples.txt', newLine: true, sort: true ),
        channel.value( file(params.confirm_taxdump, checkIfExists: true) )
    )
}
