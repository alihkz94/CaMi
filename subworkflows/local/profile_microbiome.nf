/*
 * PROFILE_MICROBIOME — the analysis chain
 * -----------------------------------------------------------------------------
 *   SUBSAMPLE_TRIM -> HOST_REMOVAL -> HUMAN_REMOVAL -> DEDUP
 *                  -> KRAKEN2 + BRACKEN  (abundance backbone)
 *                  -> KAIJU              (protein level, optional but advised)
 *                  -> AGGREGATE_QC
 *
 * Every output is named with the SAMPLE NAME. Where the reads came from — a
 * samplesheet, a directory, or an ENA accession — is decided once in
 * lib/helpers.nf::resolveInputs() and is invisible from here on.
 */

include { SUBSAMPLE_TRIM   } from '../../modules/local/subsample_trim'
include { HOST_REMOVAL     } from '../../modules/local/host_removal'
include { HUMAN_REMOVAL    } from '../../modules/local/human_removal'
include { DEDUP            } from '../../modules/local/dedup'
include { KRAKEN2_BRACKEN  } from '../../modules/local/kraken2_bracken'
include { KAIJU            } from '../../modules/local/kaiju'
include { AGGREGATE_QC     } from '../../modules/local/aggregate_qc'
include { resolveInputs     } from '../../lib/helpers'
include { asBool            } from '../../lib/helpers'

workflow PROFILE_MICROBIOME {

    take:
    host_idx
    human_idx

    main:
    // One list of [sample, read1, read2], whichever route produced it. On the
    // ENA route that is still two disjoint sources — the runs that had FASTQ at
    // the archive, and any BAM-only run that FETCH_BAM has already recovered —
    // merged there, with the same collision check as before.
    rows_ch  = channel.fromList( resolveInputs() )
    reads_ch = rows_ch.map { sample, r1, r2, _acc, _source -> tuple(sample, r1, r2) }

    // Provenance: which sample came from where. On the samplesheet and directory
    // routes there is no accession to record, so the column holds "-" rather
    // than being dropped — the file has one shape whatever the input was.
    rows_ch
        .map { sample, _r1, _r2, acc, source -> "${sample}\t${acc}\t${source}\n" }
        .collectFile( name: 'sample_accession_map.tsv', storeDir: "${params.outdir}/08_summary",
                      seed: "sample_code\trun_accession\tsource\n", sort: true )

    SUBSAMPLE_TRIM( reads_ch )
    HOST_REMOVAL(  SUBSAMPLE_TRIM.out.trimmed, host_idx )
    HUMAN_REMOVAL( HOST_REMOVAL.out.nohost,    human_idx )
    DEDUP(         HUMAN_REMOVAL.out.nonhost )
    KRAKEN2_BRACKEN( DEDUP.out.dedup )

    // Kaiju is ON by default: it found 3,249 bacterial pairs where Kraken 2 found
    // 210 on the same reads, and that difference is the reason the step exists.
    // Turn it off with --run_kaiju false where the nr_euk index will not fit —
    // it needs about 187 GB of RAM for EACH concurrent task.
    //
    // asBool, not params.run_kaiju: from the command line the value is the string
    // "false", which is true in Groovy. See lib/helpers.nf.
    if( asBool(params.run_kaiju) )
        KAIJU( DEDUP.out.dedup )

    // Pass the reports as real inputs so the QC gate cannot race publishDir.
    profiled = KRAKEN2_BRACKEN.out.profiled
    AGGREGATE_QC(
        profiled.map { row -> row[0] }.collect(),
        profiled.map { row -> row[1..3] }.flatten().collect(),
        DEDUP.out.dedup.map { row -> row[1..2] }.flatten().collect()
    )
}
