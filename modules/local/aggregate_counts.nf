/*
 * AGGREGATE_COUNTS — every sample's classifier output as one set of tables
 * -----------------------------------------------------------------------------
 * Step 1 of --step STATS.
 *
 * THE REPORTS ARE REAL PROCESS INPUTS, not files read from the published
 * directory. This is the same rule AGGREGATE_QC follows, and for the same
 * reason: publishDir is asynchronous, so a process that does not declare those
 * files as inputs can read a half-published directory and silently aggregate
 * truncated data. Declaring them also makes -resume correct — when a new
 * sample's report appears the input hash changes and this step re-runs.
 *
 * Nextflow stages inputs FLAT, but aggregate_counts.py expects the numbered
 * step layout, so the task rebuilds that layout inside its own work directory.
 * Nothing is copied: symbolic links are enough and cost no disk.
 */
process AGGREGATE_COUNTS {
    tag "${reports.size()} report file(s)"

    input:
    path reports,  stageAs: 'staged/*'
    path manifests, stageAs: 'manifests/*'
    path design

    output:
    path "tables/sample_metadata.tsv",     emit: metadata
    path "tables/counts_long.parquet",     emit: counts
    path "tables/integrity_report.tsv",    emit: integrity
    path "tables/aggregate_report.txt",    emit: report
    // The glob, not a list of names: which tables exist depends on the run.
    // A wide count matrix is written per (method, rank) actually present, and
    // Kaiju's are absent entirely under --run_kaiju false. A named output that
    // is sometimes missing fails the task instead of degrading.
    path "tables/*",                       emit: all

    script:
    // bin/ is already on PATH inside a Nextflow task, so the script is called by
    // name and its `#!/usr/bin/env python3` shebang resolves to whatever python
    // the container, the conda env, or the venv put first.
    //
    // params.stats_python is null everywhere except on a site that provisions the
    // scientific stack by hand — see conf/site/mpi_bremen.config. When it is set,
    // prepending its directory makes the shebang resolve to THAT interpreter
    // without hardcoding a path into the script or touching the environment the
    // profiling tools run from.
    def statsPath = params.stats_python ? "export PATH=\"\$(dirname ${params.stats_python}):\$PATH\"" : ''
    // Unset means "let the script decide": the sample-code grammar is used only
    // on the ENA route, where the codes are known to carry a tissue letter.
    // Passing it explicitly forces the choice either way.
    def schemeArg = params.sample_code_scheme ? "--sample-code-scheme ${params.sample_code_scheme}" : ''
    """
    set -euo pipefail

    # Rebuild the layout the script expects from the flat staging directory.
    mkdir -p results/04_dedup results/05_kraken2 results/06_bracken results/07_kaiju
    for f in staged/*; do
        b=\$(basename "\$f")
        case "\$b" in
            *.dedup.json)      ln -sf "../../\$f" "results/04_dedup/\$b"   ;;
            *.kraken2.report)  ln -sf "../../\$f" "results/05_kraken2/\$b" ;;
            *.bracken.*.tsv)   ln -sf "../../\$f" "results/06_bracken/\$b" ;;
            *.kaiju.*.tsv)     ln -sf "../../\$f" "results/07_kaiju/\$b"   ;;
        esac
    done

    mkdir -p project
    ln -sf ../manifests project/manifests

    ${statsPath}

    aggregate_counts.py \\
        --project project \\
        --results results \\
        --outdir tables \\
        --design ${design} \\
        ${schemeArg} \\
        --min-depth ${params.stats_min_depth}
    """
}
