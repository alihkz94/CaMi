/*
 * PAIRED_TEST — exact paired tumour vs matched-host differential abundance
 * -----------------------------------------------------------------------------
 * Step 3 of --step STATS. Runs once per (method, rank) combination, so the
 * Bracken and Kaiju views are independent tasks: one can fail or be re-run
 * without touching the other, and -resume caches them separately.
 *
 * COST: the exact test enumerates 2^n sign vectors, where n is the number of
 * pairs. At n=17 that is 131,072 vectors per taxon — seconds. The enumeration
 * is chunked, so memory stays flat as n grows, and above params.stats_exact_max_n
 * the script switches to a seeded Monte Carlo sample and says so in its report.
 */
process PAIRED_TEST {
    tag "${method}/${rank}"

    input:
    tuple val(method), val(rank)
    path tables, stageAs: 'tables/*'
    path pairs,  stageAs: 'pairs/*'

    output:
    path "paired/*", emit: all
    path "paired/paired_test_report.txt", emit: report

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
    """
    set -euo pipefail
    ${statsPath}

    paired_test.py \\
        --tables tables \\
        --pairs pairs \\
        --outdir paired \\
        --method ${method} \\
        --rank ${rank} \\
        --max-tier ${params.stats_max_tier} \\
        --min-prevalence ${params.stats_min_prevalence} \\
        --pseudocount ${params.stats_pseudocount} \\
        --exact-max-n ${params.stats_exact_max_n} \\
        --n-perm ${params.stats_n_perm} \\
        --seed ${params.seed}
    """
}
