/*
 * BUILD_PAIRS — derive the within-individual tumour vs matched-host design
 * -----------------------------------------------------------------------------
 * Step 2 of --step STATS. Cheap and deterministic: it reads one table and
 * writes two files, so it costs one CPU and a couple of minutes.
 *
 * Pairs are graded into tiers by how confounded they are (see bin/build_pairs.py).
 * Tier 4 — a native tumour against a WGA host, or the reverse — is written to
 * pairs.tsv for completeness but is excluded from the default test, because
 * inside such a pair disease is confounded with library preparation.
 */
process BUILD_PAIRS {
    tag "tier<=${params.stats_max_tier}"

    input:
    path tables, stageAs: 'tables/*'

    output:
    path "pairs/pairs.tsv",        emit: pairs
    path "pairs/pairs_report.txt", emit: report
    path "pairs/*",                emit: all

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

    build_pairs.py \\
        --tables tables \\
        --outdir pairs
    """
}
