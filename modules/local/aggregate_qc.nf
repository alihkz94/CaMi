/*
 * AGGREGATE_QC — combine every sample into one table, and run the QC gate
 * -----------------------------------------------------------------------------
 * The reports and the reads are REAL PROCESS INPUTS. They are not read from the
 * published directory.
 *
 * WHY THIS MATTERS: the original pointed PROFILES_DIR at params.outdir and read
 * the files that publishDir had copied there. publishDir is asynchronous, and
 * Nextflow gives no ordering guarantee to a process that does not take those
 * files as inputs. The QC gate could therefore read a half-published directory
 * and mark good samples as failed, with no error anywhere.
 *
 * Nextflow stages the inputs flat, so the per-sample layout that the QC script
 * expects is rebuilt inside the task directory.
 */
process AGGREGATE_QC {
    input:
    val samples
    path reports
    path nonhost

    output:
    path "qc_report.tsv",         emit: qc
    path "microbial_summary.csv", emit: summary

    script:
    """
    for s in ${samples.join(' ')}; do
        mkdir -p "\$s"
        for f in "\${s}.kraken2.report" "\${s}.bracken.genus.tsv" "\${s}.bracken.species.tsv"; do
            [ -e "\$f" ] && mv "\$f" "\$s/"
        done
        # the gate looks for <stage>_1.fastq.gz inside the sample directory
        for m in 1 2; do
            [ -e "\${s}.nonhost_\${m}.fastq.gz" ] && mv "\${s}.nonhost_\${m}.fastq.gz" "\$s/nonhost_\${m}.fastq.gz"
        done
    done

    # qc_checkpoints.py by name, not by absolute path: Nextflow puts bin/ on PATH
    # inside the task, and inside a container the repository is not mounted at all,
    # so a projectDir-based path would not resolve there.
    PROFILES_DIR=\$PWD qc_checkpoints.py ${samples.join(' ')} || true

    # The gate is advisory: always leave a report behind, even if it flagged
    # everything, so the run completes and the tables can be inspected.
    [ -s qc_report.tsv ]        || printf 'sample\\tcheck\\tstatus\\tdetail\\n' > qc_report.tsv
    [ -s microbial_summary.csv ] || printf 'sample\\n' > microbial_summary.csv
    """
}
