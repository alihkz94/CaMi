#!/usr/bin/env nextflow
/*
 * =============================================================================
 * CaMi — Cancer–Microbiome interaction
 * -----------------------------------------------------------------------------
 * Recovers the microbial reads that came incidentally with host whole-genome
 * sequencing, and tests that community against cancer status within the same
 * individual. Developed on the transmissible neoplasia of Cerastoderma edule,
 * ENA study PRJEB58149 (563 runs), but the contrast is not specific to cockles.
 *
 * STEPS
 *   --step REFERENCE   download the 2 genomes, build the 3 indexes   (~40 min)
 *   --step FETCH_BAM   the runs that ENA holds only as BAM           (~3-6 h)
 *   --step PROFILE     the analysis (the indexes must already exist)
 *   --step STATS       the paired differential-abundance test
 *   --step all         REFERENCE, FETCH_BAM and PROFILE, in order
 *
 * RUN
 *   nextflow run alihkz94/CaMi -profile singularity,slurm \
 *       --dataRoot /path/to/cohort --step PROFILE -resume
 *
 * Always pass -resume. Nothing here repeats work that finished.
 *
 * DOCUMENTATION
 *   docs/usage.md    how to run it, and what the data directory must contain
 *   docs/output.md   what it writes
 *   docs/methods.md  why the settings are what they are — read this before you
 *                    change any threshold. Two of them were wrong, and each one
 *                    changed the result by more than a factor of 10.
 * =============================================================================
 */

include { CAMI } from './workflows/cami'

workflow {
    CAMI()
}
