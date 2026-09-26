# Changelog

Notable changes to CaMi. Versioning is semantic, and about reproducibility:
**MAJOR** results can change, **MINOR** new capability with the same numbers,
**PATCH** documentation and resources only.

---

## Unreleased

### Changed

- fastp is given the adapter sequences (`--adapter_r1`, `--adapter_r2`; TruSeq by
  default), trims poly-G, drops reads shorter than `--min_read_len` (50 bp) and
  applies its low-complexity filter.
- `docs/methods.md`: adapter trimming (section 8) and the classifier noise floor
  re-measured after trimming with a multi-assembly host reference (section 3).

Results change: adapter dimers and poly-G reads that were assigned to bacteria
are removed. Rerun `SUBSAMPLE_TRIM` and everything after it.

---

## 1.0.1 — 2026-09-02

### Fixed

- A submission rejected by the scheduler is now retried instead of being treated
  as a failed sample. Nextflow retries failed submissions only when the error
  text matches `executor.retry.reason`, which defaults to `Socket timed out`; any
  other wording fell through to the per-process `errorStrategy` and was ignored.
  Because the affected tasks never ran, everything downstream of them waited on
  inputs that could never arrive and the run deadlocked while still holding its
  allocation. `conf/base.config` now sets a broader `reason` pattern covering the
  common transient failures for Slurm, SGE and LSF.
- `errorStrategy` for the per-sample processes now distinguishes a task that ran
  and failed from one that never started. Only the former is eligible to be
  dropped; a task with no exit status always retries.

---

## 1.0.0 — 2026-08-21

First release. Recovers the incidental microbial fraction from host whole-genome
sequencing and tests it against cancer status.

### Pipeline

- Eight profiling steps: adapter and quality trimming (fastp), host removal
  (BWA-MEM2), human removal (minimap2), deduplication, then Kraken 2, Bracken and
  Kaiju, with a per-sample QC gate.
- `--step STATS`: count aggregation with an integrity gate, within-individual
  tumour/host pairing, and an exact paired permutation test with a Wilcoxon
  cross-check and Benjamini–Hochberg correction.

### Input

Three routes, and only the last one downloads anything.

- `--input samples.csv` — a samplesheet: `sample,fastq_1,fastq_2`, plus the
  optional `individual,tissue,disease_status,dna_prep` columns that build the
  paired design for `--step STATS`.
- `--input /path/to/dir` — a directory of paired FASTQ, matched by file name
  (`<name>_R1_001`, `<name>_R1` or `<name>_1`, `.fastq.gz` or `.fq.gz`).
- ENA manifests, for reproducing a published study. See
  [scripts/download/](scripts/download/README.md).

`--sample_code_scheme` decides whether a sample name is parsed for a design.
It is off by default outside the ENA route, so one cohort's naming convention is
never applied to another's.

### Portability

- Five container images on `ghcr.io/alihkz94/cami-*`, plus a conda profile.
- Executors `local`, `slurm`, `sge`, `lsf`; engines `singularity`, `apptainer`,
  `docker`, `podman`.
- `--container_dir` to run from prebuilt images with no registry access.
- Site settings live in `conf/site/`; the defaults are portable.

### Calibration

Three defaults were measured rather than inherited, and each changed the
bacterial count by more than tenfold. Full evidence in
[docs/methods.md](docs/methods.md).

- Host reads are filtered on alignment coverage and identity, never on MAPQ.
- Kraken 2 runs at `--confidence 0`, chosen against a measured false-positive
  floor from pure host DNA.
- `Homo sapiens` calls are a database artifact, not contamination.

### Testing

83 unit tests covering the statistics against brute-force ground truth and the
sample-name parser against every edge case. CI lints every Nextflow file,
resolves every profile, checks that the two environment files agree, and
resolves all three input routes against real files.

### Images

`ghcr.io/alihkz94/cami-{fetch,qc,align,classify,stats}:1.0.0`.
Digests: [containers/README.md](containers/README.md).
