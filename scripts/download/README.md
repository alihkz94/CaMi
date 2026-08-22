# Getting reads from ENA

**Most people do not need anything in this directory.**

CaMi takes reads by three routes, and only the third involves a download:

| You have | Use |
|---|---|
| FASTQ on disk | `--input samples.csv` or `--input /path/to/fastq` |
| FASTQ from a sequencing provider | the same — point `--input` at the folder |
| only an ENA study accession | the scripts here, once, then `--input` |

Reads that are already on disk go straight into the pipeline. Nothing here runs,
no manifest is built, and no accession is needed. See
[docs/usage.md](../../docs/usage.md).

These scripts exist because CaMi was developed on ENA study **PRJEB58149**, where
563 runs had to be fetched and checked before anything could start. They are kept
so that study stays reproducible, and because fetching a different ENA study is a
reasonable thing to want.

## What each script does

| Script | Purpose |
|---|---|
| `00_build_manifest.py` | study metadata → `manifests/download_list.tsv` and `skipped_bam_only.tsv` |
| `01_download_array.sbatch` | a Slurm array that fetches the FASTQ and verifies each MD5 |
| `03_status.sh` | how much has arrived, and what failed |
| `05_chain_after_download.sbatch` | start profiling as soon as the download finishes |
| `06_fetch_bam_isolated.sbatch` | the runs ENA holds only as BAM, outside the main pipeline |
| `build_ena_metadata.py` | a submission workbook → per-sample metadata for `--step STATS` |

**`manifests/sample_name_map.csv` is the file this route actually needs**, and
nothing above writes it. It maps `run_accession,sample_code` and is produced by
`scripts/maintenance/02_name_by_sample.py`, which was written for one cluster and
carries that cluster's paths and user name near the top. Read it before running
it. The file itself is two columns and can equally be written by hand.

This is the seam in the ENA route, and it is the reason the samplesheet exists:
three scripts, a checksum convention and a hand-edited path all have to line up
before the pipeline sees a single read.

## The checksum marker

`01_download_array.sbatch` writes `<ACC>.md5.ok` beside each verified run, and the
pipeline processes only runs that have one. A file of the right size can still be
truncated, and a truncated FASTQ produces a plausible, wrong answer.

That check belongs to this route alone. It is not applied to `--input`, where the
reads did not come through this directory and CaMi has no checksum to compare
against. Use `--require_md5ok false` if you took the ENA route and want to
process runs that have not been verified.

## Metadata for the statistics

`--step STATS` needs to know which samples came from the same animal, and which
are tumour and which are host. Two ways to say so:

- **A samplesheet**, for any cohort — add `individual,tissue,disease_status,dna_prep`
  to the `--input` CSV. This is the general route.
- **`build_ena_metadata.py`**, for a study whose submission workbook you hold. It
  reads the workbook once into `manifests/ena_metadata.tsv`.

The second refuses to write a table that disagrees with PRJEB58149's known
tallies, because that workbook has **no tissue column** — tissue comes from the
sample code and nothing else can cross-check it. For a different ENA study those
tallies cannot match; pass `--skip-checks` and read the report it prints.

It needs `openpyxl`, which is deliberately not part of the pipeline environment:

```bash
pip install openpyxl==3.1.5
python3 scripts/download/build_ena_metadata.py \
    --workbook 'ENA Submission.xlsx' \
    --out <dataRoot>/manifests/ena_metadata.tsv
```
