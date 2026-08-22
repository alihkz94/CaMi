# Output

The pipeline writes one directory for each step. The steps have numbers, so the
order of the analysis is clear. Each directory holds the files of ALL samples for
that step. Every file starts with the sample code.

```
results/
├── 01_trimmed/          adapter and quality trimming (fastp)
├── 02_host_removed/     cockle reads removed (BWA-MEM2)
├── 03_human_removed/    human reads removed (minimap2)
├── 04_dedup/            duplicates removed — THE MICROBIAL READ SET
├── 05_kraken2/          taxonomic assignment (Kraken 2)
├── 06_bracken/          abundance estimates (Bracken)
├── 07_kaiju/            protein level assignment — ONLY with --run_kaiju
├── 08_summary/          the tables that combine all samples
├── 09_tables/           --step STATS: tidy counts, metadata, integrity report
├── 10_pairs/            --step STATS: the tumour/matched-host pairs and tiers
├── 11_paired_test/      --step STATS: the differential abundance results
└── pipeline_info/       run reports from Nextflow
```

Directory `07_kaiju` is absent when you do not use `--run_kaiju`. The number stays
reserved, so the same step always has the same number. Directories 09 to 11 appear
only after `--step STATS`.

## Sample names

All output files use the **sample name** — the `sample` column of the
samplesheet, or, on the ENA route, the sample code such as `EICE18_906H-wga`.

`08_summary/sample_accession_map.tsv` connects each sample to where it came from.
On the ENA route that is the run accession; on the samplesheet and directory
routes there is no accession and the column holds `-`.

On the ENA route the raw files in `fastq/` keep their accession names. That
directory is the archive copy, and its names must stay the same.

---

## 01_trimmed

| File | Content |
|---|---|
| `<sample>.fastp.json` | Trimming statistics, and the insert size histogram |
| `<sample>.fastp.html` | The same statistics as a page you can read |
| `<sample>.fastp.log` | Messages from the program |
| `<sample>.trim_1.fastq.gz` | Trimmed reads — **only with `--publish_intermediates`** |

The trimmed reads are almost as large as the raw input. To publish them for 498
samples costs about 12 TB, and the next step reads them from the work directory
anyway. Use `--publish_intermediates` when you must inspect one sample.

## 02_host_removed

| File | Content |
|---|---|
| `<sample>.bwa-mem2.cockle.log` | Alignment messages |
| `<sample>.nohost_1.fastq.gz` | Reads after cockle removal — **only with `--publish_intermediates`** |

A read is host when it aligns over at least `--host_min_cov` of its length at no
more than `--host_max_div` divergence. See `docs/methods.md`: the rule went
through two wrong versions before this one, and the wrong versions changed the
result by a factor of 14.

## 03_human_removed

Same pattern as step 02, with minimap2 against GRCh38.

## 04_dedup

**This is the microbial read set. It is the useful product of the pipeline.**

| File | Content |
|---|---|
| `<sample>.nonhost_1.fastq.gz` | Microbial reads, forward |
| `<sample>.nonhost_2.fastq.gz` | Microbial reads, reverse |
| `<sample>.dedup.json` | Duplicate statistics |
| `<sample>.dedup.log` | Messages from the program |

These files are always published. After host and human removal the set is a few
percent of the raw data, so the cost is small. Use these files for assembly, for
a protein level search with Kaiju, or for a targeted BLAST search.

## 05_kraken2

| File | Content |
|---|---|
| `<sample>.kraken2.report` | The full taxonomy tree with read counts |
| `<sample>.kraken2.log` | Messages from the program |

## 06_bracken

| File | Content |
|---|---|
| `<sample>.bracken.genus.tsv` | Abundance for each genus |
| `<sample>.bracken.species.tsv` | Abundance for each species |
| `<sample>.bracken.*.log` | Messages from the program |

Bracken can fail when a report has almost no classified reads. The pipeline then
writes an empty table and continues. The QC report shows the sample as low power.

## 07_kaiju

Written only with `--run_kaiju`.

| File | Content |
|---|---|
| `<sample>.kaiju.summary.tsv` | Read counts for each kingdom |
| `<sample>.kaiju.phylum.tsv` | Read counts for each phylum |
| `<sample>.kaiju.species.tsv` | Read counts for each species |

Kaiju searches protein sequence, so it finds organisms that have no close relative
in any nucleotide database. On the same reads it found 3,249 bacterial pairs where
Kraken 2 found 210. Use Kraken 2 and Bracken for the abundance backbone, and Kaiju
to show what the nucleotide method missed.

## 08_summary

| File | Content |
|---|---|
| `qc_report.tsv` | One row for each check, for each sample |
| `microbial_summary.csv` | One row for each sample, with the read counts at each step |
| `sample_accession_map.tsv` | Sample code and ENA accession, side by side |

### How to read the QC report

| Check | Meaning |
|---|---|
| `microbial_count` | WARN below 100 pairs, FAIL below 10. The sample has low power. |
| `dominant_taxon` | WARN when the first genus is not a marine organism. |

A WARN does not stop the pipeline. It marks a sample for you to look at.

### Important: `dominant_taxon = Homo` is expected, and it is NOT contamination

The standard-16 database holds exactly ONE eukaryote: human. Cockle is not in it,
so a cockle read can never be named correctly — it can only be unclassified, or
wrong. When the whole cockle assembly was classified, 99.998% of it was called
*Homo sapiens*.

So `Homo` in these reports means "host DNA that the aligner did not remove". It
does not mean a person touched the sample. Judge human contamination from the
alignment rate in `03_human_removed/<sample>.minimap2.human.log`, not from this
check.

## pipeline_info

Nextflow writes `report.html`, `timeline.html`, `trace.txt`, and `dag.html`.
Use `trace.txt` to find the memory and time that each task used. Use it to
correct the resource requests before the next large run.
