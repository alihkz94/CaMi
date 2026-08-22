# Usage

## Requirements

- Nextflow 24.04 or newer
- A container engine: Singularity, Apptainer, Docker or Podman (or use `-profile conda`)
- Disk: about 3× your raw FASTQ, for `work/`
- A Kraken 2 database (~16 GB for Standard-16)
- A Kaiju index, optional (~187 GB, and the same again in RAM per task)

## The reads

CaMi takes reads by three routes. **Only the third downloads anything.** If the
FASTQ are already on disk — from a sequencing provider, or from a download that
happened months ago — use route 1 or 2 and ignore `scripts/download/` entirely.

### 1. A samplesheet

```bash
nextflow run alihkz94/CaMi -profile singularity --input samples.csv --dataRoot /path/to/cohort ...
```

```csv
sample,fastq_1,fastq_2,individual,tissue,disease_status,dna_prep
T01,/data/T01_R1.fastq.gz,/data/T01_R2.fastq.gz,animal_01,haemolymph,tumor,native
H01,/data/H01_R1.fastq.gz,/data/H01_R2.fastq.gz,animal_01,foot,matched host,native
```

| Column | Required | What it is |
|---|---|---|
| `sample` | yes | names every output file. Letters, digits, `.`, `-`, `_` |
| `fastq_1`, `fastq_2` | yes | paired reads. Absolute, or relative to the samplesheet |
| `individual` | for `--step STATS` | the same value for both samples from one animal |
| `tissue` | for `--step STATS` | a name (`foot`) or a one-letter code (`F`) |
| `disease_status` | for `--step STATS` | `tumor`, `matched host`, `unmatched host` or `healthy` |
| `dna_prep` | for `--step STATS` | `native` or `wga` |

The four optional columns build the paired tumour-vs-host design. Leave them out
and steps 1–8 run exactly the same; only the paired test is skipped, and it says
so rather than reporting an empty result.

### 2. A directory of FASTQ

```bash
--input /path/to/fastq
```

Files are paired by name: `<name>_R1_001`, `<name>_R1` or `<name>_1`, with
`.fastq.gz` or `.fq.gz`. The name before the mate suffix becomes the sample name.
A read 1 with no matching read 2 stops the run — half a pair on disk means an
interrupted copy, and finding out later is worse.

This route states no design, so `--step STATS` runs the aggregation and reports
that it has no pairs to test. Add a samplesheet when you want the paired test.

### 3. An ENA study

The route CaMi was developed on, and the default when `--input` is not given.
It reads `manifests/sample_name_map.csv` and expects accession-named files in
`fastq/`. See [scripts/download/README.md](../scripts/download/README.md).

A run is processed only once `<ACC>.md5.ok` exists beside it, because a file of
the right size can still be truncated. Use `--require_md5ok false` to override.
That check applies to this route alone — on routes 1 and 2 there is no checksum
to compare against.

## The data directory

`--dataRoot` points at your data. Code and data are kept separate.

```
<dataRoot>/
├── reference/                written by --step REFERENCE
├── kraken2_db/               or point --kraken_db elsewhere
├── work/                     Nextflow work directory — do not delete
└── results/                  output
```

Everything above is created. On route 3 the directory also holds `fastq/`,
`fastq_from_bam/` and `manifests/`; on routes 1 and 2 the reads may live
anywhere, and `--dataRoot` only decides where the work and results go.

## Databases

**Kraken 2** — prebuilt indexes at <https://benlangmead.github.io/aws-indexes/k2>:

```bash
mkdir -p <dataRoot>/kraken2_db/standard_16gb && cd $_
curl -O https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20250402.tar.gz
tar xzf k2_standard_16gb_20250402.tar.gz
```

Memory and parallelism follow the database size:

| Database | `--kraken_mem_gb` | `--kraken_forks` |
|---|---|---|
| Standard-16 (16 GB) | 32 | 8 |
| PlusPF (94 GB) | 130 | 4 |

**Kaiju** — get an index from
<https://bioinformatics-centre.github.io/kaiju/downloads.html> or build one with
`kaiju-makedb -s nr_euk`, then:

```bash
--kaiju_db /db/Kaiju/kaiju_db_nr_euk.fmi \
--kaiju_nodes /db/Kaiju/nodes.dmp \
--kaiju_names /db/Kaiju/names.dmp
```

Use `--run_kaiju false` if the memory is not available. Note that on the same
reads Kaiju found 3,249 bacterial pairs where Kraken 2 found 210.

## Running

```bash
# 1. genomes and indexes (~40 min, once)
nextflow run alihkz94/CaMi -profile singularity,slurm \
    --dataRoot /path/to/cohort --step REFERENCE -resume

# 2. the analysis
nextflow run alihkz94/CaMi -profile singularity,slurm \
    --dataRoot /path/to/cohort --input samples.csv --step PROFILE \
    --kaiju_db /db/Kaiju/kaiju_db_nr_euk.fmi \
    --kaiju_nodes /db/Kaiju/nodes.dmp \
    --kaiju_names /db/Kaiju/names.dmp -resume

# 3. the statistics
nextflow run alihkz94/CaMi -profile singularity,slurm \
    --dataRoot /path/to/cohort --input samples.csv --step STATS -resume
```

`--step FETCH_BAM` sits between 1 and 2 and applies to route 3 only: it recovers
archive runs that ENA holds as BAM rather than FASTQ. With no such manifest it
reports that there is nothing to fetch and moves on.

`--step all` runs REFERENCE, FETCH_BAM and PROFILE in order. `STATS` is separate
because it is usually re-run several times; it can also run while `PROFILE` is
still going, using whatever samples have finished.

Always pass `-resume`. Do not delete `work/` — `-resume` needs it.

## Options

| Option | What it does |
|---|---|
| `--input PATH` | a samplesheet, or a directory of FASTQ. Unset = the ENA route |
| `--sample_code_scheme` | `bruzos` reads the cockle study's sample-code grammar; `none` treats a name as a name. Default: `bruzos` on the ENA route, `none` otherwise |
| `--samples A,B` | only these samples; also accepts a file, one name per line |
| `--subsample N` | use N read pairs — **testing only** |
| `--run_kaiju false` | skip the protein search |
| `--publish_intermediates false` | do not publish trimmed and host-removed FASTQ |
| `--kraken_db PATH` | another database; move `--kraken_mem_gb`/`--kraken_forks` with it |
| `--require_md5ok false` | process runs without a verified checksum |
| `--container_dir PATH` | run from prebuilt `.sif` files instead of pulling |
| `--stats_max_tier N` | how confounded a tumour/host pair may be (see statistics.md) |

Booleans need an explicit value: `--run_kaiju false`, not `--run_kaiju`.

`--subsample` is for checking that the pipeline runs. At 400,000 pairs a sample
yields about ten bacterial read pairs, which measures nothing.

## No registry access

Fetch the images somewhere with a connection, copy them over, and run offline:

```bash
bash scripts/containers/pull_images.sh /tmp/cami-images
rsync -a /tmp/cami-images/ cluster:/shared/cami-images/

nextflow run /path/to/CaMi -profile singularity,slurm \
    --container_dir /shared/cami-images --dataRoot /path/to/cohort --step PROFILE
```

Behind a proxy, Singularity needs it in its own environment:
`export HTTPS_PROXY=$https_proxy HTTP_PROXY=$http_proxy`.

## When something fails

A failed sample does not stop the run: per-sample steps retry three times with
more memory, then that sample is dropped and the rest continue. It appears in the
Nextflow summary and in `trace.txt`, and is absent from `qc_report.tsv`.

```bash
grep -E 'FAILED|ignored' <dataRoot>/results/pipeline_info/trace.txt
cat <dataRoot>/results/08_summary/qc_report.tsv
```

`trace.txt` records `peak_rss` and `realtime` per task — use it to correct the
resource requests in `conf/base.config`.

The steps that combine all samples stop the run on failure rather than continuing
with partial data.

## Adding your cluster

`conf/site/` holds one file per cluster. Copy `conf/site/mpi_bremen.config`, add
a profile line to `nextflow.config`, and open a pull request.

## The MPI Bremen installation

The data sits on the local disk of one node, `slurm-406`, which has two paths:
`/scratch/ahakimzadeh/...` on the node and `/mnt/slurm-406/ahakimzadeh/...` from
the login node. The head process must run on the node so it and its tasks use the
same path, and because `-resume` keys on that path. `run_pipeline.sh` handles it:

```bash
cd /scratch/ahakimzadeh/CaMi
./run_pipeline.sh --step PROFILE
tail -f ../Alicia_nature_data/slurm_logs/nf_<jobid>.out
```

This cluster has no container engine, so `-profile mpi_bremen` uses a conda
environment instead. `TMPDIR` must stay on `/scratch`.
