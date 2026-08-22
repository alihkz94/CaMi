# Containers

Five images cover the whole pipeline. Each is defined by three files that say the
same thing three ways:

| File | Read by | Purpose |
|---|---|---|
| `environment.yml` | both | the pinned software — **the source of truth** |
| `Dockerfile` | Docker / CI | how the published image is built |
| `Singularity.def` | Singularity / Apptainer | how to build it without Docker |

| Image | Tools | Processes |
|---|---|---|
| `cami-fetch` | aria2c, curl, samtools, bgzip, python | `FETCH_GENOME`, `FETCH_UNMAPPED_BAM` |
| `cami-qc` | fastp, seqkit, pigz, python | `SUBSAMPLE_TRIM`, `DEDUP`, `AGGREGATE_QC` |
| `cami-align` | bwa-mem2, minimap2, samtools, gawk | `BUILD_HOST_INDEX`, `BUILD_MM2_INDEX`, `HOST_REMOVAL`, `HUMAN_REMOVAL` |
| `cami-classify` | kraken2, bracken, kaiju | `KRAKEN2_BRACKEN`, `KAIJU` |
| `cami-stats` | python, polars, pyarrow, numpy, scipy, statsmodels | `AGGREGATE_COUNTS`, `BUILD_PAIRS`, `PAIRED_TEST` |

Images are grouped by pipeline stage rather than by tool, because steps like
`HOST_REMOVAL` pipe several tools together in one stream.

## Using them

Nextflow pulls the published images automatically:

```bash
nextflow run alihkz94/CaMi -profile singularity,slurm --dataRoot /path/to/cohort
```

Pre-pull them once on a cold cluster, or on a machine that has a route to the
registry when the cluster does not:

```bash
bash scripts/containers/pull_images.sh /shared/images
nextflow run alihkz94/CaMi -profile singularity,slurm --container_dir /shared/images
```

Build them yourself if you are changing one, or want to verify it:

```bash
bash scripts/containers/build_images.sh singularity          # all five
bash scripts/containers/build_images.sh docker align         # just one
```

## Changing an image

1. Edit `containers/<group>/environment.yml`.
2. Bump `container_tag` in `conf/containers.config` **in the same commit**.
3. Push. `.github/workflows/containers.yml` builds all five, checks that every
   tool actually runs, and publishes them on a `v*` tag.

Never repoint an existing tag at a different image.

## Image digests

Rebuilding from unchanged environment files does not produce identical images:
the base image moves and a conda solve can pick different builds of transitive
dependencies. The pinned versions always hold; the bytes do not. For work whose
numbers will be published, cite the digest rather than the tag.

Record each release's digests here as it is published:

```bash
for g in fetch qc align classify stats; do
    docker buildx imagetools inspect ghcr.io/alihkz94/cami-$g:<tag> \
        | awk '/^Digest:/{print "'"$g"'", $2}'
done
```

### 1.0.0

| Image | Digest |
|---|---|
| `cami-fetch`    | `sha256:3ba26b22bbdf62f8c90429203becb53d83d28267998cedfeaaf1e4f6d3becc41` |
| `cami-qc`       | `sha256:a343d3becba31d2e3950a87f8978a1c8beece0b5c0d196b0eddc5eb2e18ce280` |
| `cami-align`    | `sha256:162cc0d5ef2c8525fe537ed45854bec4a578c2c8de6c7b3e2269a3db47453bf4` |
| `cami-classify` | `sha256:336f1817ff651eec7180ebf9e23f59f4e85b723ff26fb10a9c96d463080d2836` |
| `cami-stats`    | `sha256:abe0973e988a61574e2dfb17f646c8f1e3930db51d6211edda4cbfd773a47e43` |

Each run also records the image it used in `results/pipeline_info/`.

## Databases

The databases are not in the images — Kraken 2 Standard-16 is 16 GB and the Kaiju
`nr_euk` index about 187 GB. They are bind-mounted at run time from `--kraken_db`
and `--kaiju_db`. See [docs/usage.md](../docs/usage.md).
