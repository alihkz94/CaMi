# Downstream statistics — `--step STATS`

What runs after profiling, why each choice was made, and what the numbers can
and cannot support.

Alpha and beta diversity, PERMANOVA, and the unpaired tumour-vs-healthy contrast
are deliberately **not** built. Each needs a decision about the cohort that
belongs with the person who collected it, not with a default.

---

## Run it

```bash
nextflow run alihkz94/CaMi -profile singularity,slurm \
    --dataRoot /path/to/cohort --input samples.csv --step STATS -resume
```

It is safe to run STATS while profiling is still in progress. It uses the
samples that are finished and quarantines the rest. It is also safe to re-run as
often as you like — nothing it does can affect the profiling run.

---

## The chain

| Process | Reads | Writes |
|---|---|---|
| `AGGREGATE_COUNTS` | every Bracken/Kaiju/Kraken/fastp output | `results/09_tables/` |
| `BUILD_PAIRS` | `sample_metadata.tsv` | `results/10_pairs/` |
| `PAIRED_TEST` | counts + pairs | `results/11_paired_test/` |

`PAIRED_TEST` runs once per view — `bracken/genus`, `bracken/species`,
`kaiju/phylum`, `kaiju/species` — as four independent cached tasks.

---

## Why a separate Python environment

The statistics run in their own image, `cami-stats`, not alongside the profiling
tools. Changing a statistics library must never be able to invalidate a profiling
result, and the reverse.

That separation matters most where there is no container engine. On such a site
the profiling tools come from a conda environment, and installing a scientific
Python stack into it while a 563-sample run depends on it could upgrade a shared
library underneath a running job. `scripts/analysis/setup_stats_env.sh` builds a
venv *from* that interpreter instead, so the Python version matches and the conda
environment is never modified. Point `--stats_python` at it.

Pins are in `env/requirements-stats.txt`, and CI fails if they drift from
`env/stats-conda.yml`.

---

## Aggregation: the integrity gate

Aggregation is where silent corruption happens, so `aggregate_counts.py` refuses
to trust a sample it cannot reconcile.

Kraken 2 and Kaiju both count **pairs**, and on a complete sample all three
independent counts agree exactly:

```
kraken2(classified + unclassified) == kaiju(summary total) == fastp(total_reads)/2
```

Measured on `ASCE17_1940F`: 3,142,213 by all three routes.

A sample that fails this has a half-written file — `publishDir` copies
asynchronously while the pipeline runs — and is **quarantined**, not averaged
into the cohort. The verdict for every sample is in
`09_tables/integrity_report.tsv`.

Three further traps, each guarded:

- **Taxon-name collisions.** Wide matrices pivot on `taxon_id`, never on
  `taxon_name`, so two distinct taxa that share a name cannot be silently
  merged. Names live in `taxa_<method>_<rank>.tsv`.
- **Bracken does not sum to Kraken 2.** Bracken only re-estimates reads it can
  push to the requested rank (49,435 of 54,919 classified pairs on
  `ASCE17_1909F`). Its own `fraction_total_reads` therefore has a per-sample,
  per-rank denominator and is **not comparable between samples**. Every fraction
  here is recomputed against the sample's own non-host pair count.
- **`Homo sapiens` is a database artifact**, not contamination: the standard-16
  database holds a single eukaryote, so cockle reads with no better match land on
  Homo. Flagged as `is_host_artifact`, kept in the tables, excluded from the test
  by default.

### Where the design comes from

The paired test needs four things per sample: which animal it came from, which
tissue, tumour or host, and how the library was prepared. There are two ways to
say so, and they produce the same table.

**A samplesheet**, for any cohort. Add these columns to the `--input` CSV:

```csv
sample,fastq_1,fastq_2,individual,tissue,disease_status,dna_prep
T01,/data/T01_R1.fastq.gz,/data/T01_R2.fastq.gz,animal_01,haemolymph,tumor,native
H01,/data/H01_R1.fastq.gz,/data/H01_R2.fastq.gz,animal_01,foot,matched host,native
```

`disease_status` is one of `tumor`, `matched host`, `unmatched host`, `healthy`.
A value outside that list is **named in the report**. An unrecognised word pairs
nothing, which on its own looks exactly like a cohort that has no tumours, so the
report says which values were not understood rather than leaving you to guess.

`tissue` takes a name or a one-letter code. Both reach tier 1 — `foot` and `F`
mean the same thing, so a cohort that writes words is not silently demoted.

**`individual` and `disease_status` create the design. `tissue` and `dna_prep`
only grade it.** Without the first two there are no pairs, and `BUILD_PAIRS`
says so and names the columns that would build them — a normal outcome, not a
failure. Without the last two the pairs still form and are still tested; the
tier simply records that the prep or the tissue was not stated. An unstated
`dna_prep` is **not** read as WGA.

**The sample-code grammar**, for the cockle study, where the design is encoded in
the sample name. Used only when `--sample_code_scheme bruzos`, which is the
default on the ENA route and nowhere else — another cohort's names must never be
read as tissue codes.

```
EPCE18_851H     ->  individual EPCE18_851, tissue H (haemolymph)
EPCE18_851M     ->  individual EPCE18_851, tissue M (mantle)
FRCE17_840F-wga ->  individual FRCE17_840, tissue F, prep WGA
ENCE17_321H-gp  ->  individual ENCE17_321, tissue H, prep suffix -gp (native)
PACE17_421H1    ->  individual PACE17_421, tissue H, replicate 1
ICCE19_359F_HC  ->  individual ICCE19_359, tissue F, high coverage (60 Gb)
ASCE17_1983     ->  no tissue letter — cannot be paired, reported as such
```

Strip in this order: `_HC`, then the lowercase prep suffix (`-wga`, `-gp`), then
the trailing capital letter is the tissue and an optional digit is the replicate.

Two forms have each caused a silent error. The trailing-digit form dropped the
`PACE17_421` pair. `_HC` is a **high-coverage annotation,
not a tissue** — read naively the code ends in `C`, which is not a tissue in the
Bruzos codebook at all but a country letter, and five healthy **foot** samples
were being recorded as an invented "unlabelled-C" tissue and dropped out of
every foot comparison.

The grammar is unit-tested against every edge case in the cohort, and against the
forms **in combination** — a fix for one suffix can short-circuit another.

#### Tissue cannot be cross-checked, so the grammar is load-bearing

The ENA submission workbook has **no tissue column**. On that route tissue is
derived from the sample code and from nothing else, so a parser regression
produces a wrong arm membership with no error anywhere. That is why the grammar
has more tests than anything else in this step, and why a samplesheet — where a
person states the tissue outright — is the better route when you have the choice.

---

## The paired design

Every animal that contributed **both** a tumour sample and a matched-host sample
becomes one pair. Pairing controls for individual genotype, sampling site, and
batch simultaneously — none of which the unpaired contrast controls for.

The numbers below are the cockle cohort's, as a worked example.

Across all 563 samples: **38 pairs**, graded by how confounded they are.

| Tier | Definition | n | Use |
|---|---|---|---|
| 1 | both native, host tissue = foot | **9** | headline design |
| 2 | both native, host tissue = mantle/adductor | **8** | included by default |
| 3 | both WGA | 5 | no prep confound *within* the pair, but WGA distorts abundance |
| 4 | prep **mismatched** (native vs WGA) | 16 | **excluded** — disease confounded with library prep |

Default is `--max-tier 2`, giving **n = 17**.

---

## The test

1. **CLR, not proportions.** Sequencing measures relative abundance only, so raw
   proportions produce spurious changes whenever an unrelated taxon moves.
   Counts are centred-log-ratio transformed; differences of CLR values are
   log-ratios, which is what a compositional design can support.
2. **Prevalence filter before the transform.** Sparsity is extreme — many
   samples carry under 1,000 microbial pairs. A taxon must be non-zero in
   `--min-prevalence` (default 0.5) of the tested samples.
3. **Pseudocount after filtering** (default 0.5), so it perturbs a small dense
   matrix instead of inventing structure across thousands of empty cells.
4. **Exact paired permutation test.** Under "tumour and host are exchangeable
   within a cockle", each pair's sign is a coin flip. With n pairs there are
   exactly 2^n sign assignments, and at n≈9–17 they are *enumerated* — exact, no
   normality assumption. Above `--exact-max-n` a seeded Monte Carlo sample is
   used and the report says so.
5. **p-value convention matches the mode.** Exhaustive enumeration uses `b/m`
   (the observed vector is always present, so p ≥ 1/2^n > 0). Monte Carlo uses
   the Phipson & Smyth `(b+1)/(m+1)` correction. Applying `+1` to the exact case
   would be conservative and wrong.
6. **Benjamini–Hochberg FDR** across taxa. Interpret `q_bh`, not `p`.

Wilcoxon signed-rank is reported alongside as a rank-based cross-check. It is
secondary — at n < 10 its p-values are granular and it is less powerful.

### The detectability ceiling — read this before interpreting a null result

A paired sign-flip test cannot return a p below **2/2^n**, even if *every* pair
moves the same way. That sets a hard limit on how many taxa can be screened
before nothing can reach q < 0.05:

| pairs | smallest possible p | max taxa where q<0.05 is reachable |
|---|---|---|
| 9 | 0.0039 | ~12 |
| 12 | 0.00049 | ~102 |
| 17 | 1.5e-5 | ~3,200 |

**This is why the default is tier ≤ 2 rather than tier 1 alone.** Going from 9 to
17 pairs moves the ceiling by more than two orders of magnitude — n matters far
more here than pair purity. `paired_test.py` prints this ceiling in every report
and warns when the taxon count exceeds it.

A null result under such a ceiling is **not** evidence of no difference. It is an
underpowered design, and must be reported that way.

### Reading CLR results correctly

CLR is compositional: if a few taxa genuinely rise, everything else *must* fall
in relative terms. In an end-to-end test with a real 8× effect planted in three
taxa, the test recovered all three (12/12 pairs, correct direction) **and** two
unrelated taxa showed significant negative shifts. Those are not false
positives — they are the compositional constraint. Never read a negative CLR
shift as "this taxon decreased" without checking what rose.

---

## The standing caveat

Pairing controls for individual, site, and batch. It does **not** control for
tissue. In the cockle cohort the tumour samples are haemolymph and the host
partners are foot, mantle or adductor, so the two arms differ by tissue as well
as by disease.

Any difference found under such a design is a **tissue** difference as much as a
disease difference, and no amount of pairing separates them. State it explicitly
in any write-up. `pairs.tsv` carries `tumour_tissue` and `host_tissue` for every
pair, so the confound is visible rather than assumed away.

Where you control the sampling, taking both arms from the same tissue removes
this entirely, and is worth more than any statistical correction.

Also on record from the pilot: *Vibrio* is genuinely present but shows **no
tumour enrichment**, and the second-largest bloom in the cohort is a *healthy*
cockle. Treat any "Vibrio causes the cancer" framing as refuted by our own data.
