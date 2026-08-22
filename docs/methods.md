# Methods, and why the settings are what they are

Read this before you change a threshold. Three inherited defaults were wrong.
Each one changed the bacterial count by more than a factor of 10, and none of them
produced an error message.

---

## 1. The problem this pipeline solves

Study PRJEB58149 sequenced the genome of 563 cockles. The goal was the cockle. But
every sample also holds the DNA of whatever lived in and on the animal. This
pipeline recovers that incidental microbial fraction, and compares cockles that
have disseminated neoplasia with healthy cockles.

The microbial part is small — well under 1% of the reads. That is why the settings
below matter so much: at that level, a filter that is slightly wrong changes the
answer completely.

---

## 2. Host removal: the rule was wrong twice

### Version 0 (inherited): `samtools view -f 12`

Keep a read pair only when BOTH mates are completely unmapped.

This looks safe and is not. `bwa-mem2` reports an alignment for almost any read
that has one seed match, at any quality. A single spurious hit on either mate threw
away the whole pair.

### Version 1: `MAPQ >= 20 AND divergence <= 10%`

Better, and wrong in the opposite direction. **MAPQ < 20 means the aligner cannot
choose between several equally good positions — a REPEAT.** It does not mean the
read is foreign. Molluscan genomes are full of repeats, so this rule labelled a
large block of ordinary cockle sequence "microbial".

### The measurement that settled it

Alignment of 397,334 trimmed pairs of ERR10680552 against the cockle assembly,
794,668 mates:

| category | mates | share |
|---|---|---|
| unmapped | 12,227 | 1.54 % |
| MAPQ ≥ 20, divergence ≤ 10 % | 632,174 | 79.55 % |
| MAPQ ≥ 20, divergence > 10 % | 17,428 | 2.19 % |
| **MAPQ < 20, coverage ≥ 50 %, divergence ≤ 10 %** | **120,559** | **15.17 %** |
| MAPQ < 20, coverage ≥ 50 %, divergence > 10 % | 3,623 | 0.46 % |
| MAPQ < 20, coverage < 50 % | 8,657 | 1.09 % |

The disputed 120,559 mates were classified directly:

| identity | pairs | share |
|---|---|---|
| unclassified (= cockle, absent from the database) | 51,814 | 99.33 % |
| *Homo sapiens* (human repeats) | 320 | 0.61 % |
| bacteria | 26 | 0.05 % |

The 26 bacterial pairs were singletons of 2 to 3 reads spread over unrelated
families — *Halomonas*, *Flavobacterium*, *Clostridium*, *Microcystis*. That is the
signature of a spurious assignment, not of a community.

**Conclusion: the disputed block is host sequence.** Version 1 leaked it into the
results, and it is what made *Homo* the first genus in the QC report.

### Version 2 (current): coverage and identity, no MAPQ

A read is host when it aligns over at least `host_min_cov` (0.50) of its length at
no more than `host_max_div` (0.10) divergence. `host_min_mapq` is 0, so MAPQ has no
effect. The knob stays for sweeps.

The effect on one sample, same reads, same database:

| rule | pairs kept | bacteria |
|---|---|---|
| `-f 12` | 4,138 | 3 |
| MAPQ ≥ 20 | 65,436 | 41 |
| coverage ≥ 50 % | 12,398 | 11 |

### One safety check first

The host filter deletes every read that aligns to the cockle assembly. If the
assembly itself held bacterial contigs — common in molluscs, because the animal is
full of bacteria when it is sequenced — the filter would delete real microbes, and
nothing in the results would show it.

`scripts/calibration/09_check_assembly_contamination.sbatch` classified the
assembly. **No bacterial contigs.** Run it again for any new reference.

---

## 3. The classifier: confidence 0.1 was reporting noise

The inherited pipeline used `--confidence 0.1`. The published pilot numbers were
made at `--confidence 0`. On identical reads:

| database | confidence | classified | bacteria | human |
|---|---|---|---|---|
| standard-16 | 0.1 | 81 | 1 | 61 |
| **standard-16** | **0** | **1,682** | **210** | **1,394** |
| *pilot, published* | | *~1,685* | *~211* | *~1,393* |
| PlusPF | 0.1 | 2,138 | 8 | 1,493 |
| PlusPF | 0.05 | 3,518 | 55 | 2,863 |
| PlusPF | 0 | 7,597 | 1,470 | 5,500 |
| Kaiju nr_euk (protein) | | 6,103 | **3,249** | n/a |

### The false positive floor

More bacteria is not automatically better. To know which of those numbers mean
anything, 2,000,000 fragments of PURE COCKLE DNA were classified. No microbe is
present, so every bacterial call is false by construction.

| setting | bacteria | human | unclassified | false positive rate |
|---|---|---|---|---|
| standard-16 conf 0 | 202 | 6,613 | 1,993,104 | 0.0101 % |
| standard-16 conf 0.05 | 106 | 5,243 | 1,994,344 | 0.0053 % |
| standard-16 conf 0.1 | 19 | 1,953 | 1,997,804 | 0.0009 % |
| PlusPF conf 0 | 2,594 | 20,988 | 1,974,981 | 0.1297 % |
| PlusPF conf 0.05 | 1,336 | 18,201 | 1,977,799 | 0.0668 % |
| PlusPF conf 0.1 | 104 | 13,050 | 1,985,090 | 0.0052 % |

Put the two tables together:

| setting | real reads | host-only floor | signal / noise |
|---|---|---|---|
| **standard-16 conf 0** | 0.0529 % | 0.0101 % | **5.2 ×** |
| PlusPF conf 0 | 0.3700 % | 0.1297 % | 2.9 × |
| PlusPF conf 0.05 | 0.0138 % | 0.0668 % | 0.21 × |
| PlusPF conf 0.1 | 0.0020 % | 0.0052 % | 0.38 × |
| standard-16 conf 0.1 | 0.0003 % | 0.0009 % | 0.33 × |

**Every setting with confidence ≥ 0.05 reports fewer bacteria than pure host DNA
produces by itself.** Those settings do not measure the microbiome. They measure
noise. The inherited `--confidence 0.1` was one of them.

The default is therefore **standard-16 at confidence 0**: the best ratio, and the
only configuration that reproduces the pilot.

PlusPF at confidence 0 recovers more real signal in absolute terms (0.241 % above
its floor, against 0.043 %) but with a worse ratio. Use it as a sensitivity
analysis, not as the primary result:

```
--kraken_db /slurm-databases/Kraken2/PlusPF_20250402 --kraken_mem_gb 130 --kraken_forks 4
```

---

## 4. Why `Homo` appears, and why it is not contamination

The standard-16 database holds exactly ONE eukaryote: human. Cockle is not in it.
A cockle read therefore cannot be named correctly — it can only be unclassified, or
wrong, and when it is wrong the nearest eukaryote is human.

Classifying the whole cockle assembly showed this directly: **99.998 % of the
assembly was called *Homo sapiens*.**

So `Homo` in a report means "host DNA the aligner did not remove". Judge real human
contamination from the alignment rate in `03_human_removed/<sample>.minimap2.human.log`.

---

## 5. Kaiju is not redundant

Kraken 2 matches exact 31-base nucleotide k-mers. A marine bacterium with no close
relative in the database matches nothing and is reported as unclassified. Kaiju
translates the read and searches protein space, which is far more conserved.

On the same 397,334 pairs: Kraken 2 found 210 bacterial pairs, Kaiju found 3,249.
Kaiju finding more than the largest nucleotide database is the expected result for
a marine sample. Use Kraken 2 and Bracken for the abundance backbone, and Kaiju to
show what the nucleotide method missed.

Kaiju is off by default because its index needs about 187 GB of RAM for EACH
concurrent task — it loads the index instead of memory-mapping it. Turn it on with
`--run_kaiju true`.

---

## 6. Read merging is not used

Merging overlapping mates helps amplicon data and assembly. It does not help here:

- Only 24.6 % of pairs overlap. The peak insert is 268 bp and the reads span
  2 × 149 bp, so most pairs have no overlap at all.
- Kraken 2 `--paired` already joins the mates and counts k-mers from both, so
  merging adds no new k-mer.
- Bracken's abundance model takes a FIXED read length (`-r 150`). Merged reads have
  variable length, which breaks that assumption.

`fastp` can merge with `--merge` if a later step needs it — for Kaiju, where a
longer read gives a longer protein to search, or for assembly. No extra tool is
needed.

---

## 7. Depth

At `--subsample 400000` a sample yields about 10 bacterial pairs. Every one is a
singleton, and several belong to organisms that cannot live in a cockle — a 92 °C
hyperthermophile appeared in one test. **A subsample of that size measures nothing.**

`--subsample` is for pipeline tests only. Production uses 0, which keeps every read.
A full run holds roughly 30 million pairs, about 75 times more.
