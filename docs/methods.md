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

**With one host assembly and untrimmed adapters, every setting with confidence
≥ 0.05 reports fewer bacteria than pure host DNA produces by itself.** Those settings do not measure the microbiome. They measure
noise. The inherited `--confidence 0.1` was one of them.

The default is therefore **standard-16 at confidence 0**: the best ratio, and the
only configuration that reproduces the pilot.

PlusPF at confidence 0 recovers more real signal in absolute terms (0.241 % above
its floor, against 0.043 %) but with a worse ratio. Use it as a sensitivity
analysis, not as the primary result:

```
--kraken_db /slurm-databases/Kraken2/PlusPF_20250402 --kraken_mem_gb 130 --kraken_forks 4
```

### Re-measured with adapter and poly-G trimming and a stronger host reference

The floor above was measured before adapter and poly-G trimming (section 8) and
with a single host assembly. The measurement was repeated on the full cockle
cohort (563 samples) after that trimming and a second host screen against three *C. edule*
assemblies (GCA_947846245.1, GCA_963989375.1, GCA_963989325.1). The host-only
floor was taken from 2.2 million real cockle reads that the single-assembly screen
had missed, and from 2 million simulated reads.

| standard-16 | host-only floor | median sample / floor | samples > 5 × |
|---|---|---|---|
| conf 0 | 0.0901 % | 6 × | 59 % |
| conf 0.05 | 0.0086 % | 31 × | 93 % |
| conf 0.1 | 0.0020 % | 81 × | 98 % |

With that host reference, confidence 0.05 is the better setting: it removes the
genera produced by host reads while keeping validated marine genera. Confidence
0.1 starts to lose real genera. The default stays at 0 because the ratio depends
on how well the host reference covers the host's diversity; measure the floor for
your own host before changing it.

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

Kaiju is on by default. Its index needs about 187 GB of RAM for EACH concurrent
task — it loads the index instead of memory-mapping it — so turn it off with
`--run_kaiju false` on a machine that cannot spare the memory.

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

---

## 8. Adapter and poly-G trimming

Up to version 1.0.1, fastp found adapters by the overlap of the two mates, which
needs an insert of at least about 30 bp. Adapter dimers (inserts of 10–30 bp)
therefore kept their adapter. Poly-G tails were not trimmed either, because fastp
turns poly-G trimming on only when it recognises a two-colour instrument from the
Illumina read headers, and ENA-renamed reads no longer carry them. Kraken 2 assigned
these reads to bacteria: in the cockle cohort they became two of the most prevalent
"species" (*Xanthomonas euvesicatoria* and *Mycobacterium canetti*). None of 420
such reads had a BLAST hit to either genus.

Since 1.1.0, fastp is given the adapters (`--adapter_r1`, `--adapter_r2`; TruSeq by
default), always trims poly-G (`--trim_poly_g`), drops reads shorter than
`--min_read_len` (50 bp) and applies its low-complexity filter. Adapter-carrying reads were about 0.1 % of the non-host
reads, but they sat in a few taxa and dominated those.

---

## 9. Triple confirmation of detections (`--run_confirm`)

A Kraken 2 call rests on shared 31-mers. When the true organism has no close
relative in the database, its reads can still share k-mers with a relative that
is present, and the call names the wrong genus. Low-biomass samples make this
worse, because a few hundred misassigned reads are then a large share of the
community. `--run_confirm` checks every detection with two methods that do not
share the detecting method's database or its algorithm.

Kaiju searches protein and finds organisms that Kraken 2 misses. Its own
detections are therefore units too, checked the same way, instead of being
discarded because Kraken 2 did not make them.

### What is checked

A **unit** is one (sample, taxon) pair that Kraken 2 or Kaiju detects: a genus or
species with at least `--confirm_min_reads` reads (10) and at least
`--confirm_min_frac` (1 %) of the sample's genus-resolved reads, each method in its
own counts. The same denominator is used at both ranks. Human (9606, 9605) is
excluded; for Kaiju only Bacteria, Archaea and Viruses count, because its
database also holds eukaryotes and unremoved host reads land there.
`--confirm_merge_family` counts the genera of one family as a single entry, for a
family the database cannot resolve.

A taxon both methods detect in a sample is one unit (`detected_by` =
`kraken+kaiju`) and is judged as a Kraken 2 unit. A taxon only Kaiju detects is a
Kaiju-only unit (`detected_by` = `kaiju`).

| leg | method | database | Kraken 2 unit | Kaiju-only unit |
|---|---|---|---|---|
| Kraken 2 | nucleotide k-mers | the profiling database | the detection | its calls on Kaiju's reads |
| Kaiju | translated protein search | nr_euk | its calls on the clade's reads | the detection |
| BLAST | megablast alignment | nt | ≤ 20 of the clade's reads | ≤ 20 of Kaiju's reads |

The Kraken 2 leg is rerun with `--output` to get the per-read calls, and the rerun
report must be byte-identical to the profiled one. The per-read calls inside each
unit's clade must add up to the report's clade count. The BLAST reads are drawn
with a fixed seed, so a rerun searches the same reads.

All taxids (Kraken 2's database, Kaiju's `nodes.dmp`, nt's taxdb) are resolved
through one current NCBI taxdump, `merged.dmp` included, and compared by lineage.
A call agrees with unit T when T is in its lineage. Names are never compared.

### Verdicts

**Kaiju** (Kraken 2 units), over the reads Kraken 2 placed in T's clade, and
**Kraken 2** (Kaiju-only units), over the reads Kaiju placed in T:

- *confirmed*: at least 5 reads agree, and they are at least half of the reads
  Kaiju resolved to T's rank
- *disagrees*: at least 5 reads resolved, less than half agree
- *insufficient*: otherwise. For a Kaiju-only unit this is the usual outcome:
  Kraken 2 left the reads unclassified, which is why only Kaiju detected the taxon.

Protein search rarely resolves species. A species unit whose Kaiju reads do not
reach species rank can still pass on its genus with the same thresholds
(`confirmed_at_genus`). Kraken 2 and BLAST must still agree at species rank.

**BLAST**: the best hits of a read are all hits at its top bitscore. A read
*agrees* at genus rank when a best hit lies inside T at ≥ 90 % identity. At
species rank it needs ≥ 97 % identity over ≥ 80 % of the read. Reads whose best
hits do not resolve to the rank ("uncultured bacterium") are *uninformative*.
The unit is:

- *confirmed*: at least 5 informative reads, and at least half of them agree
- *no_hit*: not confirmed, and at least half the reads hit nothing in nt
- *no_close_hit*: not confirmed, and most reads hit nothing, nothing resolvable,
  or nothing close. nt holds no close reference. **This is not evidence against
  the taxon**: nt has few whole-genome assemblies of environmental bacteria.
- *disagrees*: otherwise; most reads are close to another taxon

**Tier**: the detecting method, plus each other leg that confirmed → `triple`,
`double` or `single`.

NCBI moves species between genera, and databases built at different times
disagree about it. A genus unit therefore also accepts the current genus of any
species in its Kraken 2 clade that NCBI has since moved, when those species hold
at least 5 % of the unit's reads. In the cockle cohort that recovered *Moraxella*
(its *M. osloensis* reads are now *Faucicola*, which Kaiju and BLAST report).

### In the cockle cohort

Kraken 2 units: 563 samples, 7,656 units, 129,653 BLAST queries in 52 chunks:

| rank | triple | double | single |
|---|---|---|---|
| genus (4,619) | 3,687 (80 %) | 411 | 521 |
| species (3,037) | 1,637 (54 %) | 897 | 503 |

### Cost

Kaiju already runs in `PROFILE`; the confirmation reads its per-read output.
The Kraken 2 rerun takes seconds per sample. BLAST dominates: each chunk scans
nt once (~700 GB), then costs ~0.1 s per read. Measured with 12 threads, the scan
took ~24 min from a cold shared file system and ~2 min with nt in the page cache,
so a 2,500-read chunk takes ~28 min cold and ~6–7 min warm.

On Slurm, page cache is charged to the job that first reads a page. A 48 GB
BLAST task cannot hold nt, so it evicts nt while scanning it and every chunk
starts cold. A job that reads nt once and stays alive for the duration (e.g.
`--mem=760G`, `cat nt.* > /dev/null`, then `sleep`) keeps nt cached, and the
BLAST tasks read it without being charged for it. For the cockle cohort that is
~24 h of BLAST cold against ~6 h warm, before running chunks in parallel.
