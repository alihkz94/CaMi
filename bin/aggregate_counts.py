#!/usr/bin/env python3
"""
=============================================================================
aggregate_counts.py — one count table from every profiled sample
-----------------------------------------------------------------------------
PURPOSE
    Turn per-sample Bracken and Kaiju outputs into tidy, analysis-ready tables.
    Step 1 of --step STATS.

    Reads ONLY tool outputs (Bracken TSV, Kaiju TSV, Kraken 2 report, fastp
    JSON). It never re-scans a FASTQ to count reads — a project convention, and
    the difference between seconds and hours across 563 samples.

WHY THIS SCRIPT IS DEFENSIVE
    Aggregation is where silent corruption happens. Four traps were measured on
    this dataset before this script was written, and each one is guarded here:

    1. IN-FLIGHT FILES. publishDir copies asynchronously while the pipeline
       runs, so a file can be read half-written. Guarded by RECONCILIATION:
       Kraken 2 and Kaiju both count PAIRS, and on a complete sample

           kraken2(classified + unclassified) == kaiju(summary total)
                                              == fastp(total_reads) / 2

       Verified exactly on real samples (e.g. ASCE17_1940F: 3,142,213 all
       three ways). A sample that fails this is truncated, not unlucky, and is
       quarantined rather than averaged into the cohort.

    2. TAXON-NAME COLLISIONS. Pivoting on taxon_name would silently merge two
       distinct taxon_ids that share a name. This script pivots on TAXON_ID and
       writes a separate name dictionary. (Zero collisions exist today; that is
       a property of the current data, not a guarantee, so it is not relied on.)

    3. BRACKEN DOES NOT SUM TO KRAKEN. Bracken re-estimates only the reads it
       can push to the requested rank: on ASCE17_1909F, 49,435 of 54,919
       classified pairs reach species. Bracken's own fraction_total_reads is
       therefore relative to a per-sample, per-rank denominator and is NOT
       comparable between samples. This script recomputes every fraction
       against the sample's own non-host pair count instead.

    4. AMBIGUOUS SAMPLE IDENTITY. Manifests list one row per mate, and a sample
       could in principle appear in both. Deduplication is explicit and any
       overlap is reported, never silently resolved.

INPUTS (under --results, default <project>/results)
    04_dedup/<sample>.dedup.json              non-host depth (the denominator)
    05_kraken2/<sample>.kraken2.report        classified/unclassified totals
    06_bracken/<sample>.bracken.{species,genus}.tsv
    07_kaiju/<sample>.kaiju.{species,phylum}.tsv
    07_kaiju/<sample>.kaiju.summary.tsv       domain totals
  and ONE source of per-sample design, whichever the run used:
    --design <file>                           the samplesheet, via --input
    <project>/manifests/sample_name_map.csv   the ENA route
    <project>/manifests/skipped_bam_only.tsv  the ENA route, BAM-only runs

OUTPUTS (under --outdir, default <results>/09_tables)
    sample_metadata.tsv/.parquet      one row per sample + derived fields
    counts_long.parquet               tidy counts, all methods and ranks
    counts_<method>_<rank>_wide.tsv   sample x taxon_id matrix
    taxa_<method>_<rank>.tsv          taxon_id -> name dictionary
    integrity_report.tsv              per-sample reconciliation, PASS/FAIL
    aggregate_report.txt              human-readable summary

ENV     polars. RUN via Nextflow: --step STATS
VERSION 1.0 (2026-08-20)
=============================================================================
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import polars as pl

# Homo is a database artifact here, not contamination: the Kraken standard-16
# database holds a single eukaryote, so cockle reads with no better match land
# on Homo. Flagged, never silently dropped — the statistics step decides.
HOST_ARTIFACT_TAXA: Final[frozenset[int]] = frozenset({9606, 9605})

# Reconciliation tolerance. Exact equality is what a complete sample shows;
# a non-zero tolerance exists only to absorb a classifier that legitimately
# reports one rounding unit, never to excuse a truncated file.
RECONCILE_TOL: Final[int] = 0

# The Bruzos biobank codebook, confirmed against the ENA submission workbook.
# Used ONLY when --sample-code-scheme is bruzos; another cohort's sample names
# are never read as tissue codes.
#
# "A" is adductor muscle: the earlier "adductor/unlabelled" hedge is resolved —
# all 10 matched-host A samples and 1 unmatched-host A sample are adductor.
# "C" is NOT a tissue. It only ever appeared because "_HC" was not stripped;
# see _HIGH_COVERAGE. Letters the codebook defines but this cohort never uses
# (D digestive, G gonad, B gill, P foot) are listed so a new sample cannot
# arrive and be silently mislabelled.
TISSUE_NAMES: Final[dict[str, str]] = {
    "H": "haemolymph",
    "F": "foot",
    "P": "foot",
    "M": "mantle",
    "A": "adductor muscle",
    "D": "digestive/intestine",
    "G": "gonad",
    "B": "gill",
}

COUNT_SETS: Final[tuple[tuple[str, str], ...]] = (
    ("bracken", "species"),
    ("bracken", "genus"),
    ("kaiju", "species"),
    ("kaiju", "phylum"),
)


# --------------------------------------------------------------------------
# sample-code grammar
# --------------------------------------------------------------------------
# Codes look like  EPCE18_851H  or  FRCE17_840F-wga :
#     [SITE:2][CE][YEAR:2] _ [INDIVIDUAL] [TISSUE:1 letter] [-wga]
# The trailing letter is the TISSUE, and it is the ONLY thing that differs
# between a tumour sample and its matched-host partner from the same cockle.
# Strip -wga first or the tissue letter stays hidden behind it.
@dataclass(frozen=True, slots=True)
class CodeParts:
    individual: str
    tissue: str
    tissue_name: str
    site: str
    replicate: str


# Prep suffixes seen in this cohort: "-wga" (53 samples) and "-gp" (1 sample).
# Matched generically as a lowercase tail so a third one cannot slip past.
_PREP_SUFFIX = re.compile(r"-[a-z]+$")
# Tissue letter, optionally followed by a replicate number: "...851H", "...421H1".
# The replicate group is what makes PACE17_421H1 pair with PACE17_421M-wga.
# Codes that end in digits with no capital before them (ASCE17_1983) correctly
# fail this and are recorded as tissue "none" rather than being mis-split.
_CODE_GRAMMAR = re.compile(r"^(?P<ind>.*?)(?P<tissue>[A-Z])(?P<rep>\d*)$")
# "_HC" is a HIGH-COVERAGE annotation, not a tissue. Five healthy samples carry
# it (ICCE19_359F_HC and four more, all 60 Gb). Read naively the code ends in C,
# and C is not a tissue in the Bruzos codebook at all — it is a country letter.
# Left unstripped it turns five healthy FOOT samples into an invented
# "unlabelled-C" tissue and drops them out of every foot comparison. Strip it
# BEFORE the prep suffix and the tissue letter, so all three compose: the
# replicate group above must still fire for PACE17_421H1.
_HIGH_COVERAGE = re.compile(r"_HC$")


def split_code(code: str) -> CodeParts:
    """Break a sample code into individual / tissue / site. Never raises."""
    base = _PREP_SUFFIX.sub("", _HIGH_COVERAGE.sub("", code))
    if m := _CODE_GRAMMAR.match(base):
        individual, tissue, replicate = m["ind"], m["tissue"], m["rep"]
    else:
        # Some codes carry no tissue letter at all (healthy samples ending in a
        # number, and the one differently-named BNg14). Record that honestly
        # instead of forcing them into the grammar; they cannot be paired.
        individual, tissue, replicate = base, "none", ""
    return CodeParts(
        individual=individual,
        tissue=tissue,
        tissue_name=TISSUE_NAMES.get(tissue, tissue),
        site=code[:2] if len(code) >= 2 else "??",
        replicate=replicate,
    )


# --------------------------------------------------------------------------
# metadata
# --------------------------------------------------------------------------
# TWO SOURCES, ONE SHAPE.
#
#   a samplesheet   --design, written by Nextflow from --input. The general
#                   route: the person running CaMi states the design, and every
#                   column is optional.
#   ENA manifests   the route this pipeline was built on, where the design is
#                   encoded in the sample code and read by the Bruzos grammar.
#
# Both produce the same columns, so nothing downstream knows which was used.
# The samplesheet wins when it carries rows, because a person who wrote one
# meant it.
DESIGN_FIELDS: Final[tuple[str, ...]] = (
    "individual", "tissue", "disease_status", "dna_prep",
)


def tissue_label(value: str) -> str:
    """
    A tissue as a person would read it.

    A single capital is a Bruzos code and is expanded. Anything else is already
    a name — "foot", "gill", "whole animal" — and is passed through untouched,
    because inventing a code for it would only create a second vocabulary.
    """
    if len(value) == 1 and value.isupper():
        return TISSUE_NAMES.get(value, value)
    return value


def stated_values(frame: pl.DataFrame, column: str) -> bool:
    """
    Does this column actually say anything?

    A column that is absent, all null, or all whitespace states nothing. All
    three must read the same way, or a cohort gets told one field is missing
    and an equally empty one beside it is present.
    """
    if column not in frame.columns:
        return False
    return any(v is not None and str(v).strip() for v in frame[column].to_list())


def load_design(path: Path, scheme: str) -> pl.DataFrame | None:
    """
    The samplesheet design, or None when there is none.

    None is not "an empty design". It means no samplesheet was involved, and the
    caller must fall back to the manifests. A header with no rows says the same
    thing — that is exactly what Nextflow writes on the ENA route.
    """
    if not path or not path.exists():
        return None
    rows = [r for r in path.read_text().splitlines() if r.strip()]
    if len(rows) < 2:
        return None

    # infer_schema_length=0 reads EVERY column as text, and that is the point.
    # These are labels, not measurements. Left to infer, polars turns the sample
    # names "001,002" into the integers 1 and 2 — and those names are what every
    # output file is called, so each sample then looks ABSENT and the run dies
    # blaming the profiling. Numeric `individual` values are worse: they compare
    # as i64 against a string and raise from somewhere unrelated.
    df = pl.read_csv(
        path, separator="\t", infer_schema_length=0,
        null_values=["", "NA", "na", "None"],
    )
    if "sample" not in df.columns:
        sys.exit(f"ERROR: {path} has no 'sample' column. Nextflow writes this "
                 f"file; it should not be edited by hand.")

    for field in DESIGN_FIELDS:
        if field not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias(field))

    df = df.rename({"sample": "sample_code"}).unique(
        subset=["sample_code"], keep="first", maintain_order=True
    )

    # The grammar fills only what the sheet left blank, and only when asked for.
    # Silently parsing another cohort's names as tissue codes is exactly the
    # failure this option exists to prevent.
    codes = df["sample_code"].to_list()
    parts = [split_code(c) for c in codes] if scheme == "bruzos" else None

    def fallback(field: str, values: list[str]) -> pl.Series:
        col = df[field].to_list()
        return pl.Series(field, [c if c else v for c, v in zip(col, values)])

    if parts is not None:
        df = df.with_columns(
            fallback("individual", [p.individual for p in parts]),
            fallback("tissue", [p.tissue for p in parts]),
        )
        site = [p.site for p in parts]
        replicate = [p.replicate for p in parts]
    else:
        site = ["" for _ in codes]
        replicate = ["" for _ in codes]

    # An unstated tissue stays NULL, not the string "none". The two are not the
    # same: "none" is a value, and a column full of it reads to every later check
    # as "the tissue was stated", which is how an unstated column came to be
    # reported as present while dna_prep beside it was reported as missing.
    tissues = [t or None for t in df["tissue"].to_list()]
    return df.with_columns(
        pl.lit("-").alias("run_accession"),
        pl.lit("SAMPLESHEET").alias("source"),
        pl.Series("tissue", tissues, dtype=pl.Utf8),
        pl.Series("tissue_name", [tissue_label(t) if t else None for t in tissues],
                  dtype=pl.Utf8),
        pl.Series("site", [s or None for s in site], dtype=pl.Utf8),
        pl.Series("replicate", [r or None for r in replicate], dtype=pl.Utf8),
    ).select(
        "run_accession", "sample_code", "disease_status", "dna_prep", "source",
        "individual", "tissue", "tissue_name", "site", "replicate",
    )


def load_manifest(project: Path, scheme: str = "bruzos") -> tuple[pl.DataFrame, list[str]]:
    """Both manifests as one row per sample. Returns (table, overlapping codes)."""
    frames: list[pl.DataFrame] = []
    sets: dict[str, set[str]] = {}

    fastq_map = project / "manifests" / "sample_name_map.csv"
    if fastq_map.exists():
        df = (
            pl.read_csv(fastq_map)
            .select(
                "run_accession", "sample_code", "disease_status", "dna_prep",
                pl.lit("FASTQ").alias("source"),
            )
            .unique(subset=["sample_code"], keep="first", maintain_order=True)
        )
        sets["FASTQ"] = set(df["sample_code"].to_list())
        frames.append(df)

    bam_map = project / "manifests" / "skipped_bam_only.tsv"
    if bam_map.exists():
        df = (
            pl.read_csv(bam_map, separator="\t")
            .select(
                "run_accession", "sample_code", "disease_status", "dna_prep",
                pl.lit("BAM").alias("source"),
            )
            .unique(subset=["sample_code"], keep="first", maintain_order=True)
        )
        sets["BAM"] = set(df["sample_code"].to_list())
        frames.append(df)

    if not frames:
        sys.exit(
            f"ERROR: no manifest found under {project / 'manifests'}, and no "
            f"samplesheet was given.\n"
            f"  Run with --input samples.csv, or build the manifests with "
            f"scripts/download/."
        )

    overlap = sorted(sets.get("FASTQ", set()) & sets.get("BAM", set()))
    meta = pl.concat(frames, how="vertical").unique(
        subset=["sample_code"], keep="first", maintain_order=True
    )

    codes = meta["sample_code"].to_list()
    parts = ([split_code(c) for c in codes] if scheme == "bruzos"
             else [CodeParts(c, "none", "none", "", "") for c in codes])
    return (
        meta.with_columns(
            pl.Series("individual", [p.individual for p in parts]),
            pl.Series("tissue", [p.tissue for p in parts]),
            pl.Series("tissue_name", [p.tissue_name for p in parts]),
            pl.Series("site", [p.site for p in parts]),
            pl.Series("replicate", [p.replicate for p in parts]),
        ),
        overlap,
    )


# --------------------------------------------------------------------------
# per-sample integrity: the three independent pair counts must agree
# --------------------------------------------------------------------------
def fastp_pairs(path: Path) -> int | None:
    try:
        total = json.loads(path.read_text())["summary"]["after_filtering"]["total_reads"]
        return int(total) // 2
    except Exception:
        return None


def kraken_pairs(path: Path) -> int | None:
    """Classified (root, rank R) + unclassified (rank U) = every input pair."""
    try:
        classified = unclassified = None
        with path.open() as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                if f[3] == "U" and unclassified is None:
                    unclassified = int(f[1])
                elif f[3] == "R" and classified is None:
                    classified = int(f[1])
                if classified is not None and unclassified is not None:
                    break
        if classified is None and unclassified is None:
            return None
        return (classified or 0) + (unclassified or 0)
    except Exception:
        return None


def kaiju_pairs(path: Path) -> int | None:
    """kaiju.summary.tsv is headerless <label>\\t<reads>; rows sum to all pairs."""
    try:
        total = 0
        seen = False
        for line in path.read_text().splitlines():
            f = line.split("\t")
            if len(f) == 2 and f[1].strip().lstrip("-").isdigit():
                total += int(f[1].strip())
                seen = True
        return total if seen else None
    except Exception:
        return None


def build_integrity(results: Path, samples: list[str]) -> pl.DataFrame:
    rows = []
    for s in samples:
        fp = fastp_pairs(results / "04_dedup" / f"{s}.dedup.json")
        kr = kraken_pairs(results / "05_kraken2" / f"{s}.kraken2.report")
        kj = kaiju_pairs(results / "07_kaiju" / f"{s}.kaiju.summary.tsv")

        present = [v for v in (fp, kr, kj) if v is not None]
        if not present:
            status, note = "ABSENT", "no classifier output yet"
        elif fp is None:
            status, note = "NO_DEPTH", "dedup JSON missing; cannot normalise"
        else:
            bad = []
            if kr is not None and abs(kr - fp) > RECONCILE_TOL:
                bad.append(f"kraken2={kr}")
            if kj is not None and abs(kj - fp) > RECONCILE_TOL:
                bad.append(f"kaiju={kj}")
            if bad:
                status = "FAIL"
                note = f"pair counts disagree with fastp={fp}: " + ", ".join(bad)
            elif kr is None and kj is None:
                status, note = "ABSENT", "depth only, no classifier output"
            else:
                status = "PASS"
                note = "reconciled" + ("" if kr is not None else " (kraken2 pending)") \
                                    + ("" if kj is not None else " (kaiju pending)")
        rows.append(
            {
                "sample_code": s,
                "fastp_pairs": fp,
                "kraken2_pairs": kr,
                "kaiju_pairs": kj,
                "integrity": status,
                "note": note,
            }
        )
    return pl.DataFrame(
        rows,
        schema={
            "sample_code": pl.Utf8, "fastp_pairs": pl.Int64,
            "kraken2_pairs": pl.Int64, "kaiju_pairs": pl.Int64,
            "integrity": pl.Utf8, "note": pl.Utf8,
        },
    )


# --------------------------------------------------------------------------
# counts
# --------------------------------------------------------------------------
COUNTS_SCHEMA: Final[dict[str, pl.DataType]] = {
    "sample_code": pl.Utf8, "method": pl.Utf8, "rank": pl.Utf8,
    "taxon_id": pl.Int64, "taxon_name": pl.Utf8, "reads": pl.Int64,
}


def _empty_counts() -> pl.DataFrame:
    return pl.DataFrame(schema=COUNTS_SCHEMA)


def read_bracken(results: Path, rank: str, keep: set[str]) -> pl.DataFrame:
    """Bracken TSV -> tidy counts. new_est_reads is Bracken's re-estimate."""
    frames = []
    for tf in sorted((results / "06_bracken").glob(f"*.bracken.{rank}.tsv")):
        sample = tf.name.removesuffix(f".bracken.{rank}.tsv")
        if sample not in keep:
            continue
        try:
            df = pl.read_csv(tf, separator="\t")
            if df.height == 0:
                continue
            frames.append(
                df.select(
                    pl.lit(sample).alias("sample_code"),
                    pl.lit("bracken").alias("method"),
                    pl.lit(rank).alias("rank"),
                    pl.col("taxonomy_id").cast(pl.Int64).alias("taxon_id"),
                    pl.col("name").cast(pl.Utf8).alias("taxon_name"),
                    pl.col("new_est_reads").cast(pl.Int64).alias("reads"),
                )
            )
        except Exception as e:
            print(f"  WARN: unreadable {tf.name}: {e}", file=sys.stderr)
    return pl.concat(frames, how="vertical") if frames else _empty_counts()


def read_kaiju(results: Path, rank: str, keep: set[str]) -> pl.DataFrame:
    """Kaiju TSV -> tidy counts, minus its unclassified / cannot-assign rows."""
    frames = []
    for tf in sorted((results / "07_kaiju").glob(f"*.kaiju.{rank}.tsv")):
        sample = tf.name.removesuffix(f".kaiju.{rank}.tsv")
        if sample not in keep:
            continue
        try:
            # taxon_id holds the literal string "NA" on summary rows, so it must
            # be read as text and filtered before any numeric cast.
            df = pl.read_csv(tf, separator="\t", schema_overrides={"taxon_id": pl.Utf8})
            if df.height == 0:
                continue
            df = df.filter(
                pl.col("taxon_id").is_not_null()
                & (pl.col("taxon_id") != "NA")
                & pl.col("taxon_name").is_not_null()
                & ~pl.col("taxon_name").str.contains(
                    r"(?i)unclassified|cannot be assigned"
                )
            )
            if df.height == 0:
                continue
            frames.append(
                df.select(
                    pl.lit(sample).alias("sample_code"),
                    pl.lit("kaiju").alias("method"),
                    pl.lit(rank).alias("rank"),
                    pl.col("taxon_id").cast(pl.Int64, strict=False).alias("taxon_id"),
                    pl.col("taxon_name").cast(pl.Utf8).alias("taxon_name"),
                    pl.col("reads").cast(pl.Int64).alias("reads"),
                ).filter(pl.col("taxon_id").is_not_null())
            )
        except Exception as e:
            print(f"  WARN: unreadable {tf.name}: {e}", file=sys.stderr)
    return pl.concat(frames, how="vertical") if frames else _empty_counts()


def read_kaiju_domains(results: Path, keep: set[str]) -> pl.DataFrame:
    rows = []
    for tf in sorted((results / "07_kaiju").glob("*.kaiju.summary.tsv")):
        sample = tf.name.removesuffix(".kaiju.summary.tsv")
        if sample not in keep:
            continue
        try:
            for line in tf.read_text().splitlines():
                f = line.split("\t")
                if len(f) == 2 and f[1].strip().lstrip("-").isdigit():
                    rows.append(
                        {"sample_code": sample, "domain": f[0].strip(),
                         "reads": int(f[1].strip())}
                    )
        except Exception as e:
            print(f"  WARN: unreadable {tf.name}: {e}", file=sys.stderr)
    return pl.DataFrame(
        rows, schema={"sample_code": pl.Utf8, "domain": pl.Utf8, "reads": pl.Int64}
    )


# --------------------------------------------------------------------------
def main() -> int:
    here = Path(__file__).resolve()
    default_project = here.parents[1].parent  # <repo>/bin -> <repo> -> dataRoot

    ap = argparse.ArgumentParser(
        description="Aggregate per-sample classifier output into one count table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--project", type=Path, default=default_project)
    ap.add_argument("--results", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--min-depth", type=int, default=0,
                    help="drop samples with fewer non-host pairs than this")
    ap.add_argument("--include-failed", action="store_true",
                    help="include samples that FAIL reconciliation (not advised)")
    ap.add_argument("--design", type=Path, default=None,
                    help="samplesheet design written by Nextflow from --input; "
                         "a header with no rows means the ENA manifests are used")
    ap.add_argument("--sample-code-scheme", choices=("bruzos", "none"), default=None,
                    help="how to read a sample name. 'bruzos' applies the cockle "
                         "study's grammar (EPCE18_851H -> individual EPCE18_851, "
                         "tissue H); 'none' treats a name as a name. The default "
                         "is 'bruzos' on the manifest route and 'none' otherwise")
    args = ap.parse_args()

    project = args.project.resolve()
    results = (args.results or project / "results").resolve()
    outdir = (args.outdir or results / "09_tables").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not results.is_dir():
        sys.exit(f"ERROR: results directory not found: {results}")

    print("=" * 76)
    print(f" project : {project}")
    print(f" results : {results}")
    print(f" outdir  : {outdir}")
    print("=" * 76)

    # 1 -- metadata, from the samplesheet when there is one, else the manifests
    design_path = args.design
    scheme = args.sample_code_scheme
    design = load_design(design_path, scheme or "none") if design_path else None

    if design is not None:
        meta, overlap = design, []
        print(f"  samplesheet      : {meta.height} sample(s)")
        # meta[f] is a Series, and Series.filter wants a boolean mask, not an
        # expression — pl.col() here raised a TypeError on the FIRST field, so
        # no samplesheet run ever reached the aggregation at all.
        stated = [f for f in DESIGN_FIELDS if stated_values(meta, f)]
        missing = [f for f in DESIGN_FIELDS if f not in stated]
        print(f"  design columns   : {', '.join(stated) if stated else 'none'}")
        if missing:
            # Not a warning. A cohort with no tumour/host pairing is a normal
            # thing to profile, and BUILD_PAIRS says what it cannot test.
            print(f"  not stated       : {', '.join(missing)} — "
                  f"the tests that need them will report why they cannot run")
    else:
        # The ENA route. The grammar is the design here, so it is the default.
        if scheme is None:
            scheme = "bruzos"
        meta, overlap = load_manifest(project, scheme)
        if overlap:
            print(f"  WARN: {len(overlap)} sample(s) in BOTH manifests, FASTQ kept: "
                  f"{', '.join(overlap[:5])}", file=sys.stderr)
        print(f"  manifest samples : {meta.height}")

    # 2 -- integrity gate
    integrity = build_integrity(results, meta["sample_code"].to_list())
    integrity.write_csv(outdir / "integrity_report.tsv", separator="\t")
    tally = dict(
        integrity.group_by("integrity").agg(pl.len().alias("n"))
        .iter_rows()  # (status, n)
    )
    for status in ("PASS", "FAIL", "NO_DEPTH", "ABSENT"):
        if status in tally:
            print(f"  integrity {status:9s}: {tally[status]}")

    usable = set(
        integrity.filter(
            pl.col("integrity").is_in(["PASS"] + (["FAIL"] if args.include_failed else []))
        )["sample_code"].to_list()
    )
    if not usable:
        sys.exit("ERROR: no sample passed the integrity gate. Nothing to aggregate.")

    # 3 -- depth (only for usable samples; it is the normalisation denominator)
    depth = integrity.filter(pl.col("sample_code").is_in(list(usable))).select(
        "sample_code",
        pl.col("fastp_pairs").alias("nonhost_pairs"),
        (pl.col("fastp_pairs") * 2).alias("nonhost_reads"),
    )
    meta = meta.join(depth, on="sample_code", how="left").join(
        integrity.select("sample_code", "integrity"), on="sample_code", how="left"
    )

    if args.min_depth > 0:
        shallow = meta.filter(
            pl.col("nonhost_pairs").is_not_null()
            & (pl.col("nonhost_pairs") < args.min_depth)
        )["sample_code"].to_list()
        if shallow:
            print(f"  depth gate <{args.min_depth}: dropping {len(shallow)} sample(s)")
        usable -= set(shallow)
    meta = meta.with_columns(pl.col("sample_code").is_in(list(usable)).alias("analysed"))

    # 4 -- counts
    parts = [
        read_bracken(results, "species", usable),
        read_bracken(results, "genus", usable),
        read_kaiju(results, "species", usable),
        read_kaiju(results, "phylum", usable),
    ]
    counts = pl.concat(parts, how="vertical")
    if counts.height == 0:
        sys.exit("ERROR: no Bracken or Kaiju output found for any usable sample.")

    counts = (
        counts.with_columns(
            pl.col("taxon_id").is_in(list(HOST_ARTIFACT_TAXA)).alias("is_host_artifact")
        )
        .join(meta.select("sample_code", "nonhost_pairs"), on="sample_code", how="left")
        .with_columns(
            # Normalised against the sample's OWN non-host pairs — never against
            # Bracken's per-rank denominator, which is not comparable between
            # samples (see trap 3 in the header).
            pl.when(pl.col("nonhost_pairs") > 0)
            .then(pl.col("reads") / pl.col("nonhost_pairs"))
            .otherwise(None)
            .alias("frac_of_nonhost")
        )
    )

    # 5 -- per-sample microbial totals
    micro = (
        counts.filter(
            (pl.col("method") == "bracken") & (pl.col("rank") == "species")
            & ~pl.col("is_host_artifact")
        )
        .group_by("sample_code")
        .agg(pl.col("reads").sum().alias("microbial_pairs_bracken"))
    )
    meta = meta.join(micro, on="sample_code", how="left")

    domains = read_kaiju_domains(results, usable)
    if domains.height:
        (domains.pivot(on="domain", index="sample_code", values="reads")
         .fill_null(0).sort("sample_code")
         .write_csv(outdir / "kaiju_domain_totals.tsv", separator="\t"))
        meta = meta.join(
            domains.filter(pl.col("domain") == "Bacteria")
            .select("sample_code", pl.col("reads").alias("kaiju_bacteria_pairs")),
            on="sample_code", how="left",
        )

    # 6 -- write tidy + wide (wide pivots on taxon_id; see trap 2)
    meta = meta.sort("sample_code")
    meta.write_csv(outdir / "sample_metadata.tsv", separator="\t")
    meta.write_parquet(outdir / "sample_metadata.parquet")
    counts.write_parquet(outdir / "counts_long.parquet")

    for method, rank in COUNT_SETS:
        sub = counts.filter((pl.col("method") == method) & (pl.col("rank") == rank))
        if sub.height == 0:
            continue
        taxa = (
            sub.select("taxon_id", "taxon_name")
            .unique(subset=["taxon_id"], keep="first")
            .sort("taxon_id")
        )
        taxa.write_csv(outdir / f"taxa_{method}_{rank}.tsv", separator="\t")

        wide = (
            sub.group_by("sample_code", "taxon_id").agg(pl.col("reads").sum())
            .pivot(on="taxon_id", index="sample_code", values="reads")
            .fill_null(0)
            .sort("sample_code")
        )
        wide.write_csv(outdir / f"counts_{method}_{rank}_wide.tsv", separator="\t")
        print(f"  {method:8s} {rank:8s}: {sub['sample_code'].n_unique():4d} samples "
              f"x {taxa.height:6d} taxa")

    # 7 -- report
    profiled = set(counts["sample_code"].unique().to_list())
    expected = set(meta["sample_code"].to_list())
    orphans = sorted(
        {p.name.split(".")[0] for p in (results / "06_bracken").glob("*.bracken.species.tsv")}
        - expected
    )

    lines = [
        "AGGREGATION REPORT",
        "=" * 76,
        f"samples in manifests            : {len(expected)}",
        f"samples passing integrity gate  : {len(usable)}",
        f"samples in the count table      : {len(profiled)}",
        f"samples not yet profiled        : {len(expected - profiled)}",
    ]
    if overlap:
        lines.append(f"WARNING sample(s) in both manifests: {', '.join(overlap)}")
    if orphans:
        lines.append(f"WARNING output with no manifest row: {', '.join(orphans[:10])}")
    fails = integrity.filter(pl.col("integrity") == "FAIL")
    if fails.height:
        lines += ["", f"QUARANTINED — failed reconciliation ({fails.height}):"]
        lines += [f"  {r['sample_code']}: {r['note']}" for r in fails.iter_rows(named=True)]

    lines += ["", "Cohort composition (analysed samples only):"]
    grp = (
        meta.filter(pl.col("sample_code").is_in(list(profiled)))
        .group_by("disease_status", "dna_prep", "tissue_name")
        .agg(pl.len().alias("n"))
        .sort("n", descending=True)
    )
    # str(... or "-"): an unstated column is null here, and formatting None with
    # ":<16s" raises. That crashed AFTER every table had been written, so the
    # task failed and Nextflow discarded work that had in fact succeeded.
    lines += [
        f"  {str(r['disease_status'] or '-'):<16s} {str(r['dna_prep'] or '-'):<8s} "
        f"{str(r['tissue_name'] or '-'):<20s} n={r['n']}"
        for r in grp.iter_rows(named=True)
    ]

    report = "\n".join(lines)
    (outdir / "aggregate_report.txt").write_text(report + "\n")
    print()
    print(report)
    print()
    print(f"WROTE -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
