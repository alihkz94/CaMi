#!/usr/bin/env python3
"""
=============================================================================
build_pairs.py — derive the within-individual tumour vs matched-host design
-----------------------------------------------------------------------------
PURPOSE
    Find every individual that contributed BOTH a tumour sample and a
    matched-host sample, and grade each pair by how confounded it is.

    A within-individual paired design controls for individual genotype, site,
    and batch all at once — none of which an unpaired tumour-vs-healthy
    contrast controls for.

    It does NOT fix a tissue confound. Where the tumour and host samples are
    different tissues, any difference found is a tissue difference as much as a
    disease one, and the report says so.

HOW PAIRS ARE FOUND
    By the `individual` column of sample_metadata.tsv. Two samples that share a
    value are two samples from one animal:

        sample  individual  tissue      disease_status
        T01     animal_01   haemolymph  tumor
        H01     animal_01   foot        matched host

    That column comes either from the samplesheet given to --input, or, on the
    ENA route, from the sample-code grammar (EPCE18_851H -> EPCE18_851). Either
    way this script only reads the column.

    A cohort that states no `individual` has no paired design. That is reported
    plainly, with the columns that would create one — it is a normal outcome,
    not a failure.

PAIR TIERS — the whole point of this script
    A pair is only as good as its worst confound. Prep (native vs WGA) is the
    one confound a paired design can still be ruined by, so it drives the tier:

      tier 1  prep-matched (both native) AND host tissue = foot
              -> cleanest contrast this cohort can produce
      tier 2  prep-matched (both native), host tissue = anything else
              -> clean prep, but a second tissue type enters the comparison
      tier 3  prep-matched (both WGA)
              -> no prep confound WITHIN the pair, but WGA distorts abundance
      tier 4  prep-MISMATCHED (native vs WGA)
              -> disease is confounded with library prep inside the pair.
                 Reported for completeness; excluded from the default analysis.

    Tier 1 alone is the headline design. Tiers 1+2 is the sensitivity analysis.
    Tier 4 must never be pooled with the others without saying so explicitly.

INPUTS
    <tables>/sample_metadata.tsv   written by aggregate_counts.py

OUTPUTS (under --outdir, default <tables>/../10_pairs)
    pairs.tsv              one row per pair, with tier and both sample codes
    pairs_report.txt       tier counts, unpairable samples, exclusions

ENV     polars.  RUN via Nextflow: --step STATS
VERSION 1.0 (2026-08-20)
=============================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import polars as pl

TUMOUR_STATUS: Final[frozenset[str]] = frozenset({"tumor", "tumour"})
HOST_STATUS: Final[frozenset[str]] = frozenset({"matched host", "unmatched host"})
HEALTHY_STATUS: Final[frozenset[str]] = frozenset({"healthy", "control", "normal"})
# The whole documented vocabulary. Anything outside it is a typo or a cohort
# this design does not describe, and either way the person needs to be told —
# an unrecognised value used to produce zero pairs and no explanation.
KNOWN_STATUS: Final[frozenset[str]] = TUMOUR_STATUS | HOST_STATUS | HEALTHY_STATUS

# Tier 1 is "host tissue is foot". The cohort this was written for spells that
# as a one-letter code; anyone else writes the word. Both must reach tier 1, or
# the headline design silently degrades to tier 2 on somebody else's data.
FOOT_LABELS: Final[frozenset[str]] = frozenset({"f", "p", "foot"})

# The columns WITHOUT WHICH THERE IS NO PAIRED DESIGN AT ALL. Missing ones are
# named in the report rather than crashing three functions later on a null.
#
# dna_prep and tissue are deliberately NOT here. They grade a pair, they do not
# create one: without them the pairs still exist and are still testable, and the
# tier simply records that the prep was not stated. Requiring them turned an
# ordinary cohort — one that states who is paired with whom and nothing else —
# into "no design", which is a different and false statement.
REQUIRED_FOR_PAIRING: Final[tuple[str, ...]] = ("individual", "disease_status")


def is_foot(tissue: str) -> bool:
    return (tissue or "").strip().lower() in FOOT_LABELS


def assign_tier(t_prep: str, h_prep: str, h_tissue: str) -> tuple[int, str]:
    """
    Grade one pair. Prep mismatch dominates every other consideration.

    UNSTATED PREP IS NOT WGA. Treating a blank cell as "not native" reported a
    prep confound that nobody had claimed, and pushed the pair to tier 4 where
    the default test excludes it. An unknown prep is graded on the tissue and
    says it is unknown, so the person can see what to fill in.
    """
    t = (t_prep or "").strip().lower()
    h = (h_prep or "").strip().lower()

    if not t or not h:
        if is_foot(h_tissue):
            return 1, "host tissue foot; prep not stated for one or both samples"
        return 2, f"host tissue {h_tissue or 'not stated'}; prep not stated"

    if (t == "native") != (h == "native"):
        return 4, "prep mismatch inside the pair (native vs WGA)"
    if t == "native":
        if is_foot(h_tissue):
            return 1, "both native, host tissue foot"
        return 2, f"both native, host tissue {h_tissue}"
    return 3, "both WGA (no prep confound inside the pair, but WGA distorts abundance)"


def usable_design(meta) -> list[str]:
    """
    Which pairing columns this cohort actually states.

    A column that is present but entirely empty counts as absent: a samplesheet
    with an `individual` header and no values describes no pairing, and saying
    "the column is there" would be true and useless.
    """
    stated = []
    for col in REQUIRED_FOR_PAIRING:
        if col not in meta.columns:
            continue
        values = [v for v in meta[col].to_list() if v is not None and str(v).strip()]
        if values:
            stated.append(col)
    return stated


def no_design_report(outdir, missing: list[str], meta_height: int,
                     note: list[str] | None = None) -> None:
    """
    Say plainly that there is no paired design, and what would create one.

    This is a normal outcome, not a failure. Profiling a cohort that has no
    tumour/host pairs is an ordinary thing to do, and the pipeline must finish
    and say so rather than crash on a null or write an empty table that reads
    like "we looked and found nothing".
    """
    (outdir / "pairs.tsv").write_text(
        "individual\ttumour_code\thost_code\ttier\ttier_reason\n"
    )
    text = "\n".join([
        "PAIRED-DESIGN REPORT",
        "=" * 76,
        "NOT RUN — this cohort states no within-individual paired design.",
        "",
        f"samples          : {meta_height}",
    ] + ([f"columns missing  : {', '.join(missing)}"] if missing else []) + [
        "",
    ] + (note + [""] if note else []) + [
        "The paired test compares a tumour sample against a matched-host sample",
        "from the SAME individual. To build it, give these columns in the",
        "samplesheet passed to --input:",
        "",
        "    sample,fastq_1,fastq_2,individual,tissue,disease_status,dna_prep",
        "",
        "  individual      the same value for both samples from one animal",
        f"  disease_status  one of: {', '.join(sorted(KNOWN_STATUS))}",
        "  dna_prep        native or wga",
        "  tissue          a name, or a one-letter code",
        "",
        "Everything else in --step STATS ran normally. Only the paired test",
        "needs a design.",
        "",
    ])
    (outdir / "pairs_report.txt").write_text(text)
    print(text)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Derive within-individual tumour vs matched-host pairs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--tables", type=Path, required=True,
                    help="directory holding sample_metadata.tsv")
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--require-analysed", action="store_true", default=True,
                    help="only pair samples that passed the integrity gate")
    ap.add_argument("--allow-unanalysed", dest="require_analysed",
                    action="store_false",
                    help="pair every manifest sample, profiled or not (planning only)")
    args = ap.parse_args()

    tables = args.tables.resolve()
    outdir = (args.outdir or tables.parent / "10_pairs").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    meta_path = tables / "sample_metadata.tsv"
    if not meta_path.exists():
        sys.exit(f"ERROR: {meta_path} not found. Run aggregate_counts.py first.")

    meta = pl.read_csv(meta_path, separator="\t")

    print("=" * 76)
    print(f" metadata : {meta_path}")
    print(f" outdir   : {outdir}")
    print("=" * 76)

    total = meta.height
    if args.require_analysed and "analysed" in meta.columns:
        meta = meta.filter(pl.col("analysed"))
        print(f"  restricted to profiled samples: {meta.height} of {total}")
    else:
        print(f"  using every manifest sample: {total} (profiled or not)")

    # Is there a design at all? Ask before touching any of these columns: on a
    # cohort that states none they are absent or entirely null, and every
    # expression below would fail on a null or join nulls to each other.
    stated = usable_design(meta)
    missing = [c for c in REQUIRED_FOR_PAIRING if c not in stated]
    if missing:
        no_design_report(outdir, missing, meta.height)
        return 0

    status = pl.col("disease_status").str.strip_chars().str.to_lowercase()

    # An unrecognised disease_status used to produce zero pairs and no reason.
    # "Tumour"/"TUMOR" are handled by the lowercasing above; a genuinely unknown
    # word is a typo or a different vocabulary, and is worth stopping for.
    seen = {
        (v or "").strip().lower()
        for v in meta["disease_status"].to_list() if v is not None
    } - {""}
    unknown = sorted(seen - KNOWN_STATUS)

    tumours = meta.filter(status.is_in(list(TUMOUR_STATUS)))
    hosts = meta.filter(status.is_in(list(HOST_STATUS)))

    # An unrecognised value pairs nothing, which on its own looks exactly like a
    # cohort that has no tumours. Say which values were not understood — but do
    # NOT exit non-zero: that fails BUILD_PAIRS and makes Nextflow discard
    # everything AGGREGATE_COUNTS produced, over a spelling.
    if unknown and not (tumours.height and hosts.height):
        no_design_report(
            outdir, [], meta.height,
            note=[
                f"disease_status holds {len(unknown)} unrecognised value(s): "
                f"{', '.join(unknown)}",
                f"Recognised: {', '.join(sorted(KNOWN_STATUS))}",
                "Only 'tumor'/'tumour' paired with 'matched host'/'unmatched host'",
                "forms a pair. Nothing here matched, so no pair could be built.",
            ],
        )
        return 0
    if unknown:
        print(f"  WARN: ignoring unrecognised disease_status value(s): "
              f"{', '.join(unknown)}", file=sys.stderr)
    print(f"  tumour samples : {tumours.height}")
    print(f"  host samples   : {hosts.height}")

    # An individual may in principle contribute more than one sample per side.
    # Join on individual and keep every combination, then report duplicates
    # rather than silently choosing one.
    pairs = tumours.join(hosts, on="individual", how="inner", suffix="_host")
    if pairs.height == 0:
        (outdir / "pairs.tsv").write_text(
            "individual\ttumour_code\thost_code\ttier\ttier_reason\n"
        )
        (outdir / "pairs_report.txt").write_text(
            "No tumour/matched-host pair found in the analysed set yet.\n"
        )
        print("  no pairs found yet (the run is probably still in progress)")
        return 0

    graded = [
        assign_tier(r["dna_prep"], r["dna_prep_host"], r["tissue_host"])
        for r in pairs.iter_rows(named=True)
    ]
    pairs = pairs.with_columns(
        pl.Series("tier", [g[0] for g in graded], dtype=pl.Int64),
        pl.Series("tier_reason", [g[1] for g in graded]),
    )

    out = (
        pairs.select(
            "individual",
            pl.col("sample_code").alias("tumour_code"),
            pl.col("sample_code_host").alias("host_code"),
            pl.col("tissue").alias("tumour_tissue"),
            pl.col("tissue_host").alias("host_tissue"),
            pl.col("dna_prep").alias("tumour_prep"),
            pl.col("dna_prep_host").alias("host_prep"),
            pl.col("disease_status_host").alias("host_status"),
            pl.col("site"),
            pl.col("nonhost_pairs").alias("tumour_nonhost_pairs"),
            pl.col("nonhost_pairs_host").alias("host_nonhost_pairs"),
            "tier",
            "tier_reason",
        )
        .sort("tier", "individual")
    )
    out.write_csv(outdir / "pairs.tsv", separator="\t")

    # --- report -----------------------------------------------------------
    tier_counts = out.group_by("tier").agg(pl.len().alias("n")).sort("tier")
    dup_t = (
        out.group_by("individual").agg(pl.col("tumour_code").n_unique().alias("n"))
        .filter(pl.col("n") > 1)
    )
    dup_h = (
        out.group_by("individual").agg(pl.col("host_code").n_unique().alias("n"))
        .filter(pl.col("n") > 1)
    )
    unpairable = (
        tumours.filter(~pl.col("individual").is_in(out["individual"].to_list()))
        .select("sample_code", "tissue", "dna_prep")
        .sort("sample_code")
    )

    lines = [
        "PAIRED-DESIGN REPORT",
        "=" * 76,
        f"tumour samples considered : {tumours.height}",
        f"host samples considered   : {hosts.height}",
        f"pairs formed              : {out.height}",
        "",
        "Pairs by tier:",
    ]
    labels = {
        1: "tier 1  both native, host = foot        <- headline design",
        2: "tier 2  both native, host = other tissue",
        3: "tier 3  both WGA",
        4: "tier 4  prep mismatch  <- EXCLUDED from the default analysis",
    }
    for r in tier_counts.iter_rows(named=True):
        lines.append(f"  {labels.get(r['tier'], 'tier ' + str(r['tier'])):<52s} n={r['n']}")

    usable = out.filter(pl.col("tier") <= 2).height
    lines += [
        "",
        f"Usable for the default paired test (tier 1-2): n = {usable}",
    ]
    if usable < 6:
        lines.append("  CAUTION: fewer than 6 pairs. An exact paired test has very")
        lines.append("  limited power here; treat any result as descriptive.")

    if dup_t.height or dup_h.height:
        lines += ["", "AMBIGUOUS individuals (more than one sample on a side):"]
        for ind in sorted(set(dup_t["individual"].to_list()) | set(dup_h["individual"].to_list())):
            rows = out.filter(pl.col("individual") == ind)
            lines.append(f"  {ind}: " + "; ".join(
                f"{r['tumour_code']} vs {r['host_code']}" for r in rows.iter_rows(named=True)
            ))
        lines.append("  All combinations are listed in pairs.tsv. Decide explicitly")
        lines.append("  which to keep before testing — do not average them silently.")

    if unpairable.height:
        lines += ["", f"Tumour samples with NO host partner ({unpairable.height}):"]
        codes = unpairable["sample_code"].to_list()
        lines += ["  " + ", ".join(codes[i:i + 6]) for i in range(0, len(codes), 6)]

    report = "\n".join(lines)
    (outdir / "pairs_report.txt").write_text(report + "\n")
    print()
    print(report)
    print()
    print(f"WROTE -> {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
