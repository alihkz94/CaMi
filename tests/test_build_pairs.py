#!/usr/bin/env python3
"""
=============================================================================
test_build_pairs.py — the paired design on a cohort that is not the cockles
-----------------------------------------------------------------------------
build_pairs.py was written against one cohort, where the design is spelled in
one-letter codes and every column is always present. Three things then break
quietly on anybody else's data:

  no design stated   `individual` is absent or empty. Joining on a null column
                     is not an error in polars, so this produced either a crash
                     or an empty table that reads like "we looked and found no
                     pairs" — a different statement, and a false one.
  tissue as a word   tier 1 is "host tissue is foot". Comparing against the
                     letter "F" alone silently demotes every foot pair in a
                     cohort that writes the word to tier 2.
  an unknown status  a disease_status outside the vocabulary matched nothing,
                     so a typo looked exactly like a cohort with no tumours.

The cockle cohort must keep grading exactly as before; that is checked against
the real tables in the release procedure, not here.

RUN:  python3 -m unittest discover -s tests
=============================================================================
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "build_pairs.py"
sys.path.insert(0, str(REPO / "bin"))

import build_pairs as bp  # noqa: E402

COLUMNS = [
    "sample_code", "disease_status", "dna_prep", "individual", "tissue",
    "tissue_name", "site", "nonhost_pairs", "analysed",
]


def metadata(rows: list[dict]) -> Path:
    """Write a sample_metadata.tsv of the shape aggregate_counts.py produces."""
    tables = Path(tempfile.mkdtemp()) / "tables"
    tables.mkdir()
    out = ["\t".join(COLUMNS)]
    for r in rows:
        full = {"site": "XX", "nonhost_pairs": "1000", "analysed": "true",
                "tissue_name": r.get("tissue", ""), **r}
        out.append("\t".join(str(full.get(c, "")) for c in COLUMNS))
    (tables / "sample_metadata.tsv").write_text("\n".join(out) + "\n")
    return tables


def run(tables: Path) -> tuple[int, str, Path]:
    outdir = tables.parent / "pairs"
    p = subprocess.run(
        [sys.executable, str(SCRIPT), "--tables", str(tables), "--outdir", str(outdir)],
        capture_output=True, text=True,
    )
    return p.returncode, p.stdout + p.stderr, outdir


class TestFootIsFootHoweverItIsSpelled(unittest.TestCase):
    def test_the_letter_and_the_word_both_reach_tier_one(self):
        for label in ("F", "f", "foot", "Foot", "P"):
            self.assertEqual(
                bp.assign_tier("native", "native", label)[0], 1,
                f"host tissue {label!r} did not reach tier 1",
            )

    def test_another_tissue_is_tier_two(self):
        self.assertEqual(bp.assign_tier("native", "native", "mantle")[0], 2)
        self.assertEqual(bp.assign_tier("native", "native", "M")[0], 2)

    def test_prep_mismatch_dominates_the_tissue(self):
        """A foot pair with mismatched prep is still tier 4, not tier 1."""
        self.assertEqual(bp.assign_tier("native", "wga", "foot")[0], 4)

    def test_both_wga_is_tier_three(self):
        self.assertEqual(bp.assign_tier("wga", "wga", "foot")[0], 3)

    def test_an_unstated_prep_is_not_reported_as_wga(self):
        """
        A blank dna_prep used to grade as "not native", so a pair with one empty
        cell came out tier 4 "prep mismatch inside the pair (native vs WGA)" —
        excluded from the default test, asserting a confound nobody had claimed.
        Unknown is graded on the tissue and SAYS it is unknown.
        """
        for t_prep, h_prep in ((None, None), ("native", ""), ("", "native")):
            tier, reason = bp.assign_tier(t_prep, h_prep, "foot")
            self.assertEqual(tier, 1, f"{t_prep!r}/{h_prep!r} did not reach tier 1")
            self.assertIn("not stated", reason)
            self.assertNotIn("mismatch", reason)

    def test_a_real_mismatch_is_still_tier_four(self):
        """The fix must not blunt the check it sits next to."""
        tier, reason = bp.assign_tier("native", "wga", "foot")
        self.assertEqual(tier, 4)
        self.assertIn("mismatch", reason)


class TestNoDesignStated(unittest.TestCase):
    def test_pairs_still_form_without_dna_prep_or_tissue(self):
        """
        dna_prep and tissue GRADE a pair; they do not create one. Requiring them
        turned an ordinary cohort — one that says who is paired with whom and
        nothing else — into "no design", which is a different and false claim.
        """
        tables = metadata([
            {"sample_code": "T1", "disease_status": "tumor", "individual": "a1"},
            {"sample_code": "H1", "disease_status": "matched host",
             "individual": "a1"},
        ])
        code, out, outdir = run(tables)
        self.assertEqual(code, 0, out)
        self.assertNotIn("NOT RUN", out)
        rows = (outdir / "pairs.tsv").read_text().strip().split("\n")
        self.assertEqual(len(rows), 2, "expected one pair")
        self.assertIn("not stated", rows[1])

    def test_absent_individual_is_reported_not_crashed(self):
        tables = metadata([
            {"sample_code": "S1", "disease_status": "tumor", "dna_prep": "native",
             "tissue": "foot"},
            {"sample_code": "S2", "disease_status": "healthy", "dna_prep": "native",
             "tissue": "foot"},
        ])
        code, out, outdir = run(tables)
        self.assertEqual(code, 0, out)
        self.assertIn("NOT RUN", out)
        self.assertIn("individual", out)
        # The empty table still has its header, so anything reading it works.
        self.assertTrue((outdir / "pairs.tsv").exists())
        self.assertEqual(
            (outdir / "pairs.tsv").read_text().strip().split("\n"),
            ["individual\ttumour_code\thost_code\ttier\ttier_reason"],
        )

    def test_a_header_with_no_values_counts_as_absent(self):
        """
        A samplesheet with an `individual` column and nothing in it describes no
        pairing. Reporting "the column is present" would be true and useless.
        """
        tables = metadata([
            {"sample_code": "S1", "disease_status": "tumor", "dna_prep": "native",
             "tissue": "foot", "individual": ""},
            {"sample_code": "S2", "disease_status": "matched host",
             "dna_prep": "native", "tissue": "foot", "individual": ""},
        ])
        code, out, _ = run(tables)
        self.assertEqual(code, 0, out)
        self.assertIn("NOT RUN", out)

    def test_the_report_says_how_to_supply_a_design(self):
        tables = metadata([{"sample_code": "S1"}])
        _code, out, _ = run(tables)
        self.assertIn("--input", out)
        self.assertIn("sample,fastq_1,fastq_2,individual", out)


class TestUnknownDiseaseStatus(unittest.TestCase):
    def test_a_typo_is_named_rather_than_silently_pairing_nothing(self):
        """
        An unrecognised value pairs nothing, which looks exactly like a cohort
        with no tumours. It must be named. It must NOT exit non-zero: that fails
        BUILD_PAIRS and makes Nextflow throw away everything AGGREGATE_COUNTS
        produced, over a spelling.
        """
        tables = metadata([
            {"sample_code": "S1", "disease_status": "tumuor", "dna_prep": "native",
             "tissue": "foot", "individual": "a1"},
            {"sample_code": "S2", "disease_status": "hoost",
             "dna_prep": "native", "tissue": "foot", "individual": "a1"},
        ])
        code, out, outdir = run(tables)
        self.assertEqual(code, 0, out)
        self.assertIn("tumuor", out)
        self.assertIn("unrecognised", out)
        self.assertTrue((outdir / "pairs.tsv").exists())
        self.assertTrue((outdir / "pairs_report.txt").exists())

    def test_a_typo_beside_good_data_does_not_lose_the_good_pairs(self):
        tables = metadata([
            {"sample_code": "T1", "disease_status": "tumor", "dna_prep": "native",
             "tissue": "blood", "individual": "a1"},
            {"sample_code": "H1", "disease_status": "matched host",
             "dna_prep": "native", "tissue": "foot", "individual": "a1"},
            {"sample_code": "X1", "disease_status": "wat", "dna_prep": "native",
             "tissue": "foot", "individual": "a2"},
        ])
        code, out, outdir = run(tables)
        self.assertEqual(code, 0, out)
        rows = (outdir / "pairs.tsv").read_text().strip().split("\n")
        self.assertEqual(len(rows), 2, "the valid pair was lost")

    def test_case_and_spelling_variants_are_accepted(self):
        """tumour/tumor and any capitalisation are the same word."""
        tables = metadata([
            {"sample_code": "S1", "disease_status": "TUMOUR", "dna_prep": "native",
             "tissue": "foot", "individual": "a1"},
            {"sample_code": "S2", "disease_status": "Matched Host",
             "dna_prep": "native", "tissue": "foot", "individual": "a1"},
        ])
        code, out, outdir = run(tables)
        self.assertEqual(code, 0, out)
        rows = (outdir / "pairs.tsv").read_text().strip().split("\n")
        self.assertEqual(len(rows), 2, "expected exactly one pair")
        self.assertIn("\t1\t", rows[1])  # tier 1: both native, host foot


class TestAGeneralCohortPairsCorrectly(unittest.TestCase):
    def test_words_not_codes_produce_the_right_tiers(self):
        rows = []
        for i, (tissue, prep, tier) in enumerate(
            [("foot", "native", 1), ("liver", "native", 2), ("foot", "wga", 3)], 1
        ):
            rows += [
                {"sample_code": f"T{i}", "disease_status": "tumor",
                 "dna_prep": prep, "tissue": "blood", "individual": f"patient_{i}"},
                {"sample_code": f"H{i}", "disease_status": "matched host",
                 "dna_prep": prep, "tissue": tissue, "individual": f"patient_{i}"},
            ]
            del tier
        code, out, outdir = run(metadata(rows))
        self.assertEqual(code, 0, out)
        got = {}
        for line in (outdir / "pairs.tsv").read_text().strip().split("\n")[1:]:
            f = line.split("\t")
            got[f[0]] = int(f[-2])
        self.assertEqual(got, {"patient_1": 1, "patient_2": 2, "patient_3": 3})


if __name__ == "__main__":
    unittest.main(verbosity=2)
