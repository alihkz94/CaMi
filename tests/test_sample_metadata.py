#!/usr/bin/env python3
"""
=============================================================================
test_sample_metadata.py — where a sample's design comes from
-----------------------------------------------------------------------------
Two sources produce the same table, and each has a way of going wrong silently:

  the samplesheet   a column the person did not fill must stay unknown, not
                    become a guess, and a header with no rows must not read as
                    "a cohort with no samples"
  the sample code   the Bruzos grammar, which has now had three separate bugs,
                    every one of them silent

THE GRAMMAR IS LOAD-BEARING AND UNCHECKABLE. The ENA submission workbook has no
tissue column, so tissue is derived from the sample code and from nothing else.
A regression here moves a sample into the wrong arm of a comparison and no other
file disagrees. That is why the suffix forms are tested IN COMBINATION: a fix
for one suffix can short-circuit another, which is exactly how _HC broke the
replicate digit's predecessor.

RUN:  python3 -m unittest discover -s tests
=============================================================================
"""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bin"))

import aggregate_counts as ac  # noqa: E402


class TestSampleCodeGrammar(unittest.TestCase):
    """Every suffix form in the cohort, alone and combined."""

    def parts(self, code):
        p = ac.split_code(code)
        return p.individual, p.tissue, p.replicate

    def test_plain_code(self):
        self.assertEqual(self.parts("EPCE18_851H"), ("EPCE18_851", "H", ""))

    def test_prep_suffix(self):
        self.assertEqual(self.parts("FRCE17_840F-wga"), ("FRCE17_840", "F", ""))
        self.assertEqual(self.parts("ENCE17_321H-gp"), ("ENCE17_321", "H", ""))

    def test_replicate_digit(self):
        """This form silently dropped the PACE17_421 pair."""
        self.assertEqual(self.parts("PACE17_421H1"), ("PACE17_421", "H", "1"))

    def test_high_coverage_is_not_a_tissue(self):
        """
        Read naively the code ends in C, and C is a COUNTRY
        letter, not a tissue. Five healthy FOOT samples were being recorded as an
        invented tissue and dropped out of every foot comparison.
        """
        ind, tissue, _rep = self.parts("ICCE19_359F_HC")
        self.assertEqual(tissue, "F")
        self.assertEqual(ind, "ICCE19_359")
        self.assertEqual(ac.TISSUE_NAMES[tissue], "foot")

    def test_all_five_high_coverage_samples(self):
        for n in (359, 360, 363, 366, 368):
            self.assertEqual(self.parts(f"ICCE19_{n}F_HC"), (f"ICCE19_{n}", "F", ""))

    def test_suffixes_compose_and_do_not_short_circuit(self):
        """
        The point of this test. Stripping _HC must happen BEFORE the prep suffix
        and BEFORE the tissue letter, and must leave the replicate group able to
        fire. A fix that replaces the grammar instead of layering in front of it
        passes the two tests above and fails these.
        """
        self.assertEqual(self.parts("ICCE19_359F-wga_HC"), ("ICCE19_359", "F", ""))
        self.assertEqual(self.parts("PACE17_421H1_HC"), ("PACE17_421", "H", "1"))
        self.assertEqual(self.parts("PACE17_421H1-wga"), ("PACE17_421", "H", "1"))

    def test_no_tissue_letter_is_recorded_honestly(self):
        """A code with no tissue cannot be paired. Say so; do not force it."""
        self.assertEqual(self.parts("ASCE17_1983"), ("ASCE17_1983", "none", ""))

    def test_c_is_not_in_the_tissue_table(self):
        """The bogus entry that made the _HC bug invisible."""
        self.assertNotIn("C", ac.TISSUE_NAMES)

    def test_adductor_is_resolved_not_hedged(self):
        self.assertEqual(ac.TISSUE_NAMES["A"], "adductor muscle")

    def test_both_codebook_letters_for_foot(self):
        self.assertEqual(ac.TISSUE_NAMES["F"], ac.TISSUE_NAMES["P"])


class TestTissueLabel(unittest.TestCase):
    def test_a_code_is_expanded(self):
        self.assertEqual(ac.tissue_label("H"), "haemolymph")

    def test_a_name_passes_through(self):
        """Another cohort writes words. Do not invent a code for them."""
        for name in ("foot", "whole animal", "gill tissue", "Liver"):
            self.assertEqual(ac.tissue_label(name), name)

    def test_an_unknown_capital_is_kept_as_itself(self):
        self.assertEqual(ac.tissue_label("Z"), "Z")


class TestLoadDesign(unittest.TestCase):
    """The samplesheet route."""

    def write(self, text):
        """
        Written TAB-separated, which is what Nextflow emits. It used to be
        commas, and a tab-separated samplesheet may legally hold a value with a
        comma in it — which then produced a row with more fields than the header
        and an error naming this file instead of the sheet that caused it.
        """
        tmp = Path(tempfile.mkdtemp()) / "design.tsv"
        tmp.write_text(text.replace(",", "\t"))
        return tmp

    def test_header_only_means_no_samplesheet(self):
        """
        This is what Nextflow writes on the ENA route. It must read as "fall back
        to the manifests", NOT as "a cohort with zero samples" — which would
        replace a working design with nothing.
        """
        p = self.write("sample,individual,tissue,disease_status,dna_prep\n")
        self.assertIsNone(ac.load_design(p, "none"))

    def test_missing_file_means_no_samplesheet(self):
        self.assertIsNone(ac.load_design(Path("/nonexistent/design.csv"), "none"))

    def test_full_design_is_read_as_given(self):
        p = self.write(
            "sample,individual,tissue,disease_status,dna_prep\n"
            "T01,animal_01,haemolymph,tumor,native\n"
            "H01,animal_01,foot,matched host,native\n"
        )
        df = ac.load_design(p, "none")
        self.assertEqual(df.height, 2)
        self.assertEqual(df["individual"].to_list(), ["animal_01", "animal_01"])
        self.assertEqual(df["tissue_name"].to_list(), ["haemolymph", "foot"])
        self.assertEqual(df["disease_status"].to_list(), ["tumor", "matched host"])
        self.assertEqual(df["source"].to_list(), ["SAMPLESHEET", "SAMPLESHEET"])

    def test_absent_columns_stay_unknown(self):
        """Reads only. No design stated, so nothing may be invented."""
        p = self.write("sample\nS1\nS2\n")
        df = ac.load_design(p, "none")
        self.assertEqual(df.height, 2)
        for field in ("individual", "disease_status", "dna_prep"):
            self.assertTrue(all(v is None for v in df[field].to_list()),
                            f"{field} was invented")
        self.assertTrue(all(v is None for v in df["tissue"].to_list()),
                        'tissue must stay null, not become the string "none" — '
                        "a column full of a value reads as stated")

    def test_scheme_none_never_parses_a_name(self):
        """
        Another cohort's names must not be read as tissue codes. 'Patient1A'
        ends in a capital, which the grammar would happily call adductor muscle.
        """
        p = self.write("sample\nPatient1A\n")
        df = ac.load_design(p, "none")
        self.assertIsNone(df["individual"][0])
        self.assertIsNone(df["tissue"][0])

    def test_scheme_bruzos_fills_only_what_is_blank(self):
        """The grammar is a fallback, never an override."""
        p = self.write(
            "sample,individual,tissue,disease_status,dna_prep\n"
            "EPCE18_851H,,,tumor,native\n"
            "EPCE18_851M,stated_by_hand,mantle,matched host,native\n"
        )
        df = ac.load_design(p, "bruzos")
        self.assertEqual(df["individual"].to_list(), ["EPCE18_851", "stated_by_hand"])
        self.assertEqual(df["tissue"].to_list(), ["H", "mantle"])
        self.assertEqual(df["tissue_name"].to_list(), ["haemolymph", "mantle"])

    def test_duplicate_samples_are_collapsed_not_doubled(self):
        p = self.write("sample,individual\nS1,a\nS1,b\n")
        df = ac.load_design(p, "none")
        self.assertEqual(df.height, 1)
        self.assertEqual(df["individual"][0], "a")

    def test_the_column_shape_matches_the_manifest_route(self):
        """
        Everything downstream reads one shape. If these drift apart, a cohort
        works on one route and fails on the other.
        """
        p = self.write("sample\nS1\n")
        df = ac.load_design(p, "none")
        self.assertEqual(
            df.columns,
            ["run_accession", "sample_code", "disease_status", "dna_prep",
             "source", "individual", "tissue", "tissue_name", "site", "replicate"],
        )


class TestEverythingIsText(unittest.TestCase):
    """
    Sample names and animal IDs are LABELS. Letting polars infer their type
    turns "001" into the integer 1 — and that name is what every output file is
    called, so each sample then looks ABSENT and the run dies blaming the
    profiling instead of the sheet.
    """

    def write(self, text):
        tmp = Path(tempfile.mkdtemp()) / "design.tsv"
        tmp.write_text(text.replace(",", "\t"))
        return tmp

    def test_leading_zero_sample_names_survive(self):
        p = self.write("sample,individual\n001,1\n002,1\n010,2\n")
        df = ac.load_design(p, "none")
        self.assertEqual(df["sample_code"].to_list(), ["001", "002", "010"])

    def test_numeric_individual_stays_a_string(self):
        """Comparing an i64 column against "" raises from somewhere unrelated."""
        p = self.write("sample,individual\nS1,1\nS2,1\n")
        df = ac.load_design(p, "none")
        self.assertEqual(df["individual"].to_list(), ["1", "1"])
        self.assertEqual(ac.stated_values(df, "individual"), True)

    def test_a_numeric_tissue_does_not_crash_the_label(self):
        p = self.write("sample,tissue\nS1,1\n")
        self.assertEqual(ac.load_design(p, "none")["tissue_name"][0], "1")


class TestStatedValues(unittest.TestCase):
    """The one question three different places used to answer three ways."""

    def frame(self, **cols):
        import polars as pl
        return pl.DataFrame(cols)

    def test_absent_column(self):
        self.assertFalse(ac.stated_values(self.frame(a=["x"]), "b"))

    def test_all_null_column(self):
        self.assertFalse(ac.stated_values(self.frame(a=[None, None]), "a"))

    def test_all_blank_column(self):
        self.assertFalse(ac.stated_values(self.frame(a=["", "  "]), "a"))

    def test_one_value_is_enough(self):
        self.assertTrue(ac.stated_values(self.frame(a=["", "x"]), "a"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
