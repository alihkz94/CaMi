#!/usr/bin/env python3
"""
=============================================================================
test_ena_metadata.py — the sample-code parser, against every real edge case
-----------------------------------------------------------------------------
WHY THIS MATTERS MORE THAN A NORMAL PARSER TEST
    The ENA workbook has NO tissue column. Tissue is derived from the sample
    code and from nothing else, so no reconciliation report can catch a parser
    regression here — a wrong tissue letter silently moves a sample into the
    wrong arm of a tissue-matched comparison and nothing downstream complains.

    There is a precedent: stripping only "-wga" and requiring a
    bare trailing capital dropped a real pair without any error.

THE SUFFIXES MUST COMPOSE, NOT COMPETE
    Each form below is easy to handle alone. The tests that matter are the ones
    where two apply at once — a prep suffix on top of a replicate digit, a
    coverage annotation on top of a tissue letter — because that is where a fix
    for one form short-circuits the other.

RUN
    envs/stats/bin/python -m unittest discover -s cockle-microbiome/tests -v
=============================================================================
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "download"))

try:
    from build_ena_metadata import split_code  # noqa: E402
except SystemExit as exc:  # the script exits when openpyxl is absent
    raise unittest.SkipTest(str(exc)) from exc


class TestSampleCodeGrammar(unittest.TestCase):
    def assert_code(self, code, *, individual, tissue, replicate="",
                    prep="", hc=False):
        got = split_code(code)
        self.assertEqual(got["individual"], individual, f"{code}: individual")
        self.assertEqual(got["tissue_code"], tissue, f"{code}: tissue")
        self.assertEqual(got["replicate"], replicate, f"{code}: replicate")
        self.assertEqual(got["prep_suffix"], prep, f"{code}: prep suffix")
        self.assertEqual(got["is_high_coverage"], hc, f"{code}: _HC flag")

    # --- the plain forms ---------------------------------------------------
    def test_plain_tissue_letter(self):
        self.assert_code("EPCE18_851H", individual="EPCE18_851", tissue="H")
        self.assert_code("EPCE18_851M", individual="EPCE18_851", tissue="M")

    def test_prep_suffix_hides_the_tissue_letter(self):
        self.assert_code("FRCE17_840F-wga", individual="FRCE17_840",
                         tissue="F", prep="wga")

    def test_a_second_prep_suffix_behaves_the_same(self):
        """-gp is not -wga, and a rule that hard-codes -wga misses it."""
        self.assert_code("ENCE17_321H-gp", individual="ENCE17_321",
                         tissue="H", prep="gp")

    def test_trailing_replicate_digit(self):
        """This form silently dropped a real pair once."""
        self.assert_code("PACE17_421H1", individual="PACE17_421",
                         tissue="H", replicate="1")

    def test_no_tissue_letter_is_reported_not_guessed(self):
        got = split_code("ASCE17_1983")
        self.assertEqual(got["tissue_code"], "none")
        self.assertEqual(got["tissue_name"], "none")

    # --- the _HC trap this parser was written to close ---------------------
    def test_high_coverage_samples_are_foot_not_an_unknown_C_tissue(self):
        """
        Read naively, ICCE19_359F_HC ends in C — and C is not a tissue at all
        in the Bruzos codebook, it is a country letter. All five _HC samples
        are FOOT sequenced deeper. Treating them as unknown tissue would drop
        five healthy foot samples out of the tissue-matched comparison.
        """
        for code in ("ICCE19_359F_HC", "ICCE19_360F_HC", "ICCE19_363F_HC",
                     "ICCE19_366F_HC", "ICCE19_368F_HC"):
            self.assert_code(code, individual=code[:-4], tissue="F", hc=True)
            self.assertEqual(split_code(code)["tissue_name"], "foot")

    def test_C_is_never_accepted_as_a_tissue(self):
        got = split_code("ENCE17_999C")
        self.assertEqual(got["tissue_code"], "none")

    # --- the compositions: two suffixes at once ---------------------------
    def test_high_coverage_on_top_of_a_replicate_digit(self):
        self.assert_code("ICCE19_359F2_HC", individual="ICCE19_359",
                         tissue="F", replicate="2", hc=True)

    def test_prep_suffix_on_top_of_a_replicate_digit(self):
        self.assert_code("PACE17_421H1-wga", individual="PACE17_421",
                         tissue="H", replicate="1", prep="wga")

    def test_high_coverage_and_prep_suffix_together(self):
        self.assert_code("ICCE19_359F-wga_HC", individual="ICCE19_359",
                         tissue="F", prep="wga", hc=True)

    # --- the geography the haemolymph comparison depends on ---------------
    def test_decodes_country_location_species_and_year(self):
        got = split_code("ENCE16_224H")
        self.assertEqual(got["country_code"], "E")
        self.assertEqual(got["country_name"], "Spain")
        self.assertEqual(got["location_code"], "EN")
        self.assertEqual(got["location_name"], "Noia")
        self.assertEqual(got["species_code"], "CE")
        self.assertEqual(got["species_name"], "Cerastoderma edule")
        self.assertEqual(got["sampling_year"], 2016)

    def test_decodes_the_two_galician_sites_the_haemolymph_arm_comes_from(self):
        self.assertEqual(split_code("ECCE17_1H")["location_name"], "Carril")
        self.assertEqual(split_code("ELCE17_1H")["location_name"], "Placeres")

    def test_an_unrecognisable_prefix_yields_blanks_not_a_crash(self):
        got = split_code("BNg14")
        self.assertEqual(got["country_code"], "")
        self.assertEqual(got["location_name"], "")
        self.assertIsNone(got["sampling_year"])


class TestBuiltManifest(unittest.TestCase):
    """Checks against the committed manifest, if it has been generated."""

    @classmethod
    def setUpClass(cls):
        path = (Path(__file__).resolve().parent.parent.parent
                / "manifests" / "ena_metadata.tsv")
        if not path.exists():
            # A class-level skip is reported by unittest as ONE skip however
            # many tests it covers, so say the number. These check PRJEB58149's
            # own tallies against a file built from a submission workbook that
            # is not in this repository and never will be; they cannot run in
            # CI, and that is not the same as passing.
            raise unittest.SkipTest(
                f"7 cohort checks skipped: {path} does not exist. Build it with "
                "scripts/download/build_ena_metadata.py (PRJEB58149 only)."
            )
        import polars as pl
        cls.table = pl.read_csv(path, separator="\t", null_values=["NA"])

    def test_every_sample_is_present_exactly_once(self):
        self.assertEqual(self.table.height, 563)
        self.assertEqual(self.table["sample_code"].n_unique(), 563)

    def test_the_healthy_arithmetic_adds_up_without_double_counting(self):
        """410 foot (the 5 _HC included ONCE) + 17 + 6 + 29 = 462."""
        import polars as pl
        healthy = self.table.filter(pl.col("disease_status_ena") == "healthy")
        counts = dict(
            healthy.group_by("tissue_code").len()
            .select("tissue_code", "len").iter_rows()
        )
        self.assertEqual(counts, {"F": 410, "M": 17, "H": 6, "none": 29})
        self.assertEqual(sum(counts.values()), 462)
        # the _HC five are inside the 410, not additional to it
        self.assertEqual(healthy.filter(pl.col("is_high_coverage")).height, 5)

    def test_no_native_healthy_haemolymph_exists(self):
        """Every healthy haemolymph sample is WGA. A native one would mean the
        parser or the workbook is wrong, and would silently create a
        prep-matched comparison that this cohort cannot support."""
        import polars as pl
        native_h = self.table.filter(
            (pl.col("disease_status_ena") == "healthy")
            & (pl.col("tissue_code") == "H")
            & (pl.col("wga_status_ena") == "native")
        )
        self.assertEqual(native_h.height, 0)

    def test_the_foot_native_arms_are_the_sizes_the_analysis_assumes(self):
        import polars as pl
        foot_native = self.table.filter(
            (pl.col("tissue_code") == "F") & (pl.col("wga_status_ena") == "native")
        )
        self.assertEqual(
            foot_native.filter(
                pl.col("disease_status_ena") == "matched host").height, 16)
        self.assertEqual(
            foot_native.filter(
                pl.col("disease_status_ena") == "healthy").height, 410)

    def test_the_haemolymph_secondary_arms_are_18_versus_6(self):
        """Tumour WGA haemolymph at 15 Gb vs healthy WGA haemolymph at 15 Gb —
        the only contrast matched on tissue, prep AND depth at once."""
        import polars as pl
        wga_h_15 = self.table.filter(
            (pl.col("tissue_code") == "H")
            & (pl.col("wga_status_ena") == "wga")
            & (pl.col("sequencing_depth_gb") == 15.0)
        )
        self.assertEqual(
            wga_h_15.filter(pl.col("disease_status_ena") == "tumor").height, 18)
        self.assertEqual(
            wga_h_15.filter(pl.col("disease_status_ena") == "healthy").height, 6)

    def test_unavailable_fields_are_null_not_invented(self):
        for column in ("tumor_purity", "btn_lineage", "lane_id", "flowcell_id"):
            self.assertEqual(
                self.table[column].null_count(), self.table.height,
                f"{column} should be entirely NA — nothing has measured it",
            )

    def test_the_instrument_discrepancy_is_flagged_on_every_row(self):
        self.assertTrue(self.table["instrument_discrepancy"].all())


if __name__ == "__main__":
    unittest.main(verbosity=2)
