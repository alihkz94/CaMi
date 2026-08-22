#!/usr/bin/env python3
"""
=============================================================================
test_paired_test_e2e.py — paired_test.py end to end, on planted ground truth
-----------------------------------------------------------------------------
WHY THIS EXISTS
    On the live cohort the paired test currently has 0 pairs at tier <= 2, so
    running it against real data exercises only the "nothing to test" path.
    That is not a regression test — the exact enumeration, the CLR transform,
    the prevalence filter and the FDR step are all untouched by it.

    This builds a small synthetic cohort with a KNOWN answer, runs the real
    script as a subprocess exactly as Nextflow does, and checks the recovered
    effect against what was planted.

WHAT IS PLANTED
    12 pairs. Three taxa are multiplied by 8 in every tumour arm. The rest are
    drawn from one distribution for both arms.

WHAT IS ASSERTED
    - the three planted taxa are recovered, in the right direction
    - the test really used exact enumeration (mode + permutation count)
    - The compositional constraint still holds: unrelated taxa may show significant
      NEGATIVE shifts, because CLR is compositional — if some taxa rise the
      rest must fall in relative terms. That is the constraint, not a false
      positive, and the test documents it rather than pretending otherwise.

RUN
    envs/stats/bin/python -m unittest discover -s cockle-microbiome/tests -v
=============================================================================
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import polars as pl

BIN = Path(__file__).resolve().parent.parent / "bin"
SCRIPT = BIN / "paired_test.py"

N_PAIRS = 12
N_TAXA = 40
PLANTED = [101, 102, 103]     # taxon ids given a real effect
EFFECT = 8.0


def build_cohort(tmp: Path) -> tuple[Path, Path]:
    """Write counts_long.parquet and pairs.tsv the way the STATS step does."""
    rng = np.random.default_rng(20260821)
    tables, pairs_dir = tmp / "tables", tmp / "pairs"
    tables.mkdir()
    pairs_dir.mkdir()

    taxon_ids = list(range(100, 100 + N_TAXA))
    nonhost = 1_000_000
    rows = []
    for i in range(N_PAIRS):
        base = rng.integers(200, 2000, size=N_TAXA).astype(float)
        for role, code in (("tumour", f"IND{i:02d}H"), ("host", f"IND{i:02d}F")):
            counts = base.copy()
            if role == "tumour":
                for tid in PLANTED:
                    counts[taxon_ids.index(tid)] *= EFFECT
            counts = counts * rng.uniform(0.9, 1.1, size=N_TAXA)  # sampling noise
            for tid, c in zip(taxon_ids, counts):
                rows.append(
                    {
                        "sample_code": code,
                        "method": "bracken",
                        "rank": "genus",
                        "taxon_id": tid,
                        "taxon_name": f"Taxon_{tid}",
                        "reads": int(round(c)),
                        "is_host_artifact": False,
                        "nonhost_pairs": nonhost,
                        "frac_of_nonhost": c / nonhost,
                    }
                )
    pl.DataFrame(rows).write_parquet(tables / "counts_long.parquet")

    pl.DataFrame(
        {
            "individual": [f"IND{i:02d}" for i in range(N_PAIRS)],
            "tumour_code": [f"IND{i:02d}H" for i in range(N_PAIRS)],
            "host_code": [f"IND{i:02d}F" for i in range(N_PAIRS)],
            "tier": [1] * N_PAIRS,
        }
    ).write_csv(pairs_dir / "pairs.tsv", separator="\t")
    return tables, pairs_dir


class TestPairedTestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        tables, pairs_dir = build_cohort(tmp)
        cls.outdir = tmp / "out"
        proc = subprocess.run(
            [
                sys.executable, str(SCRIPT),
                "--tables", str(tables),
                "--pairs", str(pairs_dir),
                "--outdir", str(cls.outdir),
                "--method", "bracken", "--rank", "genus",
                "--min-prevalence", "0.5",
            ],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            raise AssertionError(f"paired_test.py failed:\n{proc.stderr}")
        cls.res = pl.read_csv(
            cls.outdir / "paired_results_bracken_genus.tsv", separator="\t"
        )
        cls.report = (cls.outdir / "paired_test_report.txt").read_text()

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_it_used_exact_enumeration_not_an_approximation(self):
        modes = set(self.res["test_mode"].to_list())
        self.assertEqual(modes, {"exact"})
        self.assertEqual(
            set(self.res["permutations"].to_list()), {2 ** N_PAIRS}
        )

    def test_every_planted_taxon_is_recovered_in_the_right_direction(self):
        for tid in PLANTED:
            row = self.res.filter(pl.col("taxon_id") == tid)
            self.assertEqual(row.height, 1, f"taxon {tid} missing from results")
            r = row.row(0, named=True)
            self.assertGreater(
                r["mean_clr_diff"], 0.0,
                f"taxon {tid} planted UP in tumour but shifted down",
            )
            self.assertEqual(
                r["n_tumour_higher"], N_PAIRS,
                f"taxon {tid} should be higher in all {N_PAIRS} pairs",
            )
            self.assertLess(r["q_bh"], 0.05, f"taxon {tid} not significant")

    def test_the_p_value_floor_is_the_hard_paired_bound(self):
        floors = set(self.res["p_floor"].to_list())
        self.assertEqual(len(floors), 1)
        self.assertAlmostEqual(floors.pop(), 1.0 / (1 << N_PAIRS), places=15)

    def test_compositional_constraint_is_visible_not_hidden(self):
        """If some taxa rise, others must fall in CLR terms.

        Unrelated taxa showing significant NEGATIVE shifts is the compositional
        constraint doing its job. The test asserts the direction of that
        artifact so nobody later "fixes" it into a false-positive hunt.
        """
        others = self.res.filter(~pl.col("taxon_id").is_in(PLANTED))
        sig_down = others.filter(
            (pl.col("q_bh") < 0.05) & (pl.col("mean_clr_diff") < 0)
        )
        sig_up = others.filter(
            (pl.col("q_bh") < 0.05) & (pl.col("mean_clr_diff") > 0)
        )
        self.assertEqual(
            sig_up.height, 0,
            "an unplanted taxon rose significantly — that is a real false positive",
        )
        self.assertGreater(
            sig_down.height, 0,
            "expected the compositional constraint to push unrelated taxa down",
        )

    def test_report_states_the_ceiling_and_the_tissue_caveat(self):
        self.assertIn("DETECTABILITY CEILING", self.report)
        self.assertIn("STANDING CAVEAT", self.report)


if __name__ == "__main__":
    unittest.main(verbosity=2)
