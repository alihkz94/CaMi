#!/usr/bin/env python3
"""
=============================================================================
test_confirm.py — the triple confirmation on a cohort small enough to reason about
-----------------------------------------------------------------------------
One sample, one family, two genera:

  Genus A (100)       40 reads. Kaiju calls every read genus A, BLAST hits
                      species A1 at 99 %                          -> triple
  Species A1 (101)    20 reads. Kaiju resolves them only to genus A, so the
                      species leg fails and the genus fallback confirms it;
                      BLAST hits A1 under its RETIRED taxid 999, which
                      merged.dmp maps to 101                      -> triple
  Genus B (200)       20 reads. Kaiju says genus A (disagrees), BLAST finds
                      nothing (no_hit)                            -> single
  Homo (9605/9606)    30 reads, never a detection and not in the denominator

The same run on the cockle cohort (563 samples, 7,656 detections) reproduced
the published confirmation table; that check needs the data and lives in the
release procedure, not here.

RUN:  python3 -m unittest discover -s tests
=============================================================================
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / "bin"
S = "S1"

REPORT = """\
10.00\t10\t10\tU\t0\tunclassified
90.00\t90\t0\tR\t1\troot
60.00\t60\t0\tD\t2\t  Bacteria
60.00\t60\t0\tF\t10\t    Family X
40.00\t40\t20\tG\t100\t      Genus A
20.00\t20\t20\tS\t101\t        Genus A sp1
20.00\t20\t20\tG\t200\t      Genus B
30.00\t30\t0\tD\t2759\t  Eukaryota
30.00\t30\t0\tG\t9605\t    Homo
30.00\t30\t30\tS\t9606\t      Homo sapiens
"""

NODES = [("1", "1", "no rank"), ("2", "1", "superkingdom"), ("10", "2", "family"),
         ("100", "10", "genus"), ("101", "100", "species"), ("200", "10", "genus"),
         ("2759", "1", "superkingdom"), ("9605", "2759", "genus"), ("9606", "9605", "species")]
NAMES = {"1": "root", "2": "Bacteria", "10": "Family X", "100": "Genus A", "101": "Genus A sp1",
         "200": "Genus B", "2759": "Eukaryota", "9605": "Homo", "9606": "Homo sapiens"}

READS = ([(f"a{i:02d}", "100") for i in range(20)] + [(f"s{i:02d}", "101") for i in range(20)]
         + [(f"b{i:02d}", "200") for i in range(20)] + [(f"h{i:02d}", "9606") for i in range(30)])
UNCLASSIFIED = [f"u{i:02d}" for i in range(10)]
SEQ = "ACGT" * 37 + "AC"


def run(*cmd, cwd):
    return subprocess.run([sys.executable, *map(str, cmd)], cwd=cwd, capture_output=True, text=True)


class TripleConfirmation(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        tx = self.d / "taxdump"
        tx.mkdir()
        # The real file has 13 columns; the rank must not be the last one.
        (tx / "nodes.dmp").write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n" for t, p, r in NODES))
        (tx / "names.dmp").write_text("".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in NAMES.items()))
        (tx / "merged.dmp").write_text("999\t|\t101\t|\n")
        (self.d / f"{S}.kraken2.report").write_text(REPORT)
        with gzip.open(self.d / f"{S}.kraken2.calls.tsv.gz", "wt") as fh:
            fh.write("".join(f"{r}/1\t{t}\n" for r, t in READS))
        with gzip.open(self.d / f"{S}_1.fastq.gz", "wt") as fh:
            for r in [r for r, _ in READS] + UNCLASSIFIED:
                fh.write(f"@{r}/1\n{SEQ}\n+\n{'I' * len(SEQ)}\n")
        (self.d / f"{S}.pairs").write_text(f"{len(READS) + len(UNCLASSIFIED)}\n")
        with gzip.open(self.d / f"{S}.kaiju.out.gz", "wt") as fh:
            for r, t in READS:
                call = "0" if t == "9606" else "100"
                fh.write(f"{'U' if call == '0' else 'C'}\t{r}\t{call}\n")
            fh.write("".join(f"U\t{r}\t0\n" for r in UNCLASSIFIED))

    def tearDown(self):
        self.tmp.cleanup()

    def draw(self):
        return run(BIN / "confirm_draw.py", "--sample", S, "--report", f"{S}.kraken2.report",
                   "--calls", f"{S}.kraken2.calls.tsv.gz", "--mate1", f"{S}_1.fastq.gz", cwd=self.d)

    def blast(self):
        r = run(BIN / "confirm_blast_prep.py", "--chunk", "25", cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(self.d / "queries.tsv") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        by_chunk = {}
        for q in rows:
            if q["read_id"].startswith(("a", "s")):
                tax = "999" if q["read_id"].startswith("s") else "101"
                by_chunk.setdefault(q["chunk"], []).append(f"{q['qid']}\tACC1\t{tax}\t99.0\t150\t150\t1e-70\t270\n")
            else:
                by_chunk.setdefault(q["chunk"], [])
        for c, lines in by_chunk.items():
            with gzip.open(self.d / f"{c}.blast.tsv.gz", "wt") as fh:
                fh.write("".join(lines))
        return len(by_chunk)

    def combine(self):
        (self.d / "samples.txt").write_text(f"{S}\n")
        return run(BIN / "confirm_combine.py", "--samples", "samples.txt", "--taxdump", "taxdump",
                   "--jobs", "1", cwd=self.d)

    def table(self):
        with open(self.d / "confirmation_long.tsv") as fh:
            return {(r["rank"], r["taxid"]): r for r in csv.DictReader(fh, delimiter="\t")}

    def test_detections_exclude_human_and_use_genus_denominator(self):
        self.assertEqual(self.draw().returncode, 0)
        with open(self.d / f"{S}.units.tsv") as fh:
            units = {(r["rank"], r["taxid"]): float(r["rel_abundance"]) for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(set(units), {("genus", "100"), ("species", "101"), ("genus", "200")})
        self.assertAlmostEqual(units[("genus", "100")], 40 / 60, places=5)

    def test_draw_is_capped_seeded_and_strips_mate_suffix(self):
        self.assertEqual(self.draw().returncode, 0)
        first = (self.d / f"{S}.draw.tsv").read_text()
        self.assertEqual(self.draw().returncode, 0)
        self.assertEqual(first, (self.d / f"{S}.draw.tsv").read_text())
        genus_a = [l for l in first.splitlines() if "\tgenus\t100\t" in l]
        self.assertEqual(len(genus_a), 20)
        self.assertNotIn("/1", first)
        self.assertTrue((self.d / f"{S}.draw.fa").read_text().startswith(">"))

    def test_clade_gate_refuses_incomplete_calls(self):
        with gzip.open(self.d / f"{S}.kraken2.calls.tsv.gz", "wt") as fh:
            fh.write("".join(f"{r}\t{t}\n" for r, t in READS[1:]))
        r = self.draw()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("differ from the report", r.stderr)

    def test_tiers(self):
        self.assertEqual(self.draw().returncode, 0)
        self.assertGreater(self.blast(), 1)
        r = self.combine()
        self.assertEqual(r.returncode, 0, r.stderr)
        t = self.table()
        self.assertEqual(t[("genus", "100")]["tier"], "triple")
        self.assertEqual(t[("genus", "100")]["kaiju_status"], "confirmed")
        sp = t[("species", "101")]
        self.assertEqual((sp["tier"], sp["kaiju_status"], sp["kaiju_level"], sp["blast_status"]),
                         ("triple", "confirmed_at_genus", "genus", "confirmed"))
        self.assertEqual(sp["genus_unit_tier"], "triple")
        b = t[("genus", "200")]
        self.assertEqual((b["tier"], b["kaiju_status"], b["blast_status"]), ("single", "disagrees", "no_hit"))

    def test_refuses_when_kaiju_and_kraken_saw_different_pairs(self):
        self.assertEqual(self.draw().returncode, 0)
        self.blast()
        (self.d / f"{S}.pairs").write_text("99\n")
        r = self.combine()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("REFUSING", r.stderr)

    def test_refuses_when_a_blast_chunk_is_missing(self):
        self.assertEqual(self.draw().returncode, 0)
        self.blast()
        next(self.d.glob("chunk_*.blast.tsv.gz")).unlink()
        r = self.combine()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("chunk results missing", r.stderr)


if __name__ == "__main__":
    unittest.main()
