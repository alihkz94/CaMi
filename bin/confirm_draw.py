#!/usr/bin/env python3
"""confirm_draw.py — one sample's detections, their Kraken 2 reads, and the BLAST draw.

  1. Units: every genus and species the report detects under the rule
     (>= --min-reads reads and >= --min-frac of the genus-resolved reads).
  2. The per-read Kraken 2 calls inside each unit's clade. Gate: the per-read
     count of every unit equals the report's clade count, or the script exits 1.
  3. Up to --draw reads per unit, seeded with "<seed>|<sample>|<rank>|<taxid>"
     over the sorted read ids, so a rerun draws the same reads.
  4. Mate 1 of every drawn read, as FASTA. Every drawn read must be found once.

Usage : confirm_draw.py --sample S --report S.kraken2.report --calls S.calls.tsv.gz
                        --mate1 S_1.fastq.gz [--draw 20 --seed 11 rule options]
Output: S.units.tsv     sample rank taxid name reads rel_abundance
        S.k2reads.tsv.gz read_id kraken_taxid (reads in any unit clade)
        S.k2check.tsv   rank taxid report_clade_reads perread_clade_reads n_drawn
        S.draw.tsv      sample rank taxid read_id
        S.draw.fa       mate 1 of the drawn reads
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import confirm_common as cc  # noqa: E402


def rule_args(p):
    p.add_argument("--min-reads", type=int, default=10)
    p.add_argument("--min-frac", type=float, default=0.01)
    p.add_argument("--exclude", default="9606,9605", help="taxids never counted as a detection")
    p.add_argument("--merge-family", default="", help="family taxid whose genera count as one entry")


def rule_from(a):
    return cc.Rule(min_reads=a.min_reads, min_frac=a.min_frac,
                   exclude=frozenset(t for t in a.exclude.split(",") if t),
                   merge_family=a.merge_family or None)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--calls", required=True)
    p.add_argument("--mate1", required=True)
    p.add_argument("--draw", type=int, default=20)
    p.add_argument("--seed", default="11")
    p.add_argument("--outdir", default=".")
    rule_args(p)
    a = p.parse_args()
    S, out = a.sample, a.outdir

    rows, _ = cc.detections(a.report, rule_from(a))
    units = [(r[0], r[1], r[2], r[3], r[4]) for r in rows if r[5]]
    with open(f"{out}/{S}.units.tsv", "w") as fh:
        fh.write("sample\trank\ttaxid\tname\treads\trel_abundance\n")
        for rank, tid, name, n, ra in units:
            fh.write(f"{S}\t{rank}\t{tid}\t{name}\t{n}\t{ra:.6g}\n")

    clades = cc.clade_taxids(a.report, [u[1] for u in units])
    tax2units = {}
    for rank, tid, *_ in units:
        for t in clades[tid]:
            tax2units.setdefault(t, []).append((rank, tid))
    members = {(u[0], u[1]): [] for u in units}
    with cc.open_text(a.calls) as fh, gzip.open(f"{out}/{S}.k2reads.tsv.gz", "wt", compresslevel=6) as kr:
        for line in fh:
            rid, tid = line.rstrip("\n").split("\t")
            us = tax2units.get(tid)
            if us:
                rid = cc.read_id(rid)
                kr.write(f"{rid}\t{tid}\n")
                for u in us:
                    members[u].append(rid)

    bad, draw = [], {}
    with open(f"{out}/{S}.k2check.tsv", "w") as ck:
        ck.write("rank\ttaxid\treport_clade_reads\tperread_clade_reads\tn_drawn\n")
        for rank, tid, _name, n_rep, _ra in units:
            ids = sorted(members[(rank, tid)])
            rng = random.Random(f"{a.seed}|{S}|{rank}|{tid}")
            pick = ids if len(ids) <= a.draw else rng.sample(ids, a.draw)
            draw[(rank, tid)] = sorted(pick)
            ck.write(f"{rank}\t{tid}\t{n_rep}\t{len(ids)}\t{len(pick)}\n")
            if len(ids) != n_rep:
                bad.append(f"{rank}:{tid} report={n_rep} per-read={len(ids)}")
    if bad:
        sys.exit(f"{S}: per-read clade counts differ from the report: {bad[:5]}")

    with open(f"{out}/{S}.draw.tsv", "w") as mp:
        mp.write("sample\trank\ttaxid\tread_id\n")
        for (rank, tid), ids in draw.items():
            for r in ids:
                mp.write(f"{S}\t{rank}\t{tid}\t{r}\n")

    want = {r for ids in draw.values() for r in ids}
    found = {}
    if want:
        with cc.open_text(a.mate1) as fh:
            while True:
                head = fh.readline()
                if not head:
                    break
                seq = fh.readline().rstrip("\n")
                fh.readline()
                fh.readline()
                rid = cc.read_id(head[1:].rstrip("\n"))
                if rid in want:
                    if rid in found:
                        sys.exit(f"{S}: read {rid} occurs twice in {a.mate1}")
                    found[rid] = seq
    if len(found) != len(want):
        sys.exit(f"{S}: drew {len(want)} reads, found {len(found)} in {a.mate1}")
    with open(f"{out}/{S}.draw.fa", "w") as fh:
        for r in sorted(found):
            fh.write(f">{r}\n{found[r]}\n")
    print(f"{S}: {len(units)} units, {len(want)} reads drawn")


if __name__ == "__main__":
    main()
