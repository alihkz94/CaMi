#!/usr/bin/env python3
"""confirm_draw.py — one sample's detections, their reads, and the BLAST draw.

Detections come from both classifiers, each under the same rule:
  1. Kraken 2 units: every genus and species the report detects (>= --min-reads
     reads and >= --min-frac of the genus-resolved reads). The per-read Kraken 2
     calls inside each unit's clade must equal the report's clade count, or the
     script exits 1.
  2. Kaiju units (with --kaiju and --taxdump): every genus and species Kaiju's own
     calls detect under the same rule (Bacteria, Archaea and Viruses only). A Kaiju
     unit that matches a Kraken 2 unit (same rank and taxid after both are resolved
     through the taxdump) is the same detection. The others are Kaiju-only: their
     Kaiju reads, and Kraken 2's call on each, go to <S>.kk2reads.tsv.gz.
  3. Up to --draw reads per unit (Kraken 2 and Kaiju-only), seeded with
     "<seed>|<sample>|<rank>|<taxid>" (Kaiju-only: "<seed>|kaiju|...") over the
     sorted read ids, so a rerun draws the same reads.
  4. Mate 1 of every drawn read, as FASTA. Every drawn read must be found once.

Usage : confirm_draw.py --sample S --report S.kraken2.report --calls S.calls.tsv.gz
                        --mate1 S_1.fastq.gz [--kaiju S.kaiju.out.gz --taxdump DIR]
                        [--draw 20 --seed 11 rule options]
Output: S.units.tsv      sample rank taxid name reads rel_abundance (Kraken 2 units)
        S.k2reads.tsv.gz read_id kraken_taxid (reads in any Kraken 2 unit clade)
        S.k2check.tsv    rank taxid report_clade_reads perread_clade_reads n_drawn
        S.kunits.tsv     sample rank taxid name reads rel_abundance kraken_unit (Kaiju units)
        S.kk2reads.tsv.gz read_id kaiju_taxid kraken_taxid (reads of Kaiju-only units;
                         kraken_taxid 0 = unclassified)
        S.draw.tsv       sample rank taxid read_id origin (kraken | kaiju)
        S.draw.fa        mate 1 of the drawn reads
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


def pick(ids, n, seed):
    ids = sorted(ids)
    return ids if len(ids) <= n else sorted(random.Random(seed).sample(ids, n))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--calls", required=True)
    p.add_argument("--mate1", required=True)
    p.add_argument("--kaiju", default="", help="Kaiju per-read output; adds Kaiju's own detections")
    p.add_argument("--taxdump", default="", help="NCBI taxdump directory (needed with --kaiju)")
    p.add_argument("--draw", type=int, default=20)
    p.add_argument("--seed", default="11")
    p.add_argument("--outdir", default=".")
    rule_args(p)
    a = p.parse_args()
    S, out = a.sample, a.outdir
    if a.kaiju and not a.taxdump:
        sys.exit("--kaiju needs --taxdump")
    rule = rule_from(a)

    rows, _ = cc.detections(a.report, rule)
    units = [(r[0], r[1], r[2], r[3], r[4]) for r in rows if r[5]]
    with open(f"{out}/{S}.units.tsv", "w") as fh:
        fh.write("sample\trank\ttaxid\tname\treads\trel_abundance\n")
        for rank, tid, name, n, ra in units:
            fh.write(f"{S}\t{rank}\t{tid}\t{name}\t{n}\t{ra:.6g}\n")

    # Kaiju's own detections, and which of them Kraken 2 does not make
    kunits, konly, kmembers, kaiju_of = [], [], {}, {}
    if a.kaiju:
        tx = cc.Taxonomy(a.taxdump)
        kraken_keys = {(rank, tid if tid == rule.merge_family else (tx.resolve(tid) or tid)) for rank, tid, *_ in units}
        calls = cc.kaiju_calls(a.kaiju)
        krows, _ = cc.kaiju_detections(calls, tx, rule)
        kunits = [(r[0], r[1], r[2], r[3], r[4], (r[0], r[1]) in kraken_keys) for r in krows if r[5]]
        konly = [u for u in kunits if not u[5]]
        want = {u[1] for u in konly}
        hit = {}
        for rid, t in calls.items():
            us = hit.get(t)
            if us is None:
                lin = set(tx.lineage(t))
                us = hit[t] = [(u[0], u[1]) for u in konly if u[1] in lin] if want & lin else []
            if us:
                kaiju_of[rid] = t
                for u in us:
                    kmembers.setdefault(u, []).append(rid)
        del calls
    with open(f"{out}/{S}.kunits.tsv", "w") as fh:
        fh.write("sample\trank\ttaxid\tname\treads\trel_abundance\tkraken_unit\n")
        for rank, tid, name, n, ra, both in kunits:
            fh.write(f"{S}\t{rank}\t{tid}\t{name}\t{n}\t{ra:.6g}\t{'yes' if both else 'no'}\n")

    clades = cc.clade_taxids(a.report, [u[1] for u in units])
    tax2units = {}
    for rank, tid, *_ in units:
        for t in clades[tid]:
            tax2units.setdefault(t, []).append((rank, tid))
    members = {(u[0], u[1]): [] for u in units}
    kraken_of = {}
    with cc.open_text(a.calls) as fh, gzip.open(f"{out}/{S}.k2reads.tsv.gz", "wt", compresslevel=6) as kr:
        for line in fh:
            rid, tid = line.rstrip("\n").split("\t")
            rid = cc.read_id(rid)
            if rid in kaiju_of:
                kraken_of[rid] = tid
            us = tax2units.get(tid)
            if us:
                kr.write(f"{rid}\t{tid}\n")
                for u in us:
                    members[u].append(rid)
    with gzip.open(f"{out}/{S}.kk2reads.tsv.gz", "wt", compresslevel=6) as fh:
        for rid in sorted(kaiju_of):
            fh.write(f"{rid}\t{kaiju_of[rid]}\t{kraken_of.get(rid, '0')}\n")

    bad, draw = [], {}
    with open(f"{out}/{S}.k2check.tsv", "w") as ck:
        ck.write("rank\ttaxid\treport_clade_reads\tperread_clade_reads\tn_drawn\n")
        for rank, tid, _name, n_rep, _ra in units:
            ids = members[(rank, tid)]
            draw[("kraken", rank, tid)] = pick(ids, a.draw, f"{a.seed}|{S}|{rank}|{tid}")
            ck.write(f"{rank}\t{tid}\t{n_rep}\t{len(ids)}\t{len(draw[('kraken', rank, tid)])}\n")
            if len(ids) != n_rep:
                bad.append(f"{rank}:{tid} report={n_rep} per-read={len(ids)}")
    if bad:
        sys.exit(f"{S}: per-read clade counts differ from the report: {bad[:5]}")
    for rank, tid, *_ in konly:
        draw[("kaiju", rank, tid)] = pick(kmembers.get((rank, tid), []), a.draw, f"{a.seed}|kaiju|{S}|{rank}|{tid}")

    with open(f"{out}/{S}.draw.tsv", "w") as mp:
        mp.write("sample\trank\ttaxid\tread_id\torigin\n")
        for (origin, rank, tid), ids in draw.items():
            for r in ids:
                mp.write(f"{S}\t{rank}\t{tid}\t{r}\t{origin}\n")

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
    print(f"{S}: {len(units)} Kraken 2 units, {len(konly)} Kaiju-only units, {len(want)} reads drawn")


if __name__ == "__main__":
    main()
