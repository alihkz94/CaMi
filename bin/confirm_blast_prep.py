#!/usr/bin/env python3
"""confirm_blast_prep.py — one BLAST query set for the cohort, in fixed-size chunks.

Reads every <S>.draw.tsv and <S>.draw.fa in the working directory.
  1. Distinct (sample, read) pairs -> query ids Q00000001.. (a read drawn for both
     a genus and a species unit is searched once).
  2. The queries are shuffled with a fixed seed across samples and cut into chunks
     of --chunk reads, so one sample rich in a single taxon cannot make one slow chunk.

Usage : confirm_blast_prep.py [--chunk 2500 --seed 11]
Output: queries.tsv          qid sample read_id chunk
        chunk_NNNN.fa        headers >Qxxxxxxxx
"""
from __future__ import annotations

import argparse
import glob
import random
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chunk", type=int, default=2500)
    p.add_argument("--seed", default="11")
    a = p.parse_args()

    samples = sorted(f[:-len(".draw.tsv")] for f in glob.glob("*.draw.tsv"))
    if not samples:
        sys.exit("no <sample>.draw.tsv in the working directory")
    reads = []
    for s in samples:
        seq, name = {}, None
        with open(f"{s}.draw.fa") as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    seq[name] = []
                else:
                    seq[name].append(line.strip())
        with open(f"{s}.draw.tsv") as fh:
            next(fh)
            ids = sorted({line.rstrip("\n").split("\t")[3] for line in fh})
        miss = [r for r in ids if r not in seq]
        if miss:
            sys.exit(f"{s}: {len(miss)} drawn reads missing from {s}.draw.fa")
        reads += [(s, r, "".join(seq[r])) for r in ids]
    reads.sort(key=lambda x: (x[0], x[1]))
    random.Random(f"{a.seed}|blast-chunks").shuffle(reads)

    with open("queries.tsv", "w") as fh:
        fh.write("qid\tsample\tread_id\tchunk\n")
        for i, (s, r, _) in enumerate(reads):
            fh.write(f"Q{i + 1:08d}\t{s}\t{r}\tchunk_{i // a.chunk + 1:04d}\n")
    for c in range((len(reads) + a.chunk - 1) // a.chunk):
        with open(f"chunk_{c + 1:04d}.fa", "w") as fh:
            for i in range(c * a.chunk, min(len(reads), (c + 1) * a.chunk)):
                fh.write(f">Q{i + 1:08d}\n{reads[i][2]}\n")
    print(f"{len(samples)} samples, {len(reads)} queries, chunks of {a.chunk}")


if __name__ == "__main__":
    main()
