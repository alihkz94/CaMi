#!/usr/bin/env python3
"""confirm_combine.py — confirmation of every detection by three independent methods.

A unit is one (sample, taxon) that Kraken 2 detects under the detection rule.
Kraken 2 is always the first leg. Kaiju (protein level, its own database) and
BLAST (nucleotide alignment against nt) are the other two. All taxids are
resolved through ONE NCBI taxdump (merged.dmp applied) and compared by lineage:
a call agrees with unit T when T is in the call's lineage. docs/methods.md,
section 9, gives the reasoning for every threshold below.

  kaiju  among the reads Kraken 2 placed in T's clade, those Kaiju resolved to the
         unit's rank (or inside T); agree = inside T.
           confirmed      agree >= min-agree and agree / resolved >= min-frac
           disagrees      resolved >= min-agree and the share is < min-frac
           insufficient   otherwise
         species units only: when the species is not confirmed, the same test on
         the species' GENUS can confirm the leg (confirmed_at_genus), because
         protein search rarely resolves species.
  blast  <= draw mate-1 reads; best hits = all hits at the read's top bitscore.
         per read: agree (genus: a best hit inside T at >= genus-pident; species:
         a best hit inside T at >= species-pident over >= min-qcov of the read),
         nohit, uninformative (best hits do not resolve to the rank), far (best
         hit below the identity or coverage needed) or disagree_close.
           confirmed      informative >= min-agree and agree / informative >= min-frac
           no_hit         not confirmed, >= half the reads have no hit at all
           no_close_hit   not confirmed, nohit + uninformative + far >= half the
                          reads: nt holds no close reference (NOT evidence against)
           disagrees      otherwise
  tier   kraken, + kaiju if confirmed, + blast if confirmed -> triple / double / single

Genus units also accept the current genus of species that NCBI has moved out of
the unit since the Kraken 2 database was built (e.g. Moraxella osloensis ->
Faucicola osloensis), when those species hold >= moved-min-share of the unit's reads.

Inputs, in the working directory, for every sample of --samples:
  S.kraken2.report  S.units.tsv  S.k2reads.tsv.gz  S.draw.tsv  S.pairs
  S.kaiju.out.gz    (full Kaiju per-read output: C/U, read, taxid)
and queries.tsv plus every chunk_NNNN.blast.tsv.gz.
Output: confirmation_long.tsv, per_sample/S.confirmation.tsv, blast_reads.tsv.gz,
        summary.tsv, params.json
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import statistics
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import confirm_common as cc  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--samples", required=True, help="file, one sample per line")
p.add_argument("--taxdump", required=True, help="directory with nodes.dmp, names.dmp, merged.dmp")
p.add_argument("--jobs", type=int, default=4)
p.add_argument("--merge-family", default="", help="the family entry the detection rule merged, if any")
p.add_argument("--min-agree", type=int, default=5)
p.add_argument("--min-frac", type=float, default=0.5)
p.add_argument("--genus-pident", type=float, default=90.0)
p.add_argument("--species-pident", type=float, default=97.0)
p.add_argument("--min-qcov", type=float, default=0.8)
p.add_argument("--moved-min-share", type=float, default=0.05)
args = p.parse_args()


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def first_int(path):
    with open(path) as fh:
        return int(fh.read().split()[0])


samples = [s for s in open(args.samples).read().split() if s]
units = defaultdict(list)
problems = []
for s in samples:
    for ext in ("kraken2.report", "units.tsv", "k2reads.tsv.gz", "draw.tsv", "pairs", "kaiju.out.gz"):
        if not os.path.isfile(f"{s}.{ext}"):
            problems.append(f"{s}: no {s}.{ext}")
    if os.path.isfile(f"{s}.units.tsv"):
        with open(f"{s}.units.tsv") as fh:
            next(fh)
            for line in fh:
                _, rank, tid, name, n, ra = line.rstrip("\n").split("\t")
                units[s].append((rank, tid, name, int(n), float(ra)))
qid_of = {}
with open("queries.tsv") as fh:
    next(fh)
    for line in fh:
        q, s, r, _c = line.rstrip("\n").split("\t")
        qid_of[(s, r)] = q
chunks = sorted(glob.glob("chunk_*.blast.tsv.gz"))
want_chunks = {f"{c}.blast.tsv.gz" for c in {line.split("\t")[3].strip() for line in list(open("queries.tsv"))[1:]}}
if want_chunks != set(chunks):
    problems.append(f"BLAST: {len(want_chunks - set(chunks))} chunk results missing")
if problems:
    sys.exit("REFUSING: missing pieces:\n  " + "\n  ".join(problems[:50]))

TX = cc.Taxonomy(args.taxdump)
best = {}
hit_rows = defaultdict(list)
queries_in_chunks = set(qid_of.values())
for c in chunks:
    with gzip.open(c, "rt") as fh:
        for line in fh:
            q, _acc, stax, pid, ln, ql, _ev, bits = line.rstrip("\n").split("\t")
            if q not in queries_in_chunks:
                sys.exit(f"REFUSING: {c} holds a hit for {q}, which is not in queries.tsv")
            hit_rows[q].append((float(bits), stax, float(pid), int(ln) / int(ql)))
for q, hs in hit_rows.items():
    top = max(h[0] for h in hs)
    b = []
    for bits, stax, pid, qc in hs:
        if bits == top:
            for t in stax.split(";"):
                rt = TX.resolve(t)
                if rt:
                    b.append((rt, pid, qc))
    tops = [h for h in hs if h[0] == top]
    best[q] = (max(h[2] for h in tops), max(h[3] for h in tops), b)
del hit_rows


def target_rank(rank, tid):
    if args.merge_family and tid == args.merge_family:
        return "family"
    return "genus" if rank == "genus" else "species"


def resolved_at(call, T, R):
    lin = TX.lineage(call)
    return bool(lin) and (T in lin or any(TX.rank.get(a) == R for a in lin))


def label(call, R):
    a = TX.at_rank(call, R)
    return TX.name.get(a, "?") if a else TX.name.get(TX.resolve(call) or "", "?")


def blast_read(q, T, R, species, Ts):
    if q not in best:
        return "nohit", None
    bp, bq, hits = best[q]
    if species:
        agree = any(T in TX.lineage(t) and pi >= args.species_pident and qc >= args.min_qcov for t, pi, qc in hits)
    else:
        agree = any(pi >= args.genus_pident and not Ts.isdisjoint(TX.lineage(t)) for t, pi, _ in hits)
    if agree:
        return "agree", bp
    if not any(resolved_at(t, T, R) for t, _, _ in hits):
        return "uninformative", bp
    close = args.species_pident if species else args.genus_pident
    if bp < close or bq < args.min_qcov:
        return "far", bp
    return "disagree_close", bp


def one_sample(s):
    us = units.get(s, [])
    problems = []
    pairs = first_int(f"{s}.pairs")
    kj, kaiju_lines = {}, 0
    k2 = {}
    with gzip.open(f"{s}.k2reads.tsv.gz", "rt") as fh:
        for line in fh:
            r, t = line.rstrip("\n").split("\t")
            k2[r] = t
    with gzip.open(f"{s}.kaiju.out.gz", "rt") as fh:
        for line in fh:
            kaiju_lines += 1
            f = line.split("\t", 3)
            if f[0] == "C":
                r = cc.read_id(f[1])
                if r in k2:
                    kj[r] = f[2].strip()
    if kaiju_lines != pairs:
        problems.append(f"{s}: Kaiju classified {kaiju_lines} pairs, Kraken 2 {pairs}")
    if not us:
        return s, [], [], problems
    clades = cc.clade_taxids(f"{s}.kraken2.report", [u[1] for u in us])
    by_tax = defaultdict(list)
    for r, t in k2.items():
        by_tax[t].append(r)
    draws = defaultdict(list)
    with open(f"{s}.draw.tsv") as fh:
        next(fh)
        for line in fh:
            _, rank, tid, r = line.rstrip("\n").split("\t")
            draws[(rank, tid)].append(r)
    rows, reads_out = [], []
    genus_clade = {u[1]: clades[u[1]] for u in us if u[0] == "genus"}
    for rank, tid, name, n_rep, ra in us:
        R = target_rank(rank, tid)
        T = TX.resolve(tid)
        species = rank == "species"
        cl = clades[tid]
        members = [r for t in cl for r in by_tax.get(t, ())]
        Ts = {T} if T else set()
        if T and not species:
            moved = Counter()
            for t in cl:
                if by_tax.get(t):
                    g = TX.at_rank(t, R)
                    if g and T not in TX.lineage(t):
                        moved[g] += len(by_tax[t])
            Ts |= {g for g, n in moved.items() if n >= args.moved_min_share * n_rep}
        if len(members) != n_rep:
            problems.append(f"{s} {rank}:{tid}: {len(members)} per-read vs {n_rep} in report")
            continue
        k_cls = k_res = k_agree = k_compat = 0
        k_other = Counter()
        for r in members:
            c = kj.get(r)
            if c is None or c == "0":
                continue
            k_cls += 1
            lin = TX.lineage(c)
            if not lin or T is None:
                continue
            if not Ts.isdisjoint(lin) or c in TX.lineage(T):
                k_compat += 1
            if resolved_at(c, T, R):
                k_res += 1
                if not Ts.isdisjoint(lin):
                    k_agree += 1
                else:
                    k_other[label(c, R)] += 1
        k_frac = k_agree / k_res if k_res else float("nan")
        if k_agree >= args.min_agree and k_frac >= args.min_frac:
            k_status = "confirmed"
        elif k_res >= args.min_agree and k_frac < args.min_frac:
            k_status = "disagrees"
        else:
            k_status = "insufficient"
        k_level = "species" if (species and k_status == "confirmed") else ("genus" if k_status == "confirmed" else "")
        kg_res = kg_agree = 0
        if species and k_status != "confirmed" and T is not None:
            G = TX.at_rank(T, "genus")
            if G:
                for r in members:
                    c = kj.get(r)
                    if c is None or c == "0" or not TX.lineage(c):
                        continue
                    if resolved_at(c, G, "genus"):
                        kg_res += 1
                        if G in TX.lineage(c):
                            kg_agree += 1
                if kg_agree >= args.min_agree and kg_agree / kg_res >= args.min_frac:
                    k_status, k_level = "confirmed_at_genus", "genus"
        cats, pids, b_other = Counter(), [], Counter()
        for r in draws[(rank, tid)]:
            q = qid_of[(s, r)]
            cat, bp = blast_read(q, T, R, species, Ts) if T else ("uninformative", None)
            cats[cat] += 1
            if bp is not None:
                pids.append(bp)
            if cat in ("disagree_close", "far") and q in best:
                labs = {label(t, R) for t, _, _ in best[q][2]}
                b_other["/".join(sorted(labs)[:2])] += 1
            reads_out.append((s, rank, tid, r, q, cat, bp))
        nq = sum(cats.values())
        n_inf = cats["agree"] + cats["far"] + cats["disagree_close"]
        b_frac = cats["agree"] / n_inf if n_inf else float("nan")
        if n_inf >= args.min_agree and b_frac >= args.min_frac:
            b_status = "confirmed"
        elif cats["nohit"] >= 0.5 * nq:
            b_status = "no_hit"
        elif cats["nohit"] + cats["uninformative"] + cats["far"] >= 0.5 * nq:
            b_status = "no_close_hit"
        else:
            b_status = "disagrees"
        legs = (["kraken"] + (["kaiju"] if k_status in ("confirmed", "confirmed_at_genus") else [])
                + (["blast"] if b_status == "confirmed" else []))
        rows.append(dict(
            sample=s, rank=rank, taxid=tid, kraken_name=name, ncbi_taxid=T or "", ncbi_name=TX.name.get(T or "", ""),
            compare_rank=R, also_accepted=";".join(sorted(TX.name.get(x, x) for x in Ts - {T})),
            kraken_reads=n_rep, rel_abundance=ra,
            kaiju_classified=k_cls, kaiju_resolved=k_res, kaiju_agree=k_agree,
            kaiju_agree_frac=k_frac, kaiju_compatible_frac=(k_compat / k_cls if k_cls else float("nan")),
            kaiju_top_other=top2(k_other), kaiju_status=k_status,
            kaiju_level=k_level, kaiju_genus_resolved=kg_res, kaiju_genus_agree=kg_agree,
            blast_queries=nq, blast_with_hit=nq - cats["nohit"], blast_informative=n_inf, blast_agree=cats["agree"],
            blast_agree_frac=b_frac, blast_far=cats["far"], blast_disagree_close=cats["disagree_close"],
            blast_uninformative=cats["uninformative"], blast_nohit=cats["nohit"],
            blast_median_pident=statistics.median(pids) if pids else float("nan"),
            blast_top_other=top2(b_other), blast_status=b_status,
            legs="+".join(legs), tier={3: "triple", 2: "double", 1: "single"}[len(legs)]))
    tier_of = {r["taxid"]: r["tier"] for r in rows if r["rank"] == "genus"}
    for r in rows:
        r["genus_unit_tier"] = ""
        if r["rank"] == "species":
            for g, cl in genus_clade.items():
                if r["taxid"] in cl:
                    r["genus_unit_tier"] = tier_of.get(g, "")
    return s, rows, reads_out, problems


def top2(counter):
    return ", ".join(f"{k} {v}" for k, v in sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))[:2])


def fmt(v):
    if isinstance(v, float):
        return "NA" if v != v else f"{v:.4g}"
    return str(v)


if __name__ == "__main__":
    with Pool(args.jobs) as pool:
        results = pool.map(one_sample, samples, chunksize=4)
    problems = [x for *_, pr in results for x in pr]
    if problems:
        sys.exit(f"REFUSING: {len(problems)} inconsistencies:\n  " + "\n  ".join(problems[:50]))
    os.makedirs("per_sample", exist_ok=True)
    cols = next((rows[0].keys() for _, rows, _, _ in results if rows), None)
    cols = list(cols) if cols else ["sample"]
    long_rows = []
    with gzip.open("blast_reads.tsv.gz", "wt") as br:
        br.write("sample\trank\ttaxid\tread_id\tqid\tcategory\tbest_pident\tbest_hit_taxids\n")
        for s, rows, reads, _ in results:
            for (s_, rank, tid, r, q, cat, bp) in reads:
                tx = ";".join(sorted({t for t, _, _ in best[q][2]})) if q in best else ""
                br.write(f"{s_}\t{rank}\t{tid}\t{r}\t{q}\t{cat}\t{fmt(bp) if bp is not None else 'NA'}\t{tx}\n")
            long_rows += rows
            with open(f"per_sample/{s}.confirmation.tsv", "w") as fh:
                fh.write("\t".join(cols) + "\n")
                for r in rows:
                    fh.write("\t".join(fmt(r[c]) for c in cols) + "\n")
    n_units = sum(len(v) for v in units.values())
    if len(long_rows) != n_units:
        sys.exit(f"REFUSING: {len(long_rows)} rows for {n_units} units")
    with open("confirmation_long.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in long_rows:
            fh.write("\t".join(fmt(r[c]) for c in cols) + "\n")
    summ = Counter((r["rank"], r["tier"], r["legs"], r["kaiju_status"], r["blast_status"]) for r in long_rows)
    with open("summary.tsv", "w") as fh:
        fh.write("rank\ttier\tlegs\tkaiju_status\tblast_status\tunits\n")
        for k, v in sorted(summ.items()):
            fh.write("\t".join(k) + f"\t{v}\n")
    taxmd5 = {f: md5(os.path.join(args.taxdump, f)) for f in ("nodes.dmp", "names.dmp", "merged.dmp")
              if os.path.isfile(os.path.join(args.taxdump, f))}
    params = {k: v for k, v in vars(args).items() if k not in ("samples", "jobs")} | dict(
        taxdump=os.path.realpath(args.taxdump), taxdump_md5=taxmd5,
        samples=len(samples), units=n_units, blast_queries=len(qid_of), blast_chunks=len(chunks),
        tiers={f"{a}:{b}": v for (a, b), v in sorted(Counter((r["rank"], r["tier"]) for r in long_rows).items())},
        outputs_md5={f: md5(f) for f in ("confirmation_long.tsv", "blast_reads.tsv.gz", "summary.tsv")},
        script_md5=md5(os.path.abspath(__file__)))
    with open("params.json", "w") as fh:
        json.dump(params, fh, indent=2)
    print(open("summary.tsv").read())
