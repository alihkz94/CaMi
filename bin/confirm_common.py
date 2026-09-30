"""Shared code for the confirmation step (confirm_*.py). Standard library only.

  detections(report, rule)  the detection rule applied to one Kraken 2 report
  clade_taxids(report, u)   every taxid inside each unit's clade in that report
  Taxonomy(dir)             one NCBI taxdump (nodes, names, merged) that resolves
                            the taxids of all three legs; lineages are compared
                            by taxid, never by name
"""
from __future__ import annotations

import gzip
import os
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Rule:
    min_reads: int = 10
    min_frac: float = 0.01
    exclude: frozenset = field(default_factory=lambda: frozenset({"9606", "9605"}))
    merge_family: str | None = None


def report_rows(path):
    """(depth, rank_code, taxid, name, clade_reads) per line of a Kraken 2 report."""
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name = f[5]
            yield len(name) - len(name.lstrip(" ")), f[3], f[4], name.strip(), int(f[1])


def detections(report, rule: Rule):
    """Rows (rank, taxid, name, reads, rel_abundance, detected) and the denominator.

    Genus and species rows of the report. The genera inside rule.merge_family are
    replaced by one entry for the family clade, kept under rank 'genus'. The
    denominator is the sum of the genus-rank rows, for both ranks."""
    rows, inside, d0 = [], False, 0
    for depth, rk, tid, name, n in report_rows(report):
        if inside and depth <= d0:
            inside = False
        if rule.merge_family and tid == rule.merge_family:
            inside, d0 = True, depth
            if n > 0:
                rows.append(("genus", tid, name, n))
            continue
        if inside and rk == "G":
            continue
        if rk in ("G", "S") and tid not in rule.exclude and n > 0:
            rows.append(("genus" if rk == "G" else "species", tid, name, n))
    total = sum(r[3] for r in rows if r[0] == "genus")
    out = []
    for rank, tid, name, n in rows:
        ra = n / total if total else 0.0
        out.append((rank, tid, name, n, ra, n >= rule.min_reads and ra >= rule.min_frac))
    return out, total


def clade_taxids(report, unit_taxids):
    """{unit_taxid: set(taxids in its clade, itself included)} from the report tree,
    so the clade read count equals the report's clade count exactly."""
    want = set(unit_taxids)
    out = {t: set() for t in want}
    open_units = []
    for depth, _rk, tid, _name, _n in report_rows(report):
        while open_units and depth <= open_units[-1][1]:
            open_units.pop()
        for u, _ in open_units:
            out[u].add(tid)
        if tid in want:
            out[tid].add(tid)
            open_units.append((tid, depth))
    return out


def read_id(name):
    """A read name as all three tools report it: first word, no /1 or /2 mate suffix."""
    name = name.split()[0] if name else name
    return name[:-2] if name.endswith(("/1", "/2")) else name


class Taxonomy:
    """NCBI taxdump. resolve() follows merged.dmp; deleted or unknown ids -> None."""

    def __init__(self, d):
        self.parent, self.rank, self.name = {}, {}, {}
        with open(os.path.join(d, "nodes.dmp")) as fh:
            for line in fh:
                f = line.split("\t|\t", 3)
                self.parent[f[0]] = f[1]
                self.rank[f[0]] = f[2]
        with open(os.path.join(d, "names.dmp")) as fh:
            for line in fh:
                if "scientific name" in line:
                    f = line.split("\t|\t", 2)
                    self.name[f[0]] = f[1]
        self.merged = {}
        merged = os.path.join(d, "merged.dmp")
        if os.path.isfile(merged):
            with open(merged) as fh:
                for line in fh:
                    f = line.split("\t|\t")
                    self.merged[f[0]] = f[1].split("\t")[0].strip()
        self._lin = {}

    def resolve(self, tid):
        tid = str(tid).strip()
        if tid in self.parent:
            return tid
        seen = 0
        while tid in self.merged and seen < 10:
            tid, seen = self.merged[tid], seen + 1
        return tid if tid in self.parent else None

    def lineage(self, tid):
        """Tuple of taxids from tid up to the root (tid first); () if unresolved."""
        t = self.resolve(tid)
        if t is None:
            return ()
        if t in self._lin:
            return self._lin[t]
        lin, cur = [], t
        while True:
            lin.append(cur)
            p = self.parent.get(cur)
            if p is None or p == cur:
                break
            cur = p
        lin = tuple(lin)
        self._lin[t] = lin
        return lin

    def at_rank(self, tid, rank):
        for a in self.lineage(tid):
            if self.rank.get(a) == rank:
                return a
        return None

    def within(self, tid, target):
        tg = self.resolve(target)
        return tg is not None and tg in self.lineage(tid)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
