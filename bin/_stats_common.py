#!/usr/bin/env python3
"""
=============================================================================
_stats_common.py — shared statistics for every STATS test
-----------------------------------------------------------------------------
PURPOSE
    One home for the transforms, permutation tests, test-universe filter and
    power arithmetic that more than one STATS script needs. Before this module
    existed, `paired_test.py` held them privately. The unpaired host-vs-healthy
    test and the WGA-vs-native descriptive comparison need the SAME code, and a
    second copy of an exact permutation test is a defect waiting to happen —
    the p-value convention bug this replaced was subtle enough the
    first time.

    The leading underscore marks this as a library, not a command. Every other
    file in bin/ is an entry point.

WHAT IS HERE
    clr()                       centred log-ratio, per sample
    filter_test_universe()      FP floor, then prevalence  (Task 6)
    exact_sign_flip_p()         PAIRED, exhaustive over 2^n sign vectors
    monte_carlo_sign_flip_p()   PAIRED, seeded sample above that
    two_sample_permutation_p()  UNPAIRED, group-label permutation
    detectability_ceiling()     what the design can possibly detect

THE ONE THING TO GET RIGHT: p-VALUE CONVENTION FOLLOWS THE MODE
    Exhaustive enumeration:  p = b/m.  The observed vector is always among the
        m enumerated, so b >= 1 and p >= 1/m > 0 by construction.
    Sampled (Monte Carlo):   p = (b+1)/(m+1), the Phipson & Smyth (2010)
        correction, because a SAMPLED reference set is not guaranteed to
        contain the observed statistic.
    Applying the +1 to the exact case is conservative and simply wrong. It was
    a real bug once (an earlier version returned 0.0909 where the true exact
    value is 0.0625) and is now covered by brute-force tests.

THE OTHER THING: TWO DIFFERENT KINDS OF "FLOOR"
    A paired sign-flip test has a HARD floor of 2/2^n set by the sample size.
    No amount of computing time lowers it. That is `kind="paired-exact"`.

    An unpaired test on groups this size cannot be enumerated at all — the
    reference set for 16 vs 405 is C(421,16), about 10^28 — so it is Monte
    Carlo from the start and its floor is 1/(B+1), a CHOICE of how many
    permutations to run. That is `kind="mc-resolution"`. Raising B lowers it.
    Small n limits POWER there, not the p-value floor. Never present the two
    as the same quantity: one is a property of the cohort, the other of the
    budget.

ENV     envs/stats (numpy)
VERSION 1.0 (2026-08-21)
=============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import comb
from typing import Final, Literal, Sequence

import numpy as np

# Sign vectors / permutations materialised at once. Bounds peak memory: at
# n=20 the full paired set is 1,048,576 vectors, handled 65,536 at a time.
CHUNK: Final[int] = 1 << 16

# Above this many distinct group assignments, an unpaired test stops trying to
# enumerate and draws a seeded sample instead. 2^20 keeps exact enumeration to
# roughly the same work as the paired test's own ceiling.
DEFAULT_EXACT_MAX_COMBINATIONS: Final[int] = 1 << 20


# ---------------------------------------------------------------------------
# transform
# ---------------------------------------------------------------------------
def clr(matrix: np.ndarray) -> np.ndarray:
    """
    Centred log-ratio, per row (one row = one sample).

    Input must be STRICTLY POSITIVE — add the pseudocount before calling, and
    add it AFTER the prevalence filter so it perturbs a small dense matrix
    instead of inventing structure across thousands of empty cells.
    """
    if np.any(matrix <= 0):
        raise ValueError(
            "clr() needs strictly positive input; add the pseudocount first"
        )
    log_x = np.log(matrix)
    return log_x - log_x.mean(axis=1, keepdims=True)


# ---------------------------------------------------------------------------
# the test universe  (Task 6)
# ---------------------------------------------------------------------------
@dataclass
class UniverseReport:
    """What the filter did, so the count that feeds the power ceiling is auditable."""

    n_taxa_before: int
    n_taxa_after_fp_floor: int
    n_taxa_after_prevalence: int
    n_samples: int
    fp_floor: float
    min_prevalence: float
    min_prevalence_samples: int
    lines: list[str] = field(default_factory=list)


def filter_test_universe(
    mat: np.ndarray,
    taxon_ids: Sequence,
    denominators: np.ndarray,
    fp_floor: float,
    min_prevalence: float,
) -> tuple[np.ndarray, list, UniverseReport]:
    """
    Cut the taxon table down to what can carry a test, in the required order.

    1. FALSE-POSITIVE FLOOR. 2,000,000 fragments of pure cockle DNA still call
       0.0101 % bacterial through this pipeline (docs/methods.md). A taxon whose
       abundance stays below that floor in EVERY sample is indistinguishable
       from that noise, so it is dropped. Abundance is measured against the
       sample's own non-host pair count, never against a per-rank denominator
       (Bracken does not sum to Kraken 2 — see docs/statistics.md).

    2. PREVALENCE. A taxon seen in a handful of samples cannot support a test
       whatever is done to it, and every one kept costs multiple-testing
       burden that the power ceiling then has to pay for.

    `min_prevalence` reads two ways, on purpose:
        < 1   a FRACTION of samples   (0.5 = half of them; the paired default)
        >= 1  an ABSOLUTE sample count (3 = "present in at least 3 samples";
              the form Task 6 asks for on the unpaired tests)

    Returns (filtered matrix, kept taxon ids, report).
    """
    n_samples, n_taxa = mat.shape
    if len(taxon_ids) != n_taxa:
        raise ValueError(
            f"taxon_ids has {len(taxon_ids)} entries for {n_taxa} matrix columns"
        )
    denominators = np.asarray(denominators, dtype=np.float64)
    if denominators.shape != (n_samples,):
        raise ValueError(
            f"denominators must be one non-host pair count per sample "
            f"({n_samples}), got {denominators.shape}"
        )
    if np.any(denominators <= 0):
        raise ValueError("a sample has a non-positive non-host pair count")

    # --- 1: false-positive floor -----------------------------------------
    fractions = mat / denominators[:, None]
    above_floor = fractions.max(axis=0) >= fp_floor
    n_after_fp = int(np.count_nonzero(above_floor))

    # --- 2: prevalence ----------------------------------------------------
    if min_prevalence >= 1:
        needed = int(min_prevalence)
    else:
        needed = int(np.ceil(min_prevalence * n_samples))
    prevalent = (mat > 0).sum(axis=0) >= needed

    keep = above_floor & prevalent
    kept_idx = np.flatnonzero(keep)
    kept_ids = [taxon_ids[i] for i in kept_idx]

    report = UniverseReport(
        n_taxa_before=n_taxa,
        n_taxa_after_fp_floor=n_after_fp,
        n_taxa_after_prevalence=len(kept_ids),
        n_samples=n_samples,
        fp_floor=fp_floor,
        min_prevalence=min_prevalence,
        min_prevalence_samples=needed,
    )
    report.lines = [
        "TEST UNIVERSE",
        f"  samples                    : {n_samples}",
        f"  taxa before any filter     : {n_taxa}",
        f"  after FP floor ({fp_floor:.3g}) : {n_after_fp}",
        f"  after prevalence (>= {needed} samples) : {len(kept_ids)}",
    ]
    return mat[:, kept_idx], kept_ids, report


# ---------------------------------------------------------------------------
# paired permutation  (sign flip)
# ---------------------------------------------------------------------------
def exact_sign_flip_p(diffs: np.ndarray) -> tuple[float, int, str]:
    """
    Two-sided exact paired p-value, enumerating every one of the 2^n sign
    vectors. Statistic is the sum of differences. Chunked so memory stays
    bounded.

    p = b/m. See the module docstring: the +1 belongs to the sampled case only.
    """
    n = diffs.size
    total = 1 << n
    observed = float(np.abs(diffs.sum()))
    idx = np.arange(n)
    at_least = 0
    for start in range(0, total, CHUNK):
        stop = min(start + CHUNK, total)
        block = np.arange(start, stop, dtype=np.int64)[:, None]
        # bit j of the block index selects the sign of pair j
        signs = np.where((block >> idx) & 1, 1.0, -1.0)
        stats = signs @ diffs
        at_least += int(np.count_nonzero(np.abs(stats) >= observed - 1e-12))
    return at_least / total, total, "exact"


def monte_carlo_sign_flip_p(
    diffs: np.ndarray, n_perm: int, rng: np.random.Generator
) -> tuple[float, int, str]:
    """Seeded paired fallback when 2^n is too large to enumerate."""
    n = diffs.size
    observed = float(np.abs(diffs.sum()))
    at_least = 0
    drawn = 0
    while drawn < n_perm:
        take = min(CHUNK, n_perm - drawn)
        signs = rng.choice((-1.0, 1.0), size=(take, n))
        stats = signs @ diffs
        at_least += int(np.count_nonzero(np.abs(stats) >= observed - 1e-12))
        drawn += take
    return (at_least + 1) / (n_perm + 1), n_perm, "monte-carlo"


# ---------------------------------------------------------------------------
# unpaired permutation  (group label)
# ---------------------------------------------------------------------------
def two_sample_permutation_p(
    values_a: np.ndarray,
    values_b: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
    exact_max_combinations: int = DEFAULT_EXACT_MAX_COMBINATIONS,
) -> tuple[float, int, str]:
    """
    Two-sided unpaired permutation p-value on the difference of group means.

    Under the null the group LABEL carries no information, so every way of
    splitting the pooled values into groups of the observed sizes is equally
    likely. Exact enumeration is used when C(n_a+n_b, n_a) is small enough;
    otherwise a seeded sample is drawn and the p-value convention switches
    accordingly.

    For the cohort this pipeline actually has — 16 vs 405 foot samples, 18 vs 6
    haemolymph — enumeration NEVER triggers: C(421,16) is about 10^28. These
    tests are Monte Carlo by necessity, and their p-value floor is 1/(B+1), a
    budget, not a bound. See detectability_ceiling().
    """
    a = np.asarray(values_a, dtype=np.float64)
    b = np.asarray(values_b, dtype=np.float64)
    n_a, n_b = a.size, b.size
    if n_a == 0 or n_b == 0:
        raise ValueError("both groups must be non-empty")
    pooled = np.concatenate((a, b))
    n = pooled.size
    observed = float(abs(a.mean() - b.mean()))

    total_splits = comb(n, n_a)
    if total_splits <= exact_max_combinations:
        from itertools import combinations

        at_least = 0
        for pick in combinations(range(n), n_a):
            mask = np.zeros(n, dtype=bool)
            mask[list(pick)] = True
            stat = abs(pooled[mask].mean() - pooled[~mask].mean())
            if stat >= observed - 1e-12:
                at_least += 1
        return at_least / total_splits, total_splits, "exact"

    at_least = 0
    drawn = 0
    while drawn < n_perm:
        take = min(CHUNK, n_perm - drawn)
        # argsort of uniform noise gives an unbiased shuffle per row
        order = np.argsort(rng.random((take, n)), axis=1)
        shuffled = pooled[order]
        stat = np.abs(
            shuffled[:, :n_a].mean(axis=1) - shuffled[:, n_a:].mean(axis=1)
        )
        at_least += int(np.count_nonzero(stat >= observed - 1e-12))
        drawn += take
    return (at_least + 1) / (n_perm + 1), n_perm, "monte-carlo"


# ---------------------------------------------------------------------------
# what the design can possibly detect
# ---------------------------------------------------------------------------
@dataclass
class Ceiling:
    kind: str
    smallest_p: float
    n_taxa: int
    max_taxa_reachable: int
    reachable: bool
    lines: list[str] = field(default_factory=list)


def detectability_ceiling(
    kind: Literal["paired-exact", "mc-resolution"],
    n_taxa: int,
    *,
    n_pairs: int | None = None,
    n_perm: int | None = None,
    n_group_a: int | None = None,
    n_group_b: int | None = None,
    alpha: float = 0.05,
) -> Ceiling:
    """
    State the best p-value the design can produce, and whether q < alpha is
    reachable at all after Benjamini-Hochberg across `n_taxa`.

    Two kinds, and conflating them misleads the reader:

    "paired-exact"   HARD limit. A sign-flip test cannot return a p below
                     2/2^n even if every pair moves the same way. This is a
                     property of the cohort. More computing time cannot help;
                     only more pairs can.

    "mc-resolution"  SOFT limit. The floor is 1/(B+1) for B permutations, a
                     choice of budget. Raising B lowers it. What a small group
                     limits here is POWER — the size of effect detectable —
                     not the floor. Reported so nobody reads it as a hard
                     bound of the paired kind.
    """
    if kind == "paired-exact":
        if not n_pairs:
            raise ValueError("paired-exact needs n_pairs")
        smallest = 2.0 / (1 << n_pairs)
        max_taxa = int(alpha / smallest) if smallest > 0 else 0
        lines = [
            "DETECTABILITY CEILING (a property of the design, not of the data):",
            f"  with {n_pairs} pairs the smallest p a perfectly consistent effect",
            f"  can reach is 2/2^{n_pairs} = {smallest:.3g}. This is a HARD limit:",
            "  it is set by the sample size and no amount of computing lowers it.",
            f"  After Benjamini-Hochberg across {n_taxa} taxa, q < {alpha} is",
            f"  attainable only while at most ~{max_taxa} taxa are tested.",
        ]
        if max_taxa < n_taxa:
            lines += [
                f"  {n_taxa} taxa are being tested, which is ABOVE that limit: no",
                "  taxon can reach q < 0.05 in this configuration regardless of",
                "  effect size. Raise the prevalence threshold to test fewer taxa,",
                "  add pairs, or report raw p-values and effect sizes as",
                "  explicitly hypothesis-generating.",
            ]
        return Ceiling(kind, smallest, n_taxa, max_taxa,
                       max_taxa >= n_taxa, lines)

    if not n_perm:
        raise ValueError("mc-resolution needs n_perm")
    smallest = 1.0 / (n_perm + 1)
    max_taxa = int(alpha / smallest) if smallest > 0 else 0
    group_txt = ""
    if n_group_a and n_group_b:
        group_txt = f" (groups of {n_group_a} and {n_group_b})"
    lines = [
        "RESOLUTION FLOOR (a property of the permutation budget, NOT a hard",
        "limit of the design):",
        f"  this test{group_txt} cannot be enumerated exactly, so it draws",
        f"  {n_perm:,} permutations and its smallest reportable p is",
        f"  1/(B+1) = {smallest:.3g}.",
        "  Raising the permutation count lowers this floor. It is a budget, not",
        "  a bound, and it must NOT be read like the 2/2^n ceiling of the",
        "  paired test.",
        f"  After Benjamini-Hochberg across {n_taxa} taxa, q < {alpha} stays",
        f"  reachable while at most ~{max_taxa} taxa are tested.",
        "  What a small group limits here is POWER — only a large effect will",
        "  be detected — not the p-value floor.",
    ]
    if max_taxa < n_taxa:
        lines += [
            f"  {n_taxa} taxa are being tested, above that limit: raise the",
            "  permutation count rather than concluding the effect is absent.",
        ]
    return Ceiling(kind, smallest, n_taxa, max_taxa, max_taxa >= n_taxa, lines)
