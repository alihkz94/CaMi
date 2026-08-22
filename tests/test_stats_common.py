#!/usr/bin/env python3
"""
=============================================================================
test_stats_common.py — brute-force ground truth for bin/_stats_common.py
-----------------------------------------------------------------------------
WHY THESE TESTS EXIST
    The first paired permutation test returned 0.0909 where
    the true exact value is 0.0625, because it applied the Phipson & Smyth +1
    correction to an EXHAUSTIVE enumeration. Nothing about the output looked
    wrong. Only a comparison against independently enumerated ground truth
    found it.

    Every permutation p-value here is therefore checked against a second,
    deliberately naive implementation written a different way — itertools and
    Python floats, no chunking, no vectorisation. If the two agree the fast
    path is right; if they drift, the fast path is wrong.

RUN
    envs/stats/bin/python -m unittest discover -s cockle-microbiome/tests -v

    Plain stdlib unittest, on purpose: a test suite must not add a dependency
    to the environment it is testing.
=============================================================================
"""

from __future__ import annotations

import sys
import unittest
from itertools import combinations, product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

import _stats_common as sc  # noqa: E402


# ---------------------------------------------------------------------------
# independent, deliberately naive reference implementations
# ---------------------------------------------------------------------------
def brute_sign_flip_p(diffs: list[float]) -> float:
    """Every sign vector, one at a time, in pure Python. p = b/m."""
    n = len(diffs)
    observed = abs(sum(diffs))
    at_least = 0
    for signs in product((-1.0, 1.0), repeat=n):
        stat = abs(sum(s * d for s, d in zip(signs, diffs)))
        if stat >= observed - 1e-12:
            at_least += 1
    return at_least / (2 ** n)


def brute_two_sample_p(a: list[float], b: list[float]) -> float:
    """Every way of splitting the pooled values, in pure Python. p = b/m."""
    pooled = list(a) + list(b)
    n, n_a = len(pooled), len(a)
    observed = abs(sum(a) / len(a) - sum(b) / len(b))
    at_least = 0
    total = 0
    for pick in combinations(range(n), n_a):
        total += 1
        left = [pooled[i] for i in pick]
        right = [pooled[i] for i in range(n) if i not in pick]
        stat = abs(sum(left) / len(left) - sum(right) / len(right))
        if stat >= observed - 1e-12:
            at_least += 1
    return at_least / total


class TestPairedSignFlip(unittest.TestCase):
    def test_the_exact_bug_from_changelog_lesson_4(self):
        """n=5, every difference the same sign: the case that caught the bug.

        Exhaustive enumeration puts 2 of the 32 sign vectors at or beyond the
        observed statistic, so p = 2/32 = 0.0625. The buggy version applied
        (b+1)/(m+1) = 3/33 = 0.0909.
        """
        diffs = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        p, m, mode = sc.exact_sign_flip_p(diffs)
        self.assertEqual(mode, "exact")
        self.assertEqual(m, 32)
        self.assertAlmostEqual(p, 0.0625, places=12)
        self.assertNotAlmostEqual(p, 3 / 33, places=4)

    def test_matches_brute_force_across_random_cases(self):
        rng = np.random.default_rng(20260821)
        for n in (4, 6, 9, 11):
            for trial in range(4):
                diffs = rng.normal(size=n) * rng.choice([0.5, 5.0])
                fast, m, mode = sc.exact_sign_flip_p(diffs)
                slow = brute_sign_flip_p(list(diffs))
                self.assertEqual(m, 2 ** n)
                self.assertAlmostEqual(
                    fast, slow, places=12,
                    msg=f"n={n} trial={trial}: fast={fast} brute={slow}",
                )

    def test_p_can_never_be_zero_and_respects_its_floor(self):
        """The observed vector is always in the enumeration, so b >= 1."""
        for n in (3, 5, 8):
            diffs = np.arange(1.0, n + 1.0) * 1000.0  # maximally consistent
            p, m, _ = sc.exact_sign_flip_p(diffs)
            self.assertGreater(p, 0.0)
            self.assertAlmostEqual(p, 2.0 / (2 ** n), places=12)

    def test_monte_carlo_uses_the_plus_one_convention(self):
        """A sampled reference set may miss the observed statistic, so +1."""
        rng = np.random.default_rng(7)
        diffs = np.array([9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0])
        p, m, mode = sc.monte_carlo_sign_flip_p(diffs, 1000, rng)
        self.assertEqual(mode, "monte-carlo")
        self.assertEqual(m, 1000)
        self.assertGreaterEqual(p, 1.0 / 1001)

    def test_monte_carlo_approaches_the_exact_value(self):
        rng = np.random.default_rng(99)
        diffs = np.array([1.5, -0.4, 2.2, 0.9, -1.1, 3.0, 0.2, 1.8])
        exact, _, _ = sc.exact_sign_flip_p(diffs)
        approx, _, _ = sc.monte_carlo_sign_flip_p(diffs, 200_000, rng)
        self.assertLess(abs(exact - approx), 0.01)


class TestTwoSamplePermutation(unittest.TestCase):
    def test_matches_brute_force_when_enumerated(self):
        rng = np.random.default_rng(4242)
        for n_a, n_b in ((3, 3), (4, 3), (5, 4), (2, 6)):
            a = rng.normal(size=n_a)
            b = rng.normal(size=n_b) + 0.8
            fast, m, mode = sc.two_sample_permutation_p(
                a, b, n_perm=1000, rng=rng, exact_max_combinations=10_000
            )
            slow = brute_two_sample_p(list(a), list(b))
            self.assertEqual(mode, "exact")
            self.assertAlmostEqual(
                fast, slow, places=12,
                msg=f"{n_a} vs {n_b}: fast={fast} brute={slow}",
            )

    def test_switches_to_monte_carlo_when_enumeration_is_too_large(self):
        """The real cohort sizes must never claim to be exact."""
        rng = np.random.default_rng(1)
        a = rng.normal(size=16)
        b = rng.normal(size=40)
        _, m, mode = sc.two_sample_permutation_p(
            a, b, n_perm=5_000, rng=rng,
            exact_max_combinations=sc.DEFAULT_EXACT_MAX_COMBINATIONS,
        )
        self.assertEqual(mode, "monte-carlo")
        self.assertEqual(m, 5_000)

    def test_recovers_a_planted_difference_and_clears_a_null(self):
        rng = np.random.default_rng(5)
        shifted_p, _, _ = sc.two_sample_permutation_p(
            rng.normal(size=12), rng.normal(size=12) + 6.0,
            n_perm=20_000, rng=rng, exact_max_combinations=0,
        )
        null_p, _, _ = sc.two_sample_permutation_p(
            rng.normal(size=12), rng.normal(size=12),
            n_perm=20_000, rng=rng, exact_max_combinations=0,
        )
        self.assertLess(shifted_p, 0.01)
        self.assertGreater(null_p, 0.05)

    def test_is_seeded_and_reproducible(self):
        a, b = np.arange(6.0), np.arange(6.0) + 1.7
        first, _, _ = sc.two_sample_permutation_p(
            a, b, 5000, np.random.default_rng(3), exact_max_combinations=0)
        again, _, _ = sc.two_sample_permutation_p(
            a, b, 5000, np.random.default_rng(3), exact_max_combinations=0)
        self.assertEqual(first, again)


class TestFilterTestUniverse(unittest.TestCase):
    def _mat(self):
        #                      t0     t1    t2     t3
        # t0 loud everywhere, t1 always below floor, t2 loud but rare,
        # t3 present everywhere but always sub-floor
        return (
            np.array(
                [
                    [1000.0, 1.0, 0.0, 2.0],
                    [900.0, 1.0, 0.0, 2.0],
                    [800.0, 1.0, 500.0, 2.0],
                    [700.0, 0.0, 0.0, 2.0],
                ]
            ),
            ["t0", "t1", "t2", "t3"],
            np.full(4, 1_000_000.0),
        )

    def test_fp_floor_drops_only_taxa_below_it_in_every_sample(self):
        mat, ids, denom = self._mat()
        # floor 0.0101 % of 1e6 pairs = 101 reads
        out, kept, rep = sc.filter_test_universe(
            mat, ids, denom, fp_floor=0.000101, min_prevalence=0.0
        )
        self.assertEqual(kept, ["t0", "t2"])       # t1, t3 never reach the floor
        self.assertEqual(rep.n_taxa_before, 4)
        self.assertEqual(rep.n_taxa_after_fp_floor, 2)

    def test_prevalence_as_a_fraction(self):
        mat, ids, denom = self._mat()
        _, kept, rep = sc.filter_test_universe(
            mat, ids, denom, fp_floor=0.0, min_prevalence=0.5
        )
        # needs >= 2 of 4 samples non-zero: t0 (4), t1 (3), t3 (4). t2 has 1.
        self.assertEqual(kept, ["t0", "t1", "t3"])
        self.assertEqual(rep.min_prevalence_samples, 2)

    def test_prevalence_as_an_absolute_count(self):
        mat, ids, denom = self._mat()
        _, kept, rep = sc.filter_test_universe(
            mat, ids, denom, fp_floor=0.0, min_prevalence=4
        )
        self.assertEqual(rep.min_prevalence_samples, 4)
        self.assertEqual(kept, ["t0", "t3"])       # only these are in all 4

    def test_both_filters_compose_in_the_required_order(self):
        mat, ids, denom = self._mat()
        out, kept, rep = sc.filter_test_universe(
            mat, ids, denom, fp_floor=0.000101, min_prevalence=0.5
        )
        self.assertEqual(kept, ["t0"])             # t2 clears the floor, fails prevalence
        self.assertEqual(out.shape, (4, 1))
        self.assertEqual(rep.n_taxa_after_fp_floor, 2)
        self.assertEqual(rep.n_taxa_after_prevalence, 1)

    def test_boundary_is_inclusive_at_exactly_the_floor(self):
        mat = np.array([[101.0], [0.0]])
        _, kept, _ = sc.filter_test_universe(
            mat, ["edge"], np.full(2, 1_000_000.0),
            fp_floor=0.000101, min_prevalence=0.0,
        )
        self.assertEqual(kept, ["edge"])

    def test_rejects_mismatched_inputs(self):
        mat = np.ones((3, 2))
        with self.assertRaises(ValueError):
            sc.filter_test_universe(mat, ["a"], np.ones(3), 0.0, 0.0)
        with self.assertRaises(ValueError):
            sc.filter_test_universe(mat, ["a", "b"], np.ones(2), 0.0, 0.0)
        with self.assertRaises(ValueError):
            sc.filter_test_universe(mat, ["a", "b"], np.zeros(3), 0.0, 0.0)


class TestClr(unittest.TestCase):
    def test_rows_sum_to_zero(self):
        rng = np.random.default_rng(2)
        out = sc.clr(rng.uniform(1, 100, size=(5, 7)))
        np.testing.assert_allclose(out.sum(axis=1), 0.0, atol=1e-10)

    def test_is_scale_invariant(self):
        """The whole point of CLR: multiplying a sample by k changes nothing."""
        x = np.array([[1.0, 2.0, 4.0, 8.0]])
        np.testing.assert_allclose(sc.clr(x), sc.clr(x * 37.0), atol=1e-12)

    def test_refuses_zeros_instead_of_returning_infinities(self):
        with self.assertRaises(ValueError):
            sc.clr(np.array([[1.0, 0.0]]))


class TestDetectabilityCeiling(unittest.TestCase):
    def test_paired_ceiling_is_the_hard_two_over_two_to_the_n(self):
        c = sc.detectability_ceiling("paired-exact", n_taxa=100, n_pairs=17)
        self.assertAlmostEqual(c.smallest_p, 2.0 / (1 << 17), places=15)
        self.assertTrue(c.reachable)             # ~3200 taxa screenable at n=17
        self.assertIn("HARD limit", "\n".join(c.lines))

    def test_paired_ceiling_warns_when_too_many_taxa_are_tested(self):
        c = sc.detectability_ceiling("paired-exact", n_taxa=5000, n_pairs=9)
        self.assertAlmostEqual(c.smallest_p, 2.0 / 512, places=15)
        self.assertFalse(c.reachable)            # ~12 taxa at n=9
        self.assertIn("ABOVE that limit", "\n".join(c.lines))

    def test_mc_floor_is_a_budget_and_says_so(self):
        c = sc.detectability_ceiling(
            "mc-resolution", n_taxa=500, n_perm=1_000_000,
            n_group_a=16, n_group_b=405,
        )
        self.assertAlmostEqual(c.smallest_p, 1.0 / 1_000_001, places=15)
        text = "\n".join(c.lines)
        self.assertIn("budget, not", text)
        self.assertIn("POWER", text)
        self.assertNotIn("HARD limit", text)

    def test_raising_the_budget_lowers_the_mc_floor(self):
        low = sc.detectability_ceiling("mc-resolution", n_taxa=10, n_perm=1_000)
        high = sc.detectability_ceiling("mc-resolution", n_taxa=10, n_perm=1_000_000)
        self.assertLess(high.smallest_p, low.smallest_p)

    def test_each_kind_requires_its_own_inputs(self):
        with self.assertRaises(ValueError):
            sc.detectability_ceiling("paired-exact", n_taxa=10)
        with self.assertRaises(ValueError):
            sc.detectability_ceiling("mc-resolution", n_taxa=10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
