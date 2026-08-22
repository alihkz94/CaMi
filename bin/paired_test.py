#!/usr/bin/env python3
"""
=============================================================================
paired_test.py — within-individual tumour vs matched-host differential test
-----------------------------------------------------------------------------
PURPOSE
    Test, per taxon, whether abundance differs between a tumour sample and the
    matched-host sample FROM THE SAME COCKLE. Step 2 of the downstream
    statistics.

    Pairing controls for individual genotype, sampling site and batch at once.
    It does NOT control for tissue. Where the two arms are different tissues,
    any difference found is a tissue difference as much as a disease one. The
    pair report prints both tissues so this stays visible; say so in any
    write-up.

METHOD, AND WHY EACH CHOICE IS FORCED BY THE DATA
    1. COMPOSITIONAL, NOT PROPORTIONAL. Sequencing gives relative abundance
       only. Raw proportions produce spurious "changes" whenever an unrelated
       taxon moves. Counts are therefore CLR-transformed (centred log-ratio):
       clr(x)_i = log(x_i) - mean_j(log(x_j)). Differences of CLR values are
       log-ratios, which are the quantity a compositional design can support.
       Required for any compositional design.

    2. PREVALENCE FILTER BEFORE CLR. Sparsity here is extreme — many samples
       carry under 1,000 microbial pairs. CLR is undefined at zero, and a taxon
       seen in 2 of 17 samples cannot support a paired test whatever is done to
       it. Taxa below --min-prevalence are removed BEFORE the transform, so the
       geometric mean is taken over taxa that are actually measured.

    3. PSEUDOCOUNT, APPLIED AFTER FILTERING. Remaining zeros get --pseudocount
       (default 0.5). Applied after filtering so it perturbs a small, dense
       matrix rather than inventing structure across thousands of empty cells.

    4. EXACT PAIRED PERMUTATION TEST. Under the null "tumour and host are
       exchangeable within a cockle", the sign of each pair's CLR difference is
       a coin flip. With n pairs there are exactly 2^n sign assignments, and
       for the n here (about 9-17) they can be ENUMERATED — no asymptotic
       approximation, no normality assumption, exact for tiny n. This is the
       right test at this sample size; a t-test is not.
       Above --exact-max-n the enumeration is replaced by a seeded Monte Carlo
       sample, and the report says which was used.

    5. p-VALUES NEVER REPORTED AS ZERO, AND THE CONVENTION MATCHES THE MODE.
       Exhaustive enumeration uses p = b/m: the observed sign vector is always
       among the m enumerated, so b >= 1 and p >= 1/2^n > 0 automatically.
       Monte Carlo uses the Phipson & Smyth (2010) correction p = (b+1)/(m+1),
       because a sampled reference set is NOT guaranteed to contain the observed
       statistic. Applying the +1 to the exact case would be conservative and
       simply wrong; the reported p_floor states the limit actually in force.

    6. FDR ACROSS TAXA. Benjamini-Hochberg (statsmodels). Both raw and adjusted
       p-values are written; interpret q, not p.

    Wilcoxon signed-rank is reported alongside as a rank-based cross-check. It
    is secondary: at n<10 its own p-value is granular and it is less powerful
    than the exact permutation test on the same data.

INPUTS
    <tables>/counts_long.parquet     from aggregate_counts.py
    <pairs>/pairs.tsv                from build_pairs.py

OUTPUTS (under --outdir, default <pairs>/../11_paired_test)
    paired_results_<method>_<rank>.tsv   per-taxon statistics, sorted by q
    paired_clr_matrix_<method>_<rank>.tsv CLR values actually tested
    paired_test_report.txt                design, filters, and what was found

ENV     envs/stats (polars, numpy, scipy, statsmodels)
RUN     via ./run_pipeline.sh --step STATS
VERSION 1.0 (2026-08-20)
=============================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

# The transform and both permutation modes moved to _stats_common so the
# unpaired host-vs-healthy test and the WGA-vs-native comparison share ONE
# implementation. A second copy of an exact permutation test is a defect
# waiting to happen: the p-value convention bug this replaced was
# subtle enough the first time. Brute-force ground truth for all of it lives in
# tests/test_stats_common.py.
#
# bin/ is on PATH inside a Nextflow task and a script's own directory is
# sys.path[0], so this import resolves the same way in both places.
from _stats_common import clr, exact_sign_flip_p, monte_carlo_sign_flip_p


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Exact paired tumour-vs-matched-host test on CLR abundances.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--tables", type=Path, required=True)
    ap.add_argument("--pairs", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument("--method", default="bracken", choices=["bracken", "kaiju"])
    ap.add_argument("--rank", default="genus",
                    help="genus/species for bracken; species/phylum for kaiju")
    ap.add_argument("--max-tier", type=int, default=2,
                    help="highest pair tier to include (4 = prep-mismatched, not advised)")
    ap.add_argument("--min-prevalence", type=float, default=0.5,
                    help="taxon must be non-zero in this fraction of tested samples")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    ap.add_argument("--exact-max-n", type=int, default=20,
                    help="enumerate all 2^n sign vectors up to this many pairs")
    ap.add_argument("--n-perm", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--include-host-artifact", action="store_true",
                    help="keep Homo (a database artifact here, normally excluded)")
    args = ap.parse_args()

    tables = args.tables.resolve()
    pairs_dir = args.pairs.resolve()
    outdir = (args.outdir or pairs_dir.parent / "11_paired_test").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    counts_path = tables / "counts_long.parquet"
    pairs_path = pairs_dir / "pairs.tsv"
    for p in (counts_path, pairs_path):
        if not p.exists():
            sys.exit(f"ERROR: {p} not found. Run the earlier STATS steps first.")

    tag = f"{args.method}_{args.rank}"
    print("=" * 76)
    print(f" counts : {counts_path}")
    print(f" pairs  : {pairs_path}")
    print(f" design : {args.method} / {args.rank}, tier <= {args.max_tier}")
    print("=" * 76)

    pairs = pl.read_csv(pairs_path, separator="\t")
    if pairs.height == 0:
        _write_empty(outdir, tag, "No pairs available yet.")
        print("  no pairs yet — nothing to test")
        return 0
    pairs = pairs.filter(pl.col("tier") <= args.max_tier)
    if pairs.height == 0:
        _write_empty(outdir, tag, f"No pairs at tier <= {args.max_tier}.")
        print(f"  no pairs at tier <= {args.max_tier}")
        return 0

    counts = pl.read_parquet(counts_path).filter(
        (pl.col("method") == args.method) & (pl.col("rank") == args.rank)
    )
    if not args.include_host_artifact:
        counts = counts.filter(~pl.col("is_host_artifact"))
    if counts.height == 0:
        _write_empty(outdir, tag, f"No counts for {args.method}/{args.rank}.")
        print("  no counts for this method/rank")
        return 0

    # Keep only pairs where BOTH members were actually profiled.
    have = set(counts["sample_code"].unique().to_list())
    pairs = pairs.filter(
        pl.col("tumour_code").is_in(list(have)) & pl.col("host_code").is_in(list(have))
    )
    if pairs.height == 0:
        _write_empty(outdir, tag,
                     "Pairs exist, but no pair has BOTH members profiled yet.")
        print("  no pair has both members profiled yet")
        return 0

    tumour_codes = pairs["tumour_code"].to_list()
    host_codes = pairs["host_code"].to_list()
    used = tumour_codes + host_codes

    wide = (
        counts.filter(pl.col("sample_code").is_in(used))
        .group_by("sample_code", "taxon_id")
        .agg(pl.col("reads").sum())
        .pivot(on="taxon_id", index="sample_code", values="reads")
        .fill_null(0)
    )
    names = (
        counts.select("taxon_id", "taxon_name")
        .unique(subset=["taxon_id"], keep="first")
    )

    # Align rows to the pair order so column i of both arrays is the same cockle.
    order = pl.DataFrame({"sample_code": used})
    wide = order.join(wide, on="sample_code", how="left").fill_null(0)

    taxa_cols = [c for c in wide.columns if c != "sample_code"]
    mat = wide.select(taxa_cols).to_numpy().astype(np.float64)
    n_pairs = pairs.height

    # --- prevalence filter, BEFORE the transform --------------------------
    prevalence = (mat > 0).mean(axis=0)
    keep = prevalence >= args.min_prevalence
    if not keep.any():
        _write_empty(outdir, tag,
                     f"No taxon reaches --min-prevalence {args.min_prevalence} "
                     f"across {mat.shape[0]} samples.")
        print("  no taxon passes the prevalence filter")
        return 0
    mat = mat[:, keep]
    kept_ids = [taxa_cols[i] for i in np.flatnonzero(keep)]
    print(f"  pairs tested            : {n_pairs}")
    print(f"  taxa before prevalence  : {len(taxa_cols)}")
    print(f"  taxa after prevalence   : {len(kept_ids)}")

    # --- CLR ---------------------------------------------------------------
    clr_mat = clr(mat + args.pseudocount)
    diffs_all = clr_mat[:n_pairs, :] - clr_mat[n_pairs:, :]  # tumour - host

    # --- per-taxon test ----------------------------------------------------
    rng = np.random.default_rng(args.seed)
    use_exact = n_pairs <= args.exact_max_n
    rows = []
    for j, tid in enumerate(kept_ids):
        d = diffs_all[:, j]
        if np.allclose(d, 0.0):
            p_perm, m_used, mode = 1.0, 0, "degenerate"
        elif use_exact:
            p_perm, m_used, mode = exact_sign_flip_p(d)
        else:
            p_perm, m_used, mode = monte_carlo_sign_flip_p(d, args.n_perm, rng)

        try:
            p_wil = float(wilcoxon(d, zero_method="wilcox", alternative="two-sided").pvalue)
        except ValueError:
            p_wil = float("nan")  # every difference zero, or n too small

        rows.append(
            {
                "taxon_id": int(tid),
                "mean_clr_diff": float(d.mean()),
                "median_clr_diff": float(np.median(d)),
                "n_pairs": n_pairs,
                "n_tumour_higher": int(np.count_nonzero(d > 0)),
                "p_permutation": float(p_perm),
                "p_wilcoxon": p_wil,
                "permutations": int(m_used),
                "test_mode": mode,
            }
        )

    res = pl.DataFrame(rows)
    p = res["p_permutation"].to_numpy()
    res = res.with_columns(
        pl.Series("q_bh", multipletests(p, method="fdr_bh")[1]),
        pl.Series("p_floor", np.full(p.shape, 1.0 / (1 << n_pairs)
                                     if use_exact else 1.0 / (args.n_perm + 1))),
    )
    res = (
        res.join(names.with_columns(pl.col("taxon_id").cast(pl.Int64)),
                 on="taxon_id", how="left")
        .select("taxon_id", "taxon_name", "n_pairs", "n_tumour_higher",
                "mean_clr_diff", "median_clr_diff", "p_permutation", "q_bh",
                "p_wilcoxon", "test_mode", "permutations", "p_floor")
        .sort("q_bh", "p_permutation")
    )
    res.write_csv(outdir / f"paired_results_{tag}.tsv", separator="\t")

    clr_out = pl.DataFrame(
        {"sample_code": used, "role": ["tumour"] * n_pairs + ["host"] * n_pairs}
    ).hstack(pl.DataFrame(clr_mat, schema=[str(t) for t in kept_ids]))
    clr_out.write_csv(outdir / f"paired_clr_matrix_{tag}.tsv", separator="\t")

    # --- report ------------------------------------------------------------
    sig = res.filter(pl.col("q_bh") < 0.05)
    floor = 1.0 / (1 << n_pairs) if use_exact else 1.0 / (args.n_perm + 1)
    lines = [
        "PAIRED TUMOUR vs MATCHED-HOST TEST",
        "=" * 76,
        f"method / rank        : {args.method} / {args.rank}",
        f"pair tiers included  : <= {args.max_tier}",
        f"pairs tested         : {n_pairs}",
        f"taxa tested          : {len(kept_ids)} "
        f"(of {len(taxa_cols)} before the prevalence filter)",
        f"prevalence threshold : {args.min_prevalence}",
        f"pseudocount          : {args.pseudocount}",
        f"test                 : {'exact enumeration of ' + str(1 << n_pairs) + ' sign vectors' if use_exact else f'Monte Carlo, {args.n_perm} draws, seed {args.seed}'}",
        f"smallest possible p  : {floor:.3g}",
        "",
        f"taxa with q < 0.05   : {sig.height}",
    ]

    # Detectability ceiling. A paired sign-flip test cannot return a p below
    # 2/2^n even when EVERY pair moves the same way, so multiple testing sets a
    # hard limit on how many taxa can be screened before nothing can reach
    # q<0.05. Stating it prevents reading a null result as evidence of absence.
    if use_exact:
        best_p = 2.0 / (1 << n_pairs)
        max_taxa = int(0.05 / best_p) if best_p > 0 else 0
        lines += [
            "",
            "DETECTABILITY CEILING (a property of the design, not of the data):",
            f"  with {n_pairs} pairs the smallest p a perfectly consistent effect can",
            f"  reach is 2/2^{n_pairs} = {best_p:.3g}.",
            f"  After Benjamini-Hochberg across {len(kept_ids)} taxa, q < 0.05 is",
            f"  attainable only while at most ~{max_taxa} taxa are tested.",
        ]
        if max_taxa < len(kept_ids):
            lines += [
                f"  {len(kept_ids)} taxa are being tested, which is ABOVE that limit:",
                "  no taxon can reach q < 0.05 in this configuration regardless of",
                "  effect size. Raise --min-prevalence to test fewer taxa, add pairs",
                "  (--max-tier 2 includes the native non-foot partners), or report",
                "  raw p-values and effect sizes as explicitly hypothesis-generating.",
            ]
    if sig.height:
        lines += ["", "Significant taxa (positive = higher in tumour):"]
        for r in sig.head(25).iter_rows(named=True):
            lines.append(
                f"  {str(r['taxon_name'])[:44]:<44s} "
                f"clr_diff={r['mean_clr_diff']:+.3f}  "
                f"p={r['p_permutation']:.4g}  q={r['q_bh']:.4g}  "
                f"({r['n_tumour_higher']}/{r['n_pairs']} pairs higher)"
            )
    else:
        lines += [
            "",
            "No taxon reaches q < 0.05.",
            "At this sample size that is the expected outcome unless an effect is",
            "large and consistent. Report it as such — it is not evidence of no",
            "difference, it is an underpowered design. The top-ranked taxa in",
            f"paired_results_{tag}.tsv are hypothesis-generating only.",
        ]
    lines += [
        "",
        "STANDING CAVEAT: pairing controls for individual, site and batch. It does",
        "NOT control for tissue — tumour samples are haemolymph and host partners",
        "are foot/mantle/adductor. Any difference found here is a tumour-vs-host",
        "TISSUE difference as much as a disease difference. State this explicitly.",
    ]
    report = "\n".join(lines)
    (outdir / "paired_test_report.txt").write_text(report + "\n")
    print()
    print(report)
    print()
    print(f"WROTE -> {outdir}")
    return 0


def _write_empty(outdir: Path, tag: str, why: str) -> None:
    (outdir / f"paired_results_{tag}.tsv").write_text(
        "taxon_id\ttaxon_name\tn_pairs\tn_tumour_higher\tmean_clr_diff\t"
        "median_clr_diff\tp_permutation\tq_bh\tp_wilcoxon\ttest_mode\t"
        "permutations\tp_floor\n"
    )
    (outdir / "paired_test_report.txt").write_text(
        f"PAIRED TUMOUR vs MATCHED-HOST TEST\n{'=' * 76}\n{why}\n"
        "This is expected while the profiling run is still in progress.\n"
    )


if __name__ == "__main__":
    sys.exit(main())
