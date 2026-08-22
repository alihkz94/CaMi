#!/usr/bin/env python3
# =============================================================================
# qc_checkpoints.py
# -----------------------------------------------------------------------------
# PURPOSE:
#   Independent quality-control gate for the per-sample microbiome profiles
#   produced by profile_microbiome.sh. It re-derives read counts directly from
#   the FASTQ files (it does NOT trust the pipeline's own summary), then applies
#   explicit PASS / WARN / FAIL checks so that unreliable samples cannot be
#   mistaken for real results. This is the "monitoring" layer: nothing is
#   reported as a community profile unless it clears these gates.
#
# CHECKS (per sample):
#   1. stage files present (sub, trim, nohost, nonhost_predup, nonhost)   [FAIL if missing]
#   2. mate pairing: R1 count == R2 count at every stage                  [FAIL if unequal]
#   3. read conservation: counts never increase downstream               [FAIL if violated]
#   4. host removal fraction within a plausible range (default 80–99.99%) [WARN if outside]
#   5. microbial read pairs >= threshold to profile reliably             [FAIL < min_fail,
#                                                                          WARN < min_ok]
#   6. dominant classified genus is a plausible marine/host-associated    [WARN if not]
#      taxon rather than a human/reagent contaminant
#   7. Kraken2 report and Bracken species table exist and are non-empty   [WARN if empty]
#
# CONDA ENV: cockleseq  (seqkit for fast counts; pure-python otherwise)
# RUN:
#   conda activate cockleseq
#   python phase3_fullscale/qc_checkpoints.py                 # all samples in profiles/
#   python phase3_fullscale/qc_checkpoints.py ERR10673701     # specific sample(s)
# OUTPUT:
#   phase3_fullscale/profiles/qc_report.tsv  + printed summary
#   Exit code 0 if no sample FAILs, 1 otherwise (so Nextflow/CI can gate on it).
#
# VERSION: 1.0 (2026-07-13)
# =============================================================================
import os, sys, subprocess, glob

HERE = os.path.dirname(os.path.abspath(__file__))
# Directory holding per-sample subfolders. Override with env var PROFILES_DIR
# (used by the Nextflow pipeline, which stages results elsewhere).
PROFILES = os.environ.get("PROFILES_DIR", os.path.join(HERE, "profiles"))
STAGES = ["sub", "trim", "nohost", "nonhost_predup", "nonhost"]
# subdirectories under profiles/ that are pipeline OUTPUTS, not samples — never QC them
NON_SAMPLE_DIRS = {"aggregated"}

# thresholds (tunable)
MIN_MICROBIAL_OK   = 100   # >= this many microbial pairs -> reliable enough to profile
MIN_MICROBIAL_FAIL = 10    # < this -> at the noise floor, not a community
HOST_PCT_LOW, HOST_PCT_HIGH = 80.0, 99.99

# taxa sanity: genera we EXPECT for a marine bivalve; and ones that signal contamination
MARINE_GENERA = {"Vibrio","Aliivibrio","Photobacterium","Halomonas","Shewanella",
    "Polaribacter","Pseudoalteromonas","Colwellia","Marinomonas","Psychrobacter",
    "Tenacibaculum","Arcobacter","Mycoplasma","Spirochaeta","Roseobacter","Sulfitobacter",
    "Allomuricauda","Muricauda","Maribacter","Zunongwangia","Winogradskyella","Aquimarina"}
CONTAMINANT_HINTS = {"Homo","Ralstonia","Bradyrhizobium","Cutibacterium","Delftia",
    "Pseudomonas","Escherichia","Sphingomonas","Acinetobacter"}


def count_reads(fq):
    """Fast read count via seqkit; returns int or None if file missing/empty."""
    if not os.path.exists(fq):
        return None
    try:
        out = subprocess.run(["seqkit","stats","-T","-j","8",fq],
                             capture_output=True, text=True, check=True).stdout
        return int(out.strip().splitlines()[1].split("\t")[3])
    except Exception:
        return None


def dominant_genus(bracken_species_tsv):
    """Return (genus, reads) of the most abundant species' genus, or (None,0)."""
    if not os.path.exists(bracken_species_tsv) or os.path.getsize(bracken_species_tsv) == 0:
        return (None, 0)
    best_name, best_reads = None, -1
    with open(bracken_species_tsv) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                reads = int(f[5])
            except ValueError:
                continue
            if reads > best_reads:
                best_reads, best_name = reads, f[0]
    genus = best_name.split()[0] if best_name else None
    return (genus, max(best_reads, 0))


def microbial_pairs(kraken_report):
    """Sum Bacteria+Archaea+Viruses (rank D) + Fungi from a kraken2 report."""
    if not os.path.exists(kraken_report):
        return None
    tot = 0
    with open(kraken_report) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            name = f[5].strip(); rank = f[3]
            if (rank == "D" and name in ("Bacteria","Archaea","Viruses")) or name == "Fungi":
                tot += int(f[1])  # col2 = reads in whole CLADE (not col3 = direct-only)
    return tot


def qc_sample(acc):
    d = os.path.join(PROFILES, acc)
    checks = []                       # (name, level, detail)  level in PASS/WARN/FAIL
    counts = {}
    # 1+2: stage presence and pairing. A stage absent from BOTH mates is skipped
    # (intermediate files are not always published, e.g. under Nextflow); a stage
    # present in only ONE mate is a real inconsistency and FAILs.
    for st in STAGES:
        r1 = count_reads(os.path.join(d, f"{st}_1.fastq.gz"))
        r2 = count_reads(os.path.join(d, f"{st}_2.fastq.gz"))
        counts[st] = r1
        if r1 is None and r2 is None:
            continue                      # stage not present here — skip
        elif r1 is None or r2 is None:
            checks.append((f"pairing:{st}", "FAIL", "only one mate file present"))
        elif r1 != r2:
            checks.append((f"pairing:{st}", "FAIL", f"R1={r1} != R2={r2}"))
        else:
            checks.append((f"pairing:{st}", "PASS", f"{r1} pairs"))
    # the final non-host set MUST exist
    if counts.get("nonhost") is None:
        checks.append(("final_output", "FAIL", "nonhost FASTQ missing"))
    # 3: conservation (counts must be non-increasing along STAGES that exist)
    seq = [counts[s] for s in STAGES if counts.get(s) is not None]
    if any(seq[i+1] > seq[i] for i in range(len(seq)-1)):
        checks.append(("conservation", "FAIL", " -> ".join(map(str, seq))))
    else:
        checks.append(("conservation", "PASS", " -> ".join(map(str, seq))))
    # 4: host removal fraction (sub -> nohost)
    if counts.get("trim") and counts.get("nohost") is not None and counts["trim"] > 0:
        host_pct = 100.0 * (counts["trim"] - counts["nohost"]) / counts["trim"]
        lvl = "PASS" if HOST_PCT_LOW <= host_pct <= HOST_PCT_HIGH else "WARN"
        checks.append(("host_removal_pct", lvl, f"{host_pct:.2f}%"))
    # 5: microbial sufficiency
    mp = microbial_pairs(os.path.join(d, f"{acc}.kraken2.report"))
    if mp is None:
        checks.append(("microbial_count", "FAIL", "no kraken2 report"))
    elif mp < MIN_MICROBIAL_FAIL:
        checks.append(("microbial_count", "FAIL", f"{mp} pairs < {MIN_MICROBIAL_FAIL} (noise floor)"))
    elif mp < MIN_MICROBIAL_OK:
        checks.append(("microbial_count", "WARN", f"{mp} pairs < {MIN_MICROBIAL_OK} (low power)"))
    else:
        checks.append(("microbial_count", "PASS", f"{mp} pairs"))
    # 6: dominant taxon sanity
    genus, greads = dominant_genus(os.path.join(d, f"{acc}.bracken.species.tsv"))
    if genus is None:
        checks.append(("dominant_taxon", "WARN", "no bracken species output"))
    elif genus in MARINE_GENERA:
        checks.append(("dominant_taxon", "PASS", f"{genus} ({greads}) — marine/expected"))
    elif genus in CONTAMINANT_HINTS:
        checks.append(("dominant_taxon", "WARN", f"{genus} ({greads}) — possible contaminant"))
    else:
        checks.append(("dominant_taxon", "WARN", f"{genus} ({greads}) — unexpected, inspect"))
    # overall verdict
    levels = [c[1] for c in checks]
    verdict = "FAIL" if "FAIL" in levels else ("WARN" if "WARN" in levels else "PASS")
    # metrics for host_fraction.csv (re-derived, independent of the pipeline logs)
    work = counts.get("trim"); nonhost = counts.get("nonhost")
    metrics = {
        "run": acc,
        "processed_pairs": work,
        "nonhost_pairs": nonhost,
        "microbial_reads": mp if mp is not None else "",
        "pct_microbial_of_processed": (f"{100.0*mp/work:.5f}" if (mp and work) else "0"),
        "dominant_taxon": genus or "",
        "verdict": verdict,
    }
    return verdict, checks, metrics


def main():
    accs = sys.argv[1:] or sorted(
        os.path.basename(p) for p in glob.glob(os.path.join(PROFILES, "*"))
        if os.path.isdir(p) and os.path.basename(p) not in NON_SAMPLE_DIRS)
    if not accs:
        print("No sample profile directories found in", PROFILES); sys.exit(1)
    rows = [("sample","verdict","check","level","detail")]
    metric_rows = []
    any_fail = False
    print(f"\n{'SAMPLE':<14} {'VERDICT':<8} details")
    print("-"*70)
    for acc in accs:
        verdict, checks, metrics = qc_sample(acc)
        any_fail |= (verdict == "FAIL")
        flag = {"PASS":"✓","WARN":"!","FAIL":"✗"}[verdict]
        print(f"{acc:<14} {flag} {verdict:<6}")
        for name, level, detail in checks:
            if level != "PASS":
                print(f"    [{level}] {name}: {detail}")
            rows.append((acc, verdict, name, level, detail))
        metric_rows.append(metrics)
    out = os.path.join(PROFILES, "qc_report.tsv")
    with open(out, "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    # also emit an independent host_fraction.csv (re-derived counts)
    cols = ["run","processed_pairs","nonhost_pairs","microbial_reads",
            "pct_microbial_of_processed","dominant_taxon","verdict"]
    with open(os.path.join(PROFILES, "microbial_summary.csv"), "w") as fh:
        fh.write(",".join(cols) + "\n")
        for m in metric_rows:
            fh.write(",".join(str(m.get(c,"")) for c in cols) + "\n")
    print("-"*70)
    print(f"QC report written: {out}")
    print("Overall:", "FAIL (>=1 sample unusable)" if any_fail else "all samples usable (PASS/WARN)")
    sys.exit(1 if any_fail else 0)


if __name__ == "__main__":
    main()
