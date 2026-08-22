#!/usr/bin/env bash
# =============================================================================
# 12_publish_intermediates.sh
# -----------------------------------------------------------------------------
# PURPOSE:
#   Put the intermediate FASTQ of every finished task into its numbered step
#   folder, so each stage's output for all samples sits together:
#
#     results/01_trimmed/<sample>.trim_[12].fastq.gz
#     results/02_host_removed/<sample>.nohost_[12].fastq.gz
#     results/03_human_removed/<sample>.nonhost_predup_[12].fastq.gz
#
#   These files are large and are NOT published by default (publish_intermediates
#   is false), so they normally stay in work/ where they are hard to find.
#
# HOW — HARD LINKS, NOT COPIES
#   work/ and results/ are on the same filesystem, so a hard link gives the file
#   a second name at ZERO extra disk cost. A copy of the trimmed reads alone
#   would be about 2.5 GB x 2 mates x 498 samples = ~2.5 TB.
#
#   A hard link is not a shortcut that can break: the data lives until the LAST
#   name is removed. So when work/ is cleaned later, these files survive in
#   results/ — which is exactly what we want. The trade is that cleaning work/
#   will free less space than expected, because this data is still referenced.
#
# SAFE TO RUN WHILE THE PIPELINE IS LIVE
#   It only reads work/ and creates links. It never deletes or moves anything.
#   It is idempotent: existing links are skipped, so run it again at any time to
#   pick up newly finished samples.
#
#   It links ONLY from task directories that exited 0, and only real files —
#   never Nextflow's staged input symlinks, which would otherwise be re-linked
#   from a downstream task's directory.
#
# USE:
#   bash 12_publish_intermediates.sh --dry-run   # show what would happen
#   bash 12_publish_intermediates.sh             # do it
#
# VERSION: 1.0 (2026-08-16)
# =============================================================================
set -uo pipefail

PROJECT="/scratch/ahakimzadeh/Alicia_nature_data"
[[ -d "$PROJECT" ]] || PROJECT="/mnt/slurm-406/ahakimzadeh/Alicia_nature_data"
WORK="$PROJECT/work"
OUT="$PROJECT/results"

DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

[[ -d "$WORK" ]] || { echo "ERROR: no work directory at $WORK" >&2; exit 1; }

echo "=============================================="
echo " source : $WORK"
echo " target : $OUT"
[[ $DRY -eq 1 ]] && echo " MODE   : dry run — nothing will be created"
echo "=============================================="

# pattern -> destination step folder
link_stage() {
    local glob="$1" dest="$2" n_new=0 n_have=0 n_skip=0
    mkdir -p "$OUT/$dest" 2>/dev/null

    # -type f excludes symlinks, so staged inputs in a DOWNSTREAM task directory
    # are never mistaken for that task's own output.
    while IFS= read -r src; do
        local d
        d="$(dirname "$src")"
        # only trust output from a task that actually succeeded
        if [[ ! -f "$d/.exitcode" ]] || [[ "$(cat "$d/.exitcode" 2>/dev/null)" != "0" ]]; then
            n_skip=$(( n_skip + 1 )); continue
        fi
        local base tgt
        base="$(basename "$src")"
        tgt="$OUT/$dest/$base"
        if [[ -e "$tgt" ]]; then
            n_have=$(( n_have + 1 )); continue
        fi
        if [[ $DRY -eq 1 ]]; then
            n_new=$(( n_new + 1 ))
        elif ln "$src" "$tgt" 2>/dev/null; then
            n_new=$(( n_new + 1 ))
        elif cp -n "$src" "$tgt" 2>/dev/null; then
            # different filesystem: fall back to a copy rather than fail
            n_new=$(( n_new + 1 ))
        else
            echo "  WARN: could not link or copy $src" >&2
        fi
    done < <(find "$WORK" -type f -name "$glob" 2>/dev/null)

    printf "  %-26s new=%-5d already=%-5d skipped(unfinished)=%d\n" "$dest" "$n_new" "$n_have" "$n_skip"
}

link_stage '*.trim_1.fastq.gz'            01_trimmed
link_stage '*.trim_2.fastq.gz'            01_trimmed
link_stage '*.nohost_1.fastq.gz'          02_host_removed
link_stage '*.nohost_2.fastq.gz'          02_host_removed
link_stage '*.nonhost_predup_1.fastq.gz'  03_human_removed
link_stage '*.nonhost_predup_2.fastq.gz'  03_human_removed

echo
echo "=== files now in each step folder ==="
for d in "$OUT"/0*/; do
    [[ -d "$d" ]] || continue
    printf "  %-26s %5s files  %8s\n" "$(basename "$d")" \
        "$(ls -1 "$d" 2>/dev/null | wc -l)" \
        "$(du -sh --apparent-size "$d" 2>/dev/null | cut -f1)"
done
echo
echo "NOTE: apparent size counts hard-linked data once per name. Real disk use is"
echo "      unchanged — the same blocks are shared with work/."
