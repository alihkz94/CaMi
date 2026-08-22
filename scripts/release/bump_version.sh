#!/usr/bin/env bash
# =============================================================================
# bump_version.sh — move CaMi to a new version, in every place at once
# -----------------------------------------------------------------------------
# USE:
#   bash scripts/release/bump_version.sh 1.1.0
#   bash scripts/release/bump_version.sh 1.1.0 --containers   # also bump the images
#
# WHAT IT TOUCHES
#   nextflow.config   manifest.version
#   CHANGELOG.md      inserts a dated skeleton heading for the new version
#   conf/containers.config   container_tag, only with --containers
#
# WHAT IT DOES NOT DO
#   It does not commit, tag or push. Those are your decision and they come after
#   you have written the changelog entry, which is the only part of a release
#   that a script cannot do for you.
#
# WHEN TO PASS --containers
#   Only when something in containers/ or env/ actually changed. The image tag is
#   deliberately independent of the pipeline version: a documentation fix should
#   not force every cluster running CaMi to re-pull several gigabytes. See
#   conf/containers.config.
#
# VERSIONING SCHEME — semantic, and about REPRODUCIBILITY rather than API:
#   MAJOR  results can change. A threshold moved, a filter rule changed, a tool
#          version changed. Anyone comparing across this boundary must re-run.
#   MINOR  new capability, same numbers. A step, a profile, a parameter.
#   PATCH  documentation, resources, error messages. Numbers cannot move.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

NEW="${1:-}"
BUMP_CONTAINERS="${2:-}"

if [[ -z "$NEW" ]]; then
    echo "USE: $0 <new-version> [--containers]" >&2
    echo "     e.g. $0 1.1.0" >&2
    exit 2
fi
if [[ ! "$NEW" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "ERROR: '$NEW' is not MAJOR.MINOR.PATCH" >&2
    exit 2
fi

OLD=$(sed -n "s/^[[:space:]]*version[[:space:]]*=[[:space:]]*'\([^']*\)'.*/\1/p" \
      "$REPO/nextflow.config" | head -1)
[[ -n "$OLD" ]] || { echo "ERROR: could not read manifest.version" >&2; exit 1; }

if [[ "$OLD" == "$NEW" ]]; then
    echo "Already at $NEW. Nothing to do."
    exit 0
fi

# Refuse to go backwards by accident. `sort -V` puts the lower version first, so
# if NEW sorts first it is older than what is already here.
if [[ "$(printf '%s\n%s\n' "$OLD" "$NEW" | sort -V | head -1)" == "$NEW" ]]; then
    echo "ERROR: $NEW is older than the current $OLD." >&2
    echo "  If that is deliberate, edit nextflow.config by hand." >&2
    exit 1
fi

echo "  $OLD  ->  $NEW"

# --- nextflow.config ---------------------------------------------------------
# Anchored to the manifest block's version line, not a bare search-and-replace:
# the string '1.0.0' could appear in a comment or a container tag, and rewriting
# those would be silent damage.
python3 - "$REPO/nextflow.config" "$OLD" "$NEW" <<'PY'
import re, sys
path, old, new = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(path).read()
pat = re.compile(r"(manifest\s*\{[^}]*?\bversion\s*=\s*')" + re.escape(old) + r"(')", re.S)
s2, n = pat.subn(r"\g<1>" + new + r"\g<2>", s)
if n != 1:
    sys.exit(f"ERROR: expected exactly one manifest.version to replace, found {n}")
open(path, "w").write(s2)
print("  nextflow.config       manifest.version updated")
PY

# --- CHANGELOG.md ------------------------------------------------------------
TODAY=$(date +%Y-%m-%d)
python3 - "$REPO/CHANGELOG.md" "$NEW" "$TODAY" <<'PY'
import sys
path, new, today = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(path).read()
anchor = "---\n\n## "
i = s.find(anchor)
if i < 0:
    sys.exit("ERROR: could not find the first '## ' heading in CHANGELOG.md")
entry = (
    "---\n\n"
    f"## {new} — {today} — TITLE ME\n\n"
    "WRITE THIS BEFORE TAGGING. A release whose changelog says 'TITLE ME' is a\n"
    "release nobody can tell apart from the one before it.\n\n"
    "### What changed\n\n"
    "-\n\n"
    "### Does this move any number?\n\n"
    "Say so explicitly, either way. Anyone comparing results across this version\n"
    "boundary needs the answer, and it is the one question a changelog must answer.\n\n"
)
open(path, "w").write(s[:i] + entry + s[i:])
print(f"  CHANGELOG.md          skeleton for {new} inserted — FILL IT IN")
PY

# --- conf/containers.config --------------------------------------------------
if [[ "$BUMP_CONTAINERS" == "--containers" ]]; then
    OLD_TAG=$(sed -n "s/^[[:space:]]*container_tag[[:space:]]*=[[:space:]]*'\([^']*\)'.*/\1/p" \
              "$REPO/conf/containers.config" | head -1)
    sed -i "s/^\([[:space:]]*container_tag[[:space:]]*=[[:space:]]*\)'${OLD_TAG}'/\1'${NEW}'/" \
        "$REPO/conf/containers.config"
    echo "  conf/containers.config  container_tag $OLD_TAG -> $NEW"
else
    echo "  conf/containers.config  UNCHANGED (pass --containers if an image changed)"
fi

echo
bash "$HERE/check_version.sh"
echo
cat <<EOF
NEXT, in this order:
  1. Write the CHANGELOG.md entry. Nothing else in a release matters as much.
  2. git add -A && git commit -m "CaMi $NEW — <what changed>"
  3. git tag -a v$NEW -m "CaMi $NEW"
  4. git push origin main && git push origin v$NEW

The tag is what triggers the container build and the GitHub Release.
EOF
