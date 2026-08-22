#!/usr/bin/env bash
# =============================================================================
# check_version.sh — one version number, in three places, that must agree
# -----------------------------------------------------------------------------
# CaMi states its version in three places, and each of them is load-bearing:
#
#   nextflow.config   manifest.version   what a RUN reports into results/
#                                        pipeline_info/ — the provenance record
#   CHANGELOG.md      the top heading    what a READER is told changed
#   the git tag       v<version>         what `nextflow run alihkz94/CaMi -r v1.0.0`
#                                        actually pulls
#
# If they drift, the provenance in a results directory names a version whose
# changelog describes different software. That is not a tidiness problem; it is
# the thing version numbers exist to prevent.
#
# USE:
#   bash scripts/release/check_version.sh              # config vs changelog
#   bash scripts/release/check_version.sh v1.0.0       # also against a tag
#
# Exit 0 when they agree, 1 when they do not.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

# manifest.version from nextflow.config
CONFIG_V=$(sed -n "s/^[[:space:]]*version[[:space:]]*=[[:space:]]*'\([^']*\)'.*/\1/p" \
           "$REPO/nextflow.config" | head -1)

# the first "## <version>" heading in the changelog
CHANGELOG_V=$(sed -n 's/^##[[:space:]]\+\([0-9][0-9.]*\).*/\1/p' \
              "$REPO/CHANGELOG.md" | head -1)

fail=0
say() { printf '  %-22s %s\n' "$1" "$2"; }

say "nextflow.config"  "${CONFIG_V:-<none found>}"
say "CHANGELOG.md"     "${CHANGELOG_V:-<none found>}"

[[ -n "$CONFIG_V" ]]    || { echo "ERROR: no manifest.version in nextflow.config" >&2; fail=1; }
[[ -n "$CHANGELOG_V" ]] || { echo "ERROR: no '## <version>' heading in CHANGELOG.md" >&2; fail=1; }

if [[ -n "$CONFIG_V" && -n "$CHANGELOG_V" && "$CONFIG_V" != "$CHANGELOG_V" ]]; then
    echo "ERROR: nextflow.config says $CONFIG_V, CHANGELOG.md says $CHANGELOG_V" >&2
    fail=1
fi

# A tag is only supplied when this runs on a tag push.
if [[ -n "${1:-}" ]]; then
    TAG="${1#refs/tags/}"
    say "git tag" "$TAG"
    if [[ "$TAG" != "v$CONFIG_V" ]]; then
        echo "ERROR: tag is $TAG but manifest.version is $CONFIG_V (expected v$CONFIG_V)" >&2
        fail=1
    fi
fi

if [[ $fail -eq 0 ]]; then
    echo "OK — version $CONFIG_V is consistent."
else
    echo >&2
    echo "Bump all of them together:  bash scripts/release/bump_version.sh <new-version>" >&2
fi
exit $fail
