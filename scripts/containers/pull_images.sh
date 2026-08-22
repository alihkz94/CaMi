#!/usr/bin/env bash
# =============================================================================
# pull_images.sh — fetch the five CaMi images once, up front
# -----------------------------------------------------------------------------
# TWO REASONS TO RUN THIS RATHER THAN LET NEXTFLOW PULL:
#
#   1. AIR-GAPPED CLUSTERS. Run it on a machine that has a route to ghcr.io,
#      copy the resulting directory across, and start the pipeline with
#      --container_dir. Nothing then needs the network.
#
#   2. A COLD CLUSTER. Without a warm cache, the first few hundred tasks all
#      reach for the same image in the same second. Pulling once first is faster
#      and avoids a half-written cache entry.
#
# USE:
#   bash scripts/containers/pull_images.sh                 # into ./sif
#   bash scripts/containers/pull_images.sh /shared/images  # somewhere shared
#
# THEN:
#   nextflow run . -profile singularity,slurm --container_dir /shared/images ...
#
# ENVIRONMENT:
#   REGISTRY  override the registry  (default: read from conf/containers.config)
#   TAG       override the tag       (default: read from conf/containers.config)
#   ENGINE    singularity | apptainer | docker  (default: whichever is installed)
#
# BEHIND A PROXY (common on an institutional cluster) singularity needs the
# proxy in its own environment, not only in yours:
#   export HTTPS_PROXY=$https_proxy HTTP_PROXY=$http_proxy
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
CONF="$REPO/conf/containers.config"
GROUPS=(fetch qc align classify stats)

read_conf() {
    local key="$1"
    sed -n "s/^[[:space:]]*${key}[[:space:]]*=[[:space:]]*'\\([^']*\\)'.*/\\1/p" "$CONF" | head -1
}
REGISTRY="${REGISTRY:-$(read_conf container_registry)}"
TAG="${TAG:-$(read_conf container_tag)}"
DEST="${1:-$REPO/sif}"

[[ -n "$REGISTRY" && -n "$TAG" ]] || {
    echo "ERROR: could not read container_registry/container_tag from $CONF" >&2
    exit 1
}

if [[ -z "${ENGINE:-}" ]]; then
    for e in singularity apptainer docker; do
        command -v "$e" >/dev/null && { ENGINE="$e"; break; }
    done
fi
[[ -n "${ENGINE:-}" ]] || {
    echo "ERROR: no singularity, apptainer or docker found." >&2
    echo "  With none of them, use -profile conda instead. It is the second-best" >&2
    echo "  option and README.md says why." >&2
    exit 1
}

mkdir -p "$DEST"
echo "=============================================="
echo " engine   : $ENGINE"
echo " registry : $REGISTRY"
echo " tag      : $TAG"
echo " into     : $DEST"
echo "=============================================="

for g in "${GROUPS[@]}"; do
    src="$REGISTRY/cami-$g:$TAG"
    if [[ "$ENGINE" == docker ]]; then
        echo "--- docker pull $src"
        docker pull "$src"
        continue
    fi
    out="$DEST/cami-$g-$TAG.sif"
    # Idempotent: an interrupted copy of a multi-gigabyte image is the thing you
    # least want to re-download, and the thing most likely to be truncated. A
    # finished file is left alone; anything else is replaced.
    if [[ -s "$out" ]] && "$ENGINE" inspect "$out" >/dev/null 2>&1; then
        echo "--- cami-$g: already present and readable, skipping"
        continue
    fi
    echo "--- $ENGINE pull cami-$g"
    rm -f "$out"
    "$ENGINE" pull "$out" "docker://$src"
done

echo
if [[ "$ENGINE" == docker ]]; then
    echo "DONE. Run with:  -profile docker"
else
    echo "DONE. Run with:  -profile $ENGINE --container_dir '$DEST'"
    echo
    echo "To move these to an air-gapped cluster, copy the whole directory:"
    echo "    rsync -a '$DEST/' cluster:/shared/images/"
    ls -lh "$DEST"/cami-*-"$TAG".sif 2>/dev/null || true
fi
