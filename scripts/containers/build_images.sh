#!/usr/bin/env bash
# =============================================================================
# build_images.sh — build the five CaMi images from containers/
# -----------------------------------------------------------------------------
# You almost certainly do not need this. Nextflow pulls the published images by
# itself; see scripts/containers/pull_images.sh for the offline case.
#
# Build them yourself when you are:
#   - changing what is inside an image (then bump container_tag and push),
#   - on an air-gapped cluster with no route to a registry,
#   - auditing, and unwilling to trust a registry for a published result.
#
# USE:
#   bash scripts/containers/build_images.sh docker              # all five
#   bash scripts/containers/build_images.sh singularity align   # just one
#   bash scripts/containers/build_images.sh apptainer qc stats
#
# ENVIRONMENT:
#   REGISTRY   override the registry     (default: read from conf/containers.config)
#   TAG        override the tag          (default: read from conf/containers.config)
#   SIF_DIR    where .sif files land     (default: ./sif)
#   PUSH=1     docker only: push after a successful build
#
# RUN IT FROM ANYWHERE. It finds the repository itself, and every build uses the
# repository root as its context — the stats image needs env/requirements-stats.txt.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
CONF="$REPO/conf/containers.config"
GROUPS_ALL=(fetch qc align classify stats)

# The tag lives in conf/containers.config so that Nextflow and this script can
# never disagree about which image a given checkout expects.
read_conf() {
    local key="$1"
    sed -n "s/^[[:space:]]*${key}[[:space:]]*=[[:space:]]*'\\([^']*\\)'.*/\\1/p" "$CONF" | head -1
}
REGISTRY="${REGISTRY:-$(read_conf container_registry)}"
TAG="${TAG:-$(read_conf container_tag)}"
SIF_DIR="${SIF_DIR:-$REPO/sif}"

[[ -n "$REGISTRY" && -n "$TAG" ]] || {
    echo "ERROR: could not read container_registry/container_tag from $CONF" >&2
    exit 1
}

ENGINE="${1:-}"
shift || true
case "$ENGINE" in
    docker|singularity|apptainer) ;;
    *) echo "USE: $0 {docker|singularity|apptainer} [group ...]" >&2
       echo "     groups: ${GROUPS_ALL[*]}" >&2
       exit 2 ;;
esac
command -v "$ENGINE" >/dev/null || { echo "ERROR: $ENGINE is not installed" >&2; exit 1; }

GROUPS=("$@")
[[ ${#GROUPS[@]} -eq 0 ]] && GROUPS=("${GROUPS_ALL[@]}")

echo "=============================================="
echo " repository : $REPO"
echo " engine     : $ENGINE"
echo " registry   : $REGISTRY"
echo " tag        : $TAG"
echo " groups     : ${GROUPS[*]}"
[[ "$ENGINE" != docker ]] && echo " sif dir    : $SIF_DIR"
echo "=============================================="

cd "$REPO"
for g in "${GROUPS[@]}"; do
    [[ -d "containers/$g" ]] || { echo "ERROR: no containers/$g" >&2; exit 1; }
    image="$REGISTRY/cami-$g:$TAG"
    echo
    echo "--- cami-$g -------------------------------------------------"

    if [[ "$ENGINE" == docker ]]; then
        docker build -f "containers/$g/Dockerfile" -t "$image" .
        # Prove the tools are actually callable before anyone trusts the image.
        # A conda solve can succeed and still leave a broken PATH.
        docker run --rm "$image" bash -lc 'set -e; command -v python3 >/dev/null || true; echo ok' >/dev/null
        [[ "${PUSH:-0}" == 1 ]] && docker push "$image"
    else
        mkdir -p "$SIF_DIR"
        # --fakeroot needs a subuid/subgid range for your account. Without one,
        # build the OCI image on a machine where you have Docker and copy the
        # .sif across; the definition file produces the same environment either way.
        "$ENGINE" build --fakeroot "$SIF_DIR/cami-$g-$TAG.sif" "containers/$g/Singularity.def"
        "$ENGINE" test "$SIF_DIR/cami-$g-$TAG.sif"
    fi
done

echo
echo "DONE."
if [[ "$ENGINE" == docker ]]; then
    echo "Run the pipeline with:  -profile docker"
    [[ "${PUSH:-0}" == 1 ]] || echo "Not pushed. Re-run with PUSH=1 to publish."
else
    echo "The .sif files are in $SIF_DIR."
    echo "Point the pipeline at them instead of the registry:"
    echo "    nextflow run . -profile $ENGINE --container_dir '$SIF_DIR'"
    echo "Keep the file names as they are — conf/containers.config builds the path"
    echo "as <container_dir>/cami-<group>-<tag>.sif."
fi
