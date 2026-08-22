# =============================================================================
# cami-qc
# -----------------------------------------------------------------------------
# BUILD FROM THE REPOSITORY ROOT, not from this directory:
#     docker build -f containers/qc/Dockerfile -t ghcr.io/alihkz94/cami-qc:$TAG .
# The stats image needs env/requirements-stats.txt, so every image uses the
# repository root as its build context and they stay consistent with each other.
#
# Prefer the build script, which reads the tag from nextflow.config:
#     bash scripts/containers/build_images.sh docker qc
# =============================================================================
FROM mambaorg/micromamba:2.3.2-debian12-slim

LABEL org.opencontainers.image.title="cami-qc" \
      org.opencontainers.image.description="CaMi pipeline — qc stage" \
      org.opencontainers.image.source="https://github.com/alihkz94/CaMi" \
      org.opencontainers.image.licenses="MIT"

USER root
COPY containers/qc/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm -f /tmp/environment.yml

# Nextflow launches a container as `<image> /bin/bash -c '...'`, which bypasses
# the micromamba entrypoint that would normally activate the environment. Putting
# the prefix on PATH is therefore required, not a convenience.
ENV PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Singularity runs as the invoking user, whoever that is, and cannot switch.
# Leaving the image on the unprivileged micromamba user is what upstream intends
# and keeps the Docker path from running the pipeline as root.
USER mambauser
