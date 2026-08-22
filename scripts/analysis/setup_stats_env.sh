#!/usr/bin/env bash
# =============================================================================
# setup_stats_env.sh — build the Python environment for --step STATS
# -----------------------------------------------------------------------------
# PURPOSE:
#   Create a venv at <dataRoot>/envs/stats holding the exact versions in
#   env/requirements-stats.txt. Run it ONCE per machine, before --step STATS.
#
# WHY NOT JUST INSTALL INTO THE `cancer` CONDA ENV:
#   That env is what bwa-mem2, Kraken 2, Kaiju, fastp and samtools run from.
#   A 563-sample profiling run uses it for days. Adding a scientific Python
#   stack to it could upgrade a shared library underneath a running job. A venv
#   is a separate prefix, so the conda env is never touched.
#
#   The venv is built FROM the conda env's interpreter, so the Python version
#   matches (3.12.13) and no second interpreter has to be installed.
#
# IDEMPOTENT: safe to re-run. It verifies and repairs rather than rebuilding.
#
# USE:
#   bash scripts/analysis/setup_stats_env.sh
#   bash scripts/analysis/setup_stats_env.sh --recreate   # discard and rebuild
#
# VERSION: 1.0 (2026-08-20)
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE="$(cd "$HERE/../.." && pwd)"          # the repository
DATA="$(cd "$CODE/.." && pwd)"             # the data root, one level up
VENV="$DATA/envs/stats"
REQ="$CODE/env/requirements-stats.txt"

BASE_PY="${BASE_PY:-/home/ahakimzadeh/miniforge3/envs/cancer/bin/python}"
[[ -x "$BASE_PY" ]] || BASE_PY="$(command -v python3)"

[[ -f "$REQ" ]] || { echo "ERROR: $REQ not found" >&2; exit 1; }

if [[ "${1:-}" == "--recreate" ]] && [[ -d "$VENV" ]]; then
    echo "removing existing $VENV"
    rm -rf "$VENV"
fi

echo "=============================================="
echo " base python : $BASE_PY  ($("$BASE_PY" -V 2>&1))"
echo " venv        : $VENV"
echo " requirements: $REQ"
echo "=============================================="

if [[ ! -x "$VENV/bin/python" ]]; then
    echo "creating venv..."
    "$BASE_PY" -m venv "$VENV"
fi

"$VENV/bin/python" -m pip install --quiet --upgrade pip
echo "installing pinned requirements (this can take a few minutes)..."
"$VENV/bin/python" -m pip install --quiet -r "$REQ"

echo
echo "--- verifying ---"
"$VENV/bin/python" - <<'PY'
import importlib.metadata as md
import sys
need = {
    "polars": "1.43.2", "pyarrow": "25.0.1", "numpy": "2.5.2",
    "scipy": "1.18.0", "statsmodels": "0.14.6",
}
bad = 0
for pkg, want in need.items():
    try:
        got = md.version(pkg)
    except Exception:
        print(f"  MISSING {pkg}")
        bad += 1
        continue
    mark = "ok " if got == want else "DIFFERS"
    if got != want:
        bad += 1
    print(f"  {mark:8s}{pkg:14s}{got:12s}(pinned {want})")

# The transform and the test must work, not merely import.
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
x = np.array([[1.0, 2.0, 3.0, 4.0]])
c = np.log(x) - np.log(x).mean(axis=1, keepdims=True)
assert abs(c.sum()) < 1e-12, "CLR failed"
assert multipletests([0.01, 0.2], method="fdr_bh")[1].size == 2, "BH failed"
wilcoxon([1.0, 2.0, 3.0, -0.5])
print("  functional check: CLR, Benjamini-Hochberg and Wilcoxon all work")
sys.exit(1 if bad else 0)
PY

echo
echo "DONE. Point the pipeline at it with:"
echo "  --stats_python $VENV/bin/python"
echo "(that is already the default in nextflow.config)"
