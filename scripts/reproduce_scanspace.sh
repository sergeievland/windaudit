#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONDONTWRITEBYTECODE=1
mkdir -p out/scanspace_upstream/spiral-fitting
cp upstream/villa/find_inconsistent_windings.py out/scanspace_upstream/spiral-fitting/find_inconsistent_windings.py
patch -d out/scanspace_upstream -p1 < upstream/solver-candidate.patch
python scripts/scanspace_real_graph.py data/scanspace_band \
  upstream/villa/find_inconsistent_windings.py \
  out/scanspace_upstream/spiral-fitting/find_inconsistent_windings.py \
  out/scanspace
python scripts/compare_scanspace.py results/scanspace out/scanspace
