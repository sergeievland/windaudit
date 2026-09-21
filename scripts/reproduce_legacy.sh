#!/usr/bin/env bash
# Reproduce every committed result from the pinned inputs.
set -euo pipefail
cd "$(dirname "$0")/.."

(cd data/paris4 && sha256sum -c SHA256SUMS)
python -m pytest -q
python -m windaudit run --inputs data/paris4 --out out/paris4
python scripts/compare_reports.py results/paris4/audit_report.json out/paris4/audit_report.json
