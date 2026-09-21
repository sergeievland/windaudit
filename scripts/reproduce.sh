#!/usr/bin/env bash
# Reproduce every committed result from the pinned inputs.
set -euo pipefail
cd "$(dirname "$0")/.."

(cd data/paris4 && sha256sum -c SHA256SUMS)
python -m pytest -q
python -m windaudit run --inputs data/paris4 --out out/paris4
python scripts/compare_reports.py results/paris4/audit_report.json out/paris4/audit_report.json
# Rebuild the additive order graph and verify its saved measurement records.
# Full, expensive order remeasurement is scripts/reproduce_order.sh.
python scripts/verify_order.py
cmp results/paris4/certified_links.json out/paris4/certified_links.json
cmp results/paris4/review_queue.json out/paris4/review_queue.json
