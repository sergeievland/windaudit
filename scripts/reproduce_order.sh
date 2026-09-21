#!/usr/bin/env bash
# Full order remeasurement; runtime-limited trials may differ in unknown status
# between machines. All proven memberships and optima must agree.
set -euo pipefail
cd "$(dirname "$0")/.."
python scripts/select_order.py out/selection_grid_order.json
python scripts/compare_reports.py results/selection_grid_order.json out/selection_grid_order.json
python -m windaudit.order_cli --margin 1 --diagnostic --sweep --out out/paris4/order_report.json
python scripts/diagnostic_order.py --fresh --out out/paris4/order_diagnostic.json
