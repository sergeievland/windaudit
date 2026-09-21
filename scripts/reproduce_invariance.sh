#!/usr/bin/env bash
# The frame-invariance test on both graph constructions. The graph as upstream
# builds it is expected to FAIL (the attachment gap); the corrected graph is
# expected to pass. Run scripts/reproduce_wide.sh first.
set -euo pipefail
cd "$(dirname "$0")/.."
python scripts/frame_invariance.py out/wide_reproduced out/frame_invariance/as_upstream_builds_it || true
WIDE_ATTACHMENT_GAP=1 python scripts/frame_invariance.py out/wide_corrected_reproduced out/frame_invariance/attachment_gap_transported
python - <<'PY'
import json
for name, expected in (("as_upstream_builds_it", False), ("attachment_gap_transported", True)):
    got = json.load(open(f"out/frame_invariance/{name}/frame_invariance.json"))["invariant"]
    ref = json.load(open(f"results/frame_invariance/{name}.json"))["invariant"]
    assert got == ref == expected, (name, got, ref)
    print(f"{name}: invariant={got}, as committed")
PY
