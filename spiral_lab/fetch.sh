#!/usr/bin/env bash
set -euo pipefail
data="${1:-$HOME/vesuvius/data/PHercParis4}"
mkdir -p "$data"
rclone copy :http: "$data" \
  --http-url https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/ \
  --transfers 4 --checkers 4 --retries 5 -P
# The public directory did not include this required file on 2026-09-17.
# Physical constants follow the official PHercParis4 tutorial. Smoke mode
# disables normals/spacing/winding models; their scales need separate checks.
python3 - "$data" <<'PYCODE'
import json,sys
from pathlib import Path
p=Path(sys.argv[1])/'spiral-scroll.json'
if not p.exists():
    with p.open('x') as f:
        json.dump(dict(schema_version=1,name='PHercParis4',voxel_size_um=9.6,
                       spiral_outward_sense='CW'),f,indent=2)
PYCODE
echo 'FETCH DONE (smoke-ready metadata; volume scales still need verification for a full fit)'
