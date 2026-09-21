#!/usr/bin/env bash
# Rebuild the real patch graph from the fetched public surfaces, both as
# upstream builds it and with the attachment gap transported, run the pinned
# upstream solver before and after the candidate patch on each, and compare
# every output with the committed results. Needs scripts/fetch_data.py once;
# runs offline afterwards.
set -euo pipefail
cd "$(dirname "$0")/.."
python scripts/fetch_data.py --check
mkdir -p out/wide_solver/spiral-fitting
cp upstream/villa/find_inconsistent_windings.py out/wide_solver/spiral-fitting/find_inconsistent_windings.py
patch --batch --quiet -p1 -d out/wide_solver < upstream/solver-candidate.patch
AFTER=out/wide_solver/spiral-fitting/find_inconsistent_windings.py
python scripts/wide_patch_audit.py data/scanspace_wide upstream/villa/find_inconsistent_windings.py "$AFTER" out/wide_reproduced
python scripts/verify_wide.py results/wide out/wide_reproduced
WIDE_ATTACHMENT_GAP=1 python scripts/wide_patch_audit.py data/scanspace_wide upstream/villa/find_inconsistent_windings.py "$AFTER" out/wide_corrected_reproduced
python scripts/verify_wide.py results/wide_corrected out/wide_corrected_reproduced
python scripts/make_certificates.py out/wide_corrected_reproduced
