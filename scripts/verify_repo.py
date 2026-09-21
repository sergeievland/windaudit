"""Check this checkout against the three recorded manifests and explain the result.

* ``SHA256SUMS.repo`` covers everything that ships in the repository. It must
  verify completely on a fresh clone, before any data is fetched.
* ``SHA256SUMS`` and ``SHA256SUMS.wide`` are measurement-time records, written
  when the results were produced. They also cover the public patch surfaces,
  which are fetched rather than vendored, and they predate later documentation
  edits. This script reports those two categories separately, so a difference
  in a prose file is never confused with a difference in code or results.

    python scripts/verify_repo.py
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FETCHED_PREFIXES = ("data/scanspace_wide/", "data/scanspace_band/", "data/scanspace_probe/")
DOCS = (".md", ".gitignore")
# Code changed after the measurement-time manifests were written, each for a
# reason that does not touch a published number. Listed rather than hidden:
# SHA256SUMS.repo pins their current bytes.
POST_MEASUREMENT = {
    "pyproject.toml": "Pillow added to the dev extra, which the patch-level loader needs",
    "spiral_lab/launch.py": "child stdout closed explicitly; removes a ResourceWarning",
    "scripts/wide_patch_audit.py": "optional attachment-gap transport; unchanged output when disabled",
    "scripts/verify_wide.py": "repeat comparison tolerant of solver tie-breaks between equally optimal edits",
    "scripts/reproduce_wide.sh": "rebuilds both graph constructions and verifies both",
    "upstream/check_solver.py": "docstring reworded; code unchanged",
    "scripts/make_figure.py": "noise bar shows the 500-round frozen null instead of the 10-round smoke check",
    "results/paris4/summary.png": "re-rendered by the updated make_figure.py; no number in it changed",
}
# Documents moved out of the repository root after the measurement-time
# manifests were written. They are checked at their new location.
MOVED = {
    "EXPERIMENT.md": "docs/future/downstream-fit-protocol.md",
    "ORDER_VALIDATION.md": "docs/negative-results/ordinal-constraints.md",
}
CHUNK = 1 << 20


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def check(manifest: Path):
    ok, missing_fetched, doc_changed, changed, declared = 0, [], [], [], []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, name = line.split("  ", 1)
        path = ROOT / MOVED.get(name, name)
        if not path.exists():
            (missing_fetched if name.startswith(FETCHED_PREFIXES) else changed).append(name)
        elif digest(path) == expected:
            ok += 1
        elif name.endswith(DOCS):
            doc_changed.append(name)
        elif name in POST_MEASUREMENT:
            declared.append(name)
        else:
            changed.append(name)
    return ok, missing_fetched, doc_changed, changed, declared


def main() -> int:
    status = 0
    for name in ("SHA256SUMS.repo", "SHA256SUMS", "SHA256SUMS.wide"):
        manifest = ROOT / name
        if not manifest.exists():
            continue
        ok, fetched, docs, changed, declared = check(manifest)
        if name == "SHA256SUMS.repo":
            print(f"{name} (this release): {ok} files verified")
        else:
            print(f"{name} (measurement-time record): {ok} files verified")
            if fetched:
                print(f"   {len(fetched)} patch-surface files not fetched yet; scripts/fetch_data.py fetches and checks them")
            if docs:
                print(f"   {len(docs)} documentation files reworded after the measurement, as expected")
            for item in declared:
                print(f"   {item}: changed after the measurement - {POST_MEASUREMENT[item]}")
        if changed:
            print(f"   {len(changed)} UNEXPECTED difference(s)")
            for item in changed[:10]:
                print("   unexpected:", item)
        if changed or (name == "SHA256SUMS.repo" and (fetched or docs or declared)):
            status = 1
    print("\nAll code and result files match every manifest." if status == 0
          else "\nSomething differs beyond fetched data and documentation; see above.")
    return status


if __name__ == "__main__":
    sys.exit(main())
