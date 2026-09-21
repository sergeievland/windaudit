"""Fetch the public patch surfaces this repository's measurements were made on.

The annotation inputs (``data/paris4``) ship with the repository; the patch
surfaces do not, because they are 236 MB of public data that is better fetched
from its source than vendored. Every file is pinned by SHA-256 in a manifest
that was written at measurement time, so a fetch either reproduces exactly the
bytes the results were computed from, or fails.

    python scripts/fetch_data.py            # download whatever is missing
    python scripts/fetch_data.py --check    # verify what is already on disk

Source: https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/
See THIRD_PARTY_DATA.md for provenance and terms.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/"
WIDE_MANIFEST = ROOT / "results/wide/dataset_manifest.json"
BAND_MANIFEST = ROOT / "results/scanspace_input_hashes.json"
CHUNK = 1 << 20


def planned_files() -> dict[Path, str]:
    """Map every expected local path to its recorded SHA-256."""
    want: dict[Path, str] = {}
    wide = json.loads(WIDE_MANIFEST.read_text(encoding="utf-8"))
    for patch in wide["patches"]:
        for name, meta in patch["files"].items():
            want[ROOT / "data/scanspace_wide" / patch["id"] / name] = meta["sha256"]
    for rel, digest in json.loads(BAND_MANIFEST.read_text(encoding="utf-8")).items():
        want[ROOT / rel] = digest
    return want


def digest_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def remote_url(path: Path) -> str:
    """The upstream URL for a local patch file: <patch id>/<file name>."""
    return BASE + path.parent.name + "/" + path.name


def fetch_one(path: Path, expected: str) -> tuple[Path, str]:
    if path.exists() and digest_of(path) == expected:
        return path, "present"
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".partial")
    last = None
    for _ in range(3):
        try:
            with urllib.request.urlopen(remote_url(path), timeout=120) as response, \
                    temp.open("wb") as out:
                while True:
                    block = response.read(CHUNK)
                    if not block:
                        break
                    out.write(block)
            got = digest_of(temp)
            if got != expected:
                raise ValueError(f"sha256 mismatch: expected {expected}, got {got}")
            temp.replace(path)
            return path, "fetched"
        except Exception as exc:  # retried; reported if the last attempt fails
            last = exc
    temp.unlink(missing_ok=True)
    return path, f"failed: {last}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true",
                        help="verify local files against the manifests; download nothing")
    parser.add_argument("--jobs", type=int, default=8, help="parallel transfers (default 8)")
    args = parser.parse_args(argv)

    want = planned_files()
    total_bytes = json.loads(WIDE_MANIFEST.read_text(encoding="utf-8"))["total_bytes"]

    if args.check:
        missing = [p for p in want if not p.exists()]
        wrong = [p for p in want if p.exists() and digest_of(p) != want[p]]
        print(f"{len(want)} files expected, {len(missing)} missing, {len(wrong)} corrupted")
        for path in (missing + wrong)[:10]:
            print("  ", path.relative_to(ROOT))
        return 1 if missing or wrong else 0

    todo = [(p, d) for p, d in want.items() if not p.exists()]
    print(f"{len(want)} files expected ({total_bytes / 1e6:.0f} MB of patch surfaces); "
          f"{len(todo)} to download")
    failures = []
    done = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(fetch_one, p, d) for p, d in want.items()]
        for future in as_completed(futures):
            path, status = future.result()
            done += 1
            if status.startswith("failed"):
                failures.append((path, status))
            if done % 500 == 0:
                print(f"  {done}/{len(want)}", flush=True)
    if failures:
        print(f"{len(failures)} file(s) failed:")
        for path, status in failures[:10]:
            print("  ", path.relative_to(ROOT), status)
        return 1
    print("complete; every file matches its recorded SHA-256")
    return 0


if __name__ == "__main__":
    sys.exit(main())
