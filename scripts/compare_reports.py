"""Compare two audit reports on every reproducible section."""

import json
import sys

VOLATILE = {"environment", "runtime_seconds", "inputs"}


def main(expected_path: str, actual_path: str) -> int:
    with open(expected_path, encoding="utf-8") as f:
        expected = json.load(f)
    with open(actual_path, encoding="utf-8") as f:
        actual = json.load(f)
    diverged = [
        key for key in sorted(set(expected) | set(actual))
        if key not in VOLATILE and expected.get(key) != actual.get(key)
    ]
    if diverged:
        print("diverged sections: " + ", ".join(diverged))
        return 1
    print(f"reproduced: {len(set(expected) - VOLATILE)} sections identical")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_reports.py EXPECTED ACTUAL")
    raise SystemExit(main(sys.argv[1], sys.argv[2]))
