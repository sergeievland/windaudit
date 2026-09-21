# Third-party data

The source code of windaudit is MIT-licensed. The Vesuvius Challenge data this
repository uses is not covered by that grant; it keeps the terms under which the
Vesuvius Challenge publishes it (CC BY-NC 4.0 unless stated otherwise for a
specific asset). Attribution: Vesuvius Challenge / ScrollPrize and the authors
of the annotations and surfaces.

Two kinds of public data are used. The winding annotations and umbilicus ship
with the repository because they are small and are the audit's entire input.
The patch surfaces do not ship: they are fetched from their source by
`scripts/fetch_data.py` and each file is checked against a hash recorded at
measurement time.

## PHercParis4 winding annotations and umbilicus

| File | SHA-256 |
|---|---|
| `relative_windings.json` | `a3243511d4eb91387a9b32f4dbff11514b08c3ae36e9b2a2b8222607b4883ac1` |
| `same_windings.json` | `d9be52c5ebb42853f75f235241cbfd159738f6f34468bcd182523bc91dc91048` |
| `abs_winding.json` | `4e566731f7cbaf8f5ec843de687b3f72f4a784c40b587ebbaf544550902172c1` |
| `umbilicus.json` | `c5f30b0d135c1d333f8170e592079a3d5a636e0071c0dcce560eecebb0ee2602` |

Origin: the public spiral-fitting dataset for PHercParis4,
`https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4`. The copies here
were taken unmodified from the pinned reproduction bundle
`https://github.com/TAUIL-Abd-Elilah/pherc-paris4-absolute-winding-repro`
(`inputs/pcl/` and `inputs/umbilicus.json`), which records the same SHA-256
values for these four files and was assembled for ScrollPrize/villa PR #1626.
`data/paris4/SHA256SUMS` lets `sha256sum -c` confirm integrity.

Provenance of the evidence: these four files are the entire input to the
annotation audit. No surface prediction, fitted spiral checkpoint or trained
model enters any of its steps, and the calibration uses only the
absolute-winding anchors, so every reported annotation check is a statement
about the annotations themselves.

On 17 September 2026 these four files were also downloaded directly from the
public data server and matched the bundled SHA-256 values exactly. The download
record is `results/server_verification.json`; no duplicate data copy is included.
The public dataset did not contain `spiral-scroll.json` on that check (HTTP 404).

## PHercParis4 verified patch surfaces

The patch-graph results additionally use public tifxyz patch surfaces from
`https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/`.
1,639 patch directories were selected and downloaded, 235,651,934 bytes in total, out of
89,237 directories in the public listing; the complete tree was never fetched.

These files are not vendored here. `results/wide/dataset_manifest.json` records
every selected patch with the SHA-256 and byte length of each of its
`meta.json`, `x.tif`, `y.tif` and `z.tif`, and
`results/scanspace_input_hashes.json` does the same for the earlier pilot
corpus. `python scripts/fetch_data.py` downloads exactly that set and refuses a
file whose hash does not match, so a reproduction runs on the same bytes the
measurements were made from, or it fails. `python scripts/fetch_data.py --check`
verifies an existing copy without downloading anything.

Surfaces are read only. Nothing in this repository modifies, redistributes or
republishes them, and no derived surface data is committed — only measurements
taken from them.
