# spiral_lab

Pinned setup, data fetch and run launcher for the downstream fit comparison designed in [`docs/future/downstream-fit-protocol.md`](../docs/future/downstream-fit-protocol.md). It is kept with the audit so that, if that comparison is ever run, its provenance and the audit's provenance live in one repository.

`setup.sh` clones ScrollPrize/villa at the pinned revision `2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7` into a fresh directory, refusing to overwrite an existing one, and verifies CUDA before reporting success. `fetch.sh` mirrors the public PHercParis4 spiral dataset and writes the `spiral-scroll.json` that the public directory did not contain on 17 September 2026, with the physical constants taken from the official tutorial. `launch.py` runs one fit with explicit overrides, records the exact configuration, input hashes and full logs, and refuses to overwrite a previous run directory. `smoke.json` is a 300-step, 1,000-slice operational check: it tests that the pipeline runs, not that it converges.

What is covered by the test suite: the launcher's orchestration, override handling, failure logging, dry run and environment isolation, and the fact that the smoke overrides are accepted by the pinned upstream `Config` class. What is not: the installation itself and the CUDA fit, neither of which has been executed here. `launch.py` hashes the annotation JSONs and `uv.lock`, not the patch trees, and records that limitation in its own output.

Nothing in this repository's measured results depends on this directory. The checkpoint-free path in [`SCANSPACE_TRANSPORT.md`](../SCANSPACE_TRANSPORT.md) reaches real patch data without a fit at all, which is why the results that exist were obtained without ever running these scripts.
