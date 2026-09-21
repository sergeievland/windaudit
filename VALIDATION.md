# Validation log

What was checked, and what was not. Last updated 21 September 2026.

## Test suite and reproduction

The test suite has been run in three environments that share no installed package versions:

| | Environment A | Environment B | Environment C |
|---|---|---|---|
| Python | 3.12.14 | 3.11.15 | 3.10.12 |
| NumPy | 2.3.5 | 2.4.4 | 2.2.6 |
| SciPy | 1.17.0 | 1.17.1 | 1.15.3 |
| Annotation and solver tests (238) | pass | pass | pass |
| Frame-identity tests added in this build (81) | not run | pass | pass |
| Suite with warnings as errors | — | pass | — |
| Full audit with sweep and planted defects | rerun | rerun | not run |

The full PHercParis4 audit — including the planted defects, the noise null and the 27-setting tolerance sweep — was rerun in environments A and B. All 17 compared report sections are identical between them, and `certified_links.json` and `review_queue.json` are byte-identical (`sha256` `e865fcaa4a02cd78…` and `067f30a026daccde…`, recorded in `results/export_reproduction.json`).

Environment C matters because its SciPy predates the HiGHS build the others use, and an earlier version of this tool did produce solver-dependent results there. Its test suite passes and its `certified_links.json` is byte-identical to the other two; the full sweep was not run to completion in that environment, so no claim is made about the sweep-derived sections there.

Runtime, environment and input basenames are excluded from the comparison by design; it is a semantic comparison of the reported quantities, not a claim of byte identity for the whole report.

## Scan-space transport and the patch graph

The checkpoint-free transport was gated before use. On real patch surfaces, 2,400 point pairs were each routed five different ways through the valid-quad graph and 9,600 closed contours were transported: zero disagreements and zero nonzero contours. Every elementary edge of the sampled centre graph was then audited for exact integrability across 1,637 retained patches, 24,287,474 edges and 6,223,592 valid quads, with zero nonzero residuals. Zero is expected, since θ around the umbilicus is single-valued and no patch encloses the axis; the audit confirms that one-cell sampling never lets an angle step reach π. A synthetic surface enclosing the axis is required to fail it.

The frame-invariance proposition was tested on the real 645-equation graph under three admissible rotations about the umbilicus (`results/frame_invariance/`). With the attachment gap transported, 17, 119 and 57 equation values change and every change is a difference of per-patch gauges; cycle sums, the inconsistent cycle, the equations retained by both solver versions and the edits they return are identical to the untransformed run. On the graph as upstream builds it, the same test fails under all three transforms at one equation, which is how the attachment gap was found. The identity behind the proposition is also tested without data in `tests/test_frame_gauge.py`.

Both constructions were rebuilt from the fetched surfaces in a second environment (Python 3.11.15, NumPy 2.4.4, SciPy 1.17.1) and compared with the committed results: every output agrees except which of two equally optimal edits the solver returns. On the graph as upstream builds it, SciPy 1.17.0 returned `col264` and SciPy 1.17.1 returns `col208`; both are optimal and both leave the graph consistent. `scripts/verify_wide.py` therefore compares everything but the chosen edit, and checks each run's own edits for optimality and consistency instead.

Independently of the shipped checks, every contradiction and every solver proposal was re-verified from the `measured_edges.json` files with a separate traversal. As upstream builds it, the graph fails to close in three independent cycles and either solver's proposal closes all of them. With the gap transported, one cycle fails; exactly two single edits repair it, `col264` or `col208`. Applying the upstream diagnostic's other recommendation, the `col106` edit, to the corrected graph raises the number of failing cycles from one to three. The rate estimate for the attachment gap in `results/wide_corrected/attachment_gap_rate.json` uses each strip's actual last valid sample, not the bilinear lift of the attachment, which falls on an invalid cell for 12 endpoints of this graph.

The patch attachment uses the CPU triangle fallback with a fixed tolerance and one-cell erosion. Bit parity with the optional native surface index is not claimed, and a reconstructed surface can affect which annotation points attach.

## Measurement noise

The ten-round null inside every report is a smoke check. `results/frozen_null.json` repeats it 500 times with the frozen configuration (hash `eb22ab42a4ad0e36…`, identical to the report's) and a seed not used in selection: 86 of 500 rounds implicate at least one collection the clean run does not, 17% with a 95% Clopper–Pearson interval of 14–21%. Only 10 of 384 collections are ever implicated, and `same:19`, `same:6` and `same:7` account for 65 of the 86 alarm rounds. The 20-round grid on which the edge-support rule was selected had seen 2 alarms in 20; the 500-round null is the calibrated figure, and it is the one the summary figure shows.

## Inputs

The four bundled input files were downloaded again from `https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/` on 17 September 2026 and match the bundled copies by SHA-256 (`results/server_verification.json`). `spiral-scroll.json` was not present in that public directory on that date.

## Solver

Both repair objectives are compared against exhaustive enumeration on 24 random graphs, covering the complete lexicographic objective — removal cardinality, then lost link support, then the canonical choice — rather than cardinality alone.

The certainty labels are checked the same way: for each random graph the tests enumerate *every* optimal quarantine by brute force and require the reported label — removed by all optima, by some, or by none — to match exactly, and require the canonical repair to be one of the enumerated optima. Agreement is exact on every graph tested.

That comparison found a real defect in the 0.1.0 objective, recorded as a regression in `results/quarantine_counterexample.json`: on that graph 0.1.0 removed two collections losing 24 units of link support, while the optimum at the same cardinality loses 23. The cause was that support was accumulated per endpoint, so an edge between two removed collections was charged twice. The edge-loss indicators now encode the OR of the endpoint indicators and each lost edge is charged exactly once.

Canonicalisation accepts only solver status 0 or proven infeasibility (status 2); a time limit or solver error now invalidates the optimality claim instead of being silently read as "this candidate must be removed". Export refuses a repair that was not proven optimal and independently rechecks the retained graph before writing. Malformed inputs — non-finite coordinates or labels, wrong coordinate dimension, colliding normalised identifiers, non-integer or negative graph weights — are rejected rather than coerced.

Correcting the objective changed no headline result. Comparing the 0.1.0 and 0.2.0 reports field by field, every reported count is identical; only two sweep-derived robustness values moved, because some secondary repairs inside the 27-setting sweep now resolve differently (`results/upgrade_comparison.json`).

## Candidate upstream patch

`upstream/solver-candidate.patch` applies cleanly to `spiral-fitting/find_inconsistent_windings.py` at the pinned revision `2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7` and at `main` as of 20 September 2026 (`23c3d759`), where the file is byte-identical. On both, `upstream/check_solver.py` reports both regressions failing before the patch, exit code 1, and both passing after it, exit code 0.

Its synthetic regressions are complemented by the real-data runs: the upstream solver retains 7 of 645 equations on the graph as upstream builds it and 2 of 645 with the attachment gap transported; the patched solver retains all 645 on both. Both versions reach optimality and propose valid repairs on these graphs, where the editable edges happen to suffice; the patched version keeps every equation, so the failure shown by the synthetic regression cannot occur in it.

## Launcher

`spiral_lab/launch.py` is covered for success and failure logging, refusal to overwrite a run directory, dry run and environment isolation, using a subprocess fixture. Its smoke overrides are accepted by the pinned upstream `Config` class.

## Outside this release

Windows or WSL installation, the CUDA fit, ink rendering, held-out checkpoint evaluation (designed in advance in [docs/future/downstream-fit-protocol.md](docs/future/downstream-fit-protocol.md)), CT review of the findings, and a numerical comparison of the scan-space transport against a fitted checkpoint — which the invariance proposition makes unnecessary for contradictions, and for which no public spiral checkpoint was available. Cross-machine identity is claimed for the exports and the semantic report, not for every byte of the full report.
