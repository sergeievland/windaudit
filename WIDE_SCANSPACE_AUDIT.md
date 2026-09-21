# The real patch graph, and what it shows about the upstream diagnostic

A patch graph was built from 1,639 public patch surfaces and the unmodified public annotations, with no fitted checkpoint and no GPU, and the pinned upstream functions were run on it unchanged. Built exactly as upstream builds it, the graph has three inconsistent cycles and the upstream diagnostic recommends two annotation edits. A frame-invariance test then showed that two of those three cycles are not properties of the annotations at all: they come from a step upstream's construction never transports, the gap between an annotation point and its attachment on the surface. With that gap transported, one contradiction remains, the graph is frame-invariant, and one of the two recommended edits turns out to break two cycles that are consistent in the corrected graph. Throughout, the upstream solver keeps only the equations on the inconsistent cycles in its model.

The method and the invariance proposition are in [SCANSPACE_TRANSPORT.md](SCANSPACE_TRANSPORT.md).

## Data

The band is `[6000, 18000)` in scan z. It contains all 7,645 annotation points, whose z range is about 6,527.57 to 17,252.55. The selection combines 689 previously inspected metadata records with 1,024 evenly spaced auto-grown catalogue entries; it is a bounded catalogue sample, not an exhaustive search of `verified_patches`. Per-file metadata supplied bounding boxes, both known axis conventions were allowed for before downloading, and actual scan-coordinate quads were checked afterwards. The patch count and a byte bound were fixed before any surface file was fetched, and nothing outside the selection was downloaded.

| Measurement | Value |
|---|---:|
| Metadata records inspected | 1,713 |
| Selected patches | 1,639 |
| Announced upper bound | 236,698,409 bytes |
| Actual selected files | 6,631 |
| Actual total size | 235,651,934 bytes |
| Patches retained after masks, erosion and actual band intersection | 1,637 |
| Valid quads audited for sampling adequacy | 6,223,592 |
| Centre-graph edges audited | 24,287,474 |
| Nonzero integrability residuals | 0 |
| Annotation points attached | 3,523 / 7,645 |

Projection uses the CPU triangle fallback, tolerance 2.5 voxels, default one-cell erosion and largest-area-then-nearest attachment. Bit parity with the optional native surface index is not claimed. No sampling step or attachment tolerance was tuned after observing outcomes.

## The graph as upstream builds it

`classify_pcl`, `build_rel_adjacency` and `solve_min_edge_fix` are extracted from villa commit `2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7` and executed as they are; the scan-space transport supplies the within-patch crossing counts. 335 patches are reached, 330 of them joined by 645 measurable equations, with 7 further edges unmeasurable. The graph has 27 components and 342 fundamental cycles, three of them inconsistent (`results/wide/`).

| | Upstream | Patched |
|---|---:|---:|
| Equations retained in the model | 7 | 645 |
| MILP inequality rows | 14 | 1,290 |
| Solver status | Optimal | Optimal |
| Edits proposed | 2 | 2 |
| Graph consistent after the edits | Yes | Yes |

Both versions propose changing `col106`, points 1530 → 1529, from −1 to −2, together with one of `col264` or `col208`. Which of those last two a solver returns is a tie-break: SciPy 1.17.0 returns `col264` for the upstream version and SciPy 1.17.1 returns `col208`, and both are optimal. Reproduction therefore compares everything except the chosen edit, and checks each run's own edits for optimality and consistency.

## Finding the attachment gap

The frame-invariance test rebuilds the graph under three admissible rotations about the umbilicus. On this graph it fails under all three, and in the same place each time: one equation changes by one winding more than any assignment of per-patch gauges can explain.

Tracing that equation leads to a single annotation point. Upstream counts branch crossings along each annotation at the annotation points' own coordinates, and begins and ends each within-patch strip at the point's attachment on the surface; the segment between the two is transported by neither. For `col109` point 1560 that segment crosses the branch ray:

| | θ (radians) | Position |
|---|---:|---|
| Annotation point `col109` #1560 | 6.28214 | just short of 2π |
| Its attachment on the surface, 1.37 voxels away | 0.00016 | just past 0 |

The two equations that pass through this point — its steps to points 1559 and 1561 — are each off by one winding, and together they create two inconsistent cycles that do not exist once the step across the ray is counted. `tests/test_frame_gauge.py` reproduces the mechanism in six numbers.

This is a defect in the construction, not in the scan frame. Every fitted transform places its branch ray somewhere, and wherever it passes between an annotation point and its attachment the same missing ±1 transport term can occur. For a ray placed uniformly around the axis, the gaps in this graph give 0.087 straddled endpoints on average — about an 8% chance per frame (`results/wide_corrected/attachment_gap_rate.json`). The scan frame hits it; so can any fitted frame, which is why the correction belongs in the construction.

## The graph with the gap transported

With `WIDE_ATTACHMENT_GAP=1` the harness adds the crossing count from each strip's actual end to the annotation point, and from the arrival point to its strip's end, to every edge. Exactly two undirected equations change, both at `col109` point 1560. Nothing else in the graph moves (`results/wide_corrected/`).

| | Upstream | Patched |
|---|---:|---:|
| Inconsistent cycles | 1 | 1 |
| Equations retained in the model | 2 | 645 |
| MILP inequality rows | 4 | 1,290 |
| Edits proposed | 1 | 1 |
| Graph consistent after the edit | Yes | Yes |

This graph passes the frame-invariance test under all three transforms: 17, 119 and 57 of the 645 equation values change, every change is a difference of patch gauges, and every cycle sum and solver output is identical (`results/frame_invariance/`).

It also settles what the earlier `col106` recommendation was. Applied to the corrected graph, that edit raises the number of inconsistent cycles from one to three. It was needed only to cancel the artifact, and following it would break an annotation that agrees with every other equation in the graph.

## The contradiction that remains

| Step | Annotation | Points | Walked | Contribution | Attachment distances |
|---|---|---|---|---:|---|
| patch A → patch B | `col264` | 2882 → 2883 | as stored | +1 | 0.75, 1.80 voxels |
| patch B → patch A | `col208` | 2335 → 2337 | reversed | −2 | 1.67, 0.60 voxels |
| | | | **sum** | **−1** | |

Patch A is `auto_grown_20260521161130783_sel_20260521_161535_29` and patch B `auto_grown_20260521225553961_sel_20260521_225725_10`; the two annotations lie at z ≈ 10,358 and 10,215. On both steps the branch correction is zero and the two strip terms cancel, so the contradiction is carried entirely by the annotations' own winding differences: one says the patches are one winding apart, the other two. It survives all three transforms and the gap correction. Exactly two single edits repair the graph — changing either annotation — so neither is `certain` and both are `possible`. Which of the two is wrong is a question for CT review in VC3D; the certificate gives both annotations, both patches and every attachment distance needed to start it.

`results/wide/certificates.json` and `results/wide_corrected/certificates.json` write out every certificate with the signed contribution of each step and the rule for adding them.

## What the solver comparison shows

In both constructions the upstream solver keeps only the equations that lie on inconsistent cycles — 7 of 645, then 2 of 645 — because `allowed_edge_keys` removes the equations of non-editable edges instead of only forbidding their edit. On these two graphs the editable edges happen to be enough, so its repairs are still valid. The synthetic regression in `upstream/` shows the failure the missing equations allow on other graphs: a repair that closes one cycle and breaks another that was never in the model. The patched solver keeps every equation, so that failure cannot occur.

## Traces 149, 150 and 151

All three traces lie at z = 12,418.1259765625, inside the band. `same:149` has 26 of 32 points attached; `same:150` and `same:151` have none. The data therefore do not decide whether 149 and 151 lie on the same local sheet, and the annotation equality between them is not used as geometric evidence. Hit inventories are in `point_surface_hits.json` and `target_geometry.json`.

## Reproduction

```bash
python scripts/fetch_data.py            # once; every file checked against its recorded SHA-256
bash scripts/reproduce_wide.sh          # both constructions, both solvers, compared with results/
bash scripts/reproduce_invariance.sh    # the frame test on both constructions
```

Set `WIDE_FORCE_RECOMPUTE=1` to bypass the derived surface and transport caches. Once the data is fetched neither command needs the network. The upstream collection identity handling is already correct; no collection-ID fix is included or claimed.
