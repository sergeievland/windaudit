# Winding transport without a checkpoint

Upstream's winding diagnostic needs a fitted spiral checkpoint, so it cannot be run until a GPU fit has already been paid for. This document shows which part of what it computes depends on the fit and which does not, proves that the part a diagnostic relies on — its contradictions — is the same in every admissible frame, tests that on real patch surfaces, and describes the CPU implementation that follows.

The runs built on this path — the 1,639-patch graph, the upstream solver comparison and the attachment-gap defect — are in [WIDE_SCANSPACE_AUDIT.md](WIDE_SCANSPACE_AUDIT.md).

## What the integer is, and what it is not

`strip_winding_delta` in `find_inconsistent_windings.py` transforms a within-patch path into spiral space, computes θ along it, and sums the seam corrections from `get_theta_crossing_step_adjustments`. Those corrections are `±1 · dr_per_winding`, and the next line divides the sum by `dr.detach()`. The factor cancels exactly: the returned `delta_windings` is the signed count of θ=0 crossings along the path. The spiral's radial scale plays no part in it.

The transform itself does. θ is measured after it, and a transform that rotates points about the axis moves the branch ray relative to the scroll. A path's crossing count in scan coordinates is therefore not, in general, the same number upstream would compute through a checkpoint, and no such equality is claimed. What does hold is the statement below, which is the one a contradiction finder needs.

## Proposition: contradictions do not depend on the frame

*Setting.* The graph's equations are built from pieces of path, each with a signed branch-crossing count: within-patch strips, steps along an annotation between consecutive points, and — once transported — the attachment gap between an annotation point and the place its strip ends on the surface. An edge from patch P to patch R walks: the strip on P from P's entry point to the end of the strip at the departure attachment, the gap to the departure annotation point, the annotation's own steps to the arrival point, the gap to the arrival attachment, and the strip on R back to R's entry point. Its equation is the annotation's winding difference plus the signed count along that walk.

*Admissible transform.* A continuous map T of the working region into spiral space that sends the umbilicus onto the spiral axis, sends no other point onto the axis, and carries a small loop around the umbilicus to a loop that winds once around the axis.

*Sampling condition.* Consecutive samples of every piece turn by less than π about the axis, both before and after T.

*Claim.* There is an integer-valued function k on points such that, for every piece from a to b, the crossing count after T equals the count before T plus k(b) − k(a). Consequently every closed walk has the same count in both frames; every edge equation changes by exactly k(entry of R) − k(entry of P), a difference of per-patch gauges; and so every cycle sum, the set of inconsistent cycles, whether any subgraph is consistent, the minimum number of edits and the set of optimal edit sets are identical in both frames. Only the suggested new values shift, by the gauge.

*Proof.* Let D be the working region without the umbilicus. The angle difference ψ(p) = θ(T(p)) − θ(p) is a continuous map from D to the circle. Its winding around the generator of π₁(D) — a small loop around the umbilicus — is the winding of T's image loop minus one, which is zero by admissibility, so ψ lifts to a continuous real function F on D. Put k(p) = (θ(p) + F(p) − θ(T(p))) / 2π, with both angles taken in [0, 2π); k is integer-valued. Under the sampling condition, the jump rule counts the crossings of the continuous path, and for a path from a to b that count is (Φ − θ(b) + θ(a)) / 2π, where Φ is the total continuous change of angle. The total changes in the two frames differ by F(b) − F(a); substituting gives count_T = count + k(b) − k(a). Along an edge's walk the k-differences of consecutive pieces telescope: strip on P, k(end) − k(entry_P); gap, k(point) − k(end); annotation steps, k(arrival point) − k(departure point); gap, k(arrival end) − k(arrival point); strip on R reversed, k(entry_R) − k(arrival end). The total is k(entry_R) − k(entry_P). Adding a potential difference to every equation leaves every cycle sum unchanged and maps the consistent subsets, and the edits that produce them, onto themselves. ∎

*When it fails.* If a piece of the walk is not transported, its k-difference is missing and the telescoping leaves a residual. Upstream does not transport the attachment gap: it counts annotation steps at the annotation points' own coordinates but begins and ends strips at their attachments. Whenever the branch ray passes between a point and its attachment in one frame and not in another, that edge's equation changes by one winding more than any gauge explains. [WIDE_SCANSPACE_AUDIT.md](WIDE_SCANSPACE_AUDIT.md) shows this happening on real data, and shows that transporting the gap removes it.

## The proposition, tested on the real graph

`scripts/frame_invariance.py` rebuilds the whole 645-equation patch graph under synthetic admissible transforms — rotations about the umbilicus by a smooth angle α(p), with |∂α/∂θ| < 1 so each circle maps bijectively — and compares every output with the untransformed run.

| Transform | α(p) | Equation values changed | Changes not explained by patch gauges | Cycle sums | Upstream and patched solver outputs |
|---|---|---:|---:|---|---|
| Varying with z | 1.1 sin(z/700 + 0.3) | 17 of 645 | **0** | identical | identical |
| Varying with radius | 1.3 cos(r/150) | 119 of 645 | **0** | identical | identical |
| Varying with z, radius and angle | 0.9 sin(z/700 + 0.3) + 0.7 cos(r/150) + 0.45 sin θ | 57 of 645 | **0** | identical | identical |

That is the graph with the attachment gap transported. The same test on the graph as upstream builds it fails under every one of the three transforms, each time on one equation: the one the attachment gap of `col109` point 1560 shifts (`results/frame_invariance/`). The identity itself is also checked directly, independent of any data, in `tests/test_frame_gauge.py`: eighty random smooth paths under random admissible rotations, and one untransported gap that is required to break closure.

The proposition is what makes a comparison against a real fitted transform unnecessary for contradictions — useful, since no public spiral checkpoint was available. Individual equation values do differ from upstream's by per-patch gauges, exactly as the proposition says.

## The sampling gate

The proposition assumes that no sampled step turns by π or more. On a single patch, a count that depends on the route rather than only its endpoints would mean exactly that assumption failing, so route independence was tested on real data before the transport was used for anything.

| Check on real patch surfaces | Result |
|---|---:|
| Endpoint pairs, each routed five different ways through the valid-quad graph | 2,400 |
| Pairs whose five routes disagreed | **0** |
| Closed walks, each required to transport to zero | 9,600 |
| Closed walks with a nonzero result | **0** |

The five route variants were medial weights 4, 0 and 20, plus two seeded positive perturbations of weight 4. Endpoints and the valid 8-neighbour quad graph stayed fixed; each alternative was compared with the baseline, and the alternative plus the reversed baseline formed the closed walk. The sampling step was the pinned upstream default of one grid cell and was never tuned.

The bulk check is stronger than the sampled one. Rather than sampling routes, it audits every elementary edge of the sampled centre graph for exact integrability, which certifies *every* closed walk on those surfaces at once. Zero is the expected answer — θ around the umbilicus is single-valued, and no patch encloses the axis — so the audit is a check that one-cell sampling is fine enough everywhere, not a discovery:

| Bulk integrability audit | Result |
|---|---:|
| Patches audited | 1,637 |
| Valid quads | 6,223,592 |
| Centre-graph edges checked | 24,287,474 |
| Nonzero residuals | **0** |
| Independent cycles certified, 46-patch pilot band alone | 1,429,048 |

Bilinear lifting and upstream's invalid-sample dropping are preserved, and dropped samples are counted explicitly rather than silently. The gate is not vacuous: the test suite includes a synthetic surface enclosing the axis whose closed walk transports to +1, and the gate is required to fail on it.

An unseen surface is not covered by a finite audit. The gate therefore runs as part of the pipeline, before projection or any solver call, rather than once as a claim.

## Implementation

The transport reuses upstream's own machinery wherever it exists: the valid-quad Dijkstra route with medial weighting, the polyline sampling at a fixed step, the bilinear lift and the invalid-sample rule. What changes is that θ is taken in scan coordinates around the umbilicus instead of in spiral space, so no checkpoint and no Torch are loaded. The code is NumPy and SciPy on one CPU.

With `WIDE_ATTACHMENT_GAP=1` the harness also transports the attachment gap: for every edge it takes the last valid sample of each strip — the point where the strip really ends, which differs from the bilinear lift of the attachment when that lift falls on an invalid cell — and adds the crossing count from there to the annotation point, and from the arrival point to its strip's end, to the edge's winding difference. BFS entry points do not depend on winding values, so the tree is unchanged and only accumulated windings are recomputed. With the variable unset the harness builds the graph exactly as upstream does.

Point attachment uses the pinned `Patch.project` triangular fallback on actual valid faces, with conservative KD-tree pruning and largest-area-then-nearest selection. Equivalence to the optional native `vc_spiral` surface index is not claimed, and a reconstructed surface can change which annotation points attach. Patch selection keeps surfaces whose *valid quad extents* intersect the requested band; this deliberately differs from upstream `load_tifxyz`'s vertex-only z prefilter, which can miss a two-voxel band falling between grid rows. Upstream's default one-cell binary erosion and per-patch overrides are respected.

The upstream `build_rel_adjacency` and `solve_min_edge_fix` functions are executed from the pinned source AST, so the comparison runs upstream's real code rather than a reimplementation of it.

## Pilot run

The method was first exercised on a deliberately small graph: a band one z-slice wide around z = 12418, 42 patches selected from 341 inspected metadata records, 10.5 MB. Of 7,645 annotation points, 357 attached; 67 collections produced 22 reached patches and 53 measurable equations, with no unmeasurable edges.

That graph is consistent, so neither solver version proposed a repair, and the pilot on its own demonstrates operation rather than a difference in outcome. It did already expose the behavioural gap that the wider run then measured properly: in the restricted mode upstream actually uses, the unpatched solver considered **0** of the 53 real equations and reported "no edges", while the patched solver considered all 53 and solved to optimality.

The pilot also failed to settle the `same:149` / `same:151` question it was aimed at: 26 of 32 points of `same:149` attached to a surface, but `same:150` and `same:151` attached none, so there was no surface path to compare. The wider run did not change that. Lying inside the z band is not the same as being covered by a reconstructed surface, and the annotation equality is not used as geometric evidence for itself.

## Reproduction

```bash
python scripts/fetch_data.py            # public patch surfaces, checked against recorded hashes
python scripts/verify_repo.py           # every shipped file against its manifest
bash scripts/reproduce.sh               # annotation audit
bash scripts/reproduce_wide.sh          # sampling gate, both graph constructions, both solvers
bash scripts/reproduce_invariance.sh    # the frame-invariance test on both constructions
python -m pytest -q tests/test_frame_gauge.py   # the proposition's identity, no data needed
```

`reproduce_wide.sh` runs the exhaustive edge gate before any projection or solver call, so a failed gate stops the run instead of producing numbers. Set `WIDE_FORCE_RECOMPUTE=1` to bypass the derived surface and transport caches. Neither reproduction command needs network access once the data is fetched.

Source: ScrollPrize/villa commit `2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7`. Public data: `https://dl.ash2txt.org/datasets/spiral_datasets/PHercParis4/verified_patches/`. Per-file hashes for every input accompany the evidence in `results/`.
