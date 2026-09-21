# Methods

This document defines every quantity in `audit_report.json` and the protocol that produced the committed PHercParis4 results.

## Inputs

windaudit consumes the four files a spiral-fitting dataset already provides:

* `relative_windings.json`: collections whose points carry integer `wind_a` labels, valid up to one additive constant per collection;
* `same_windings.json`: collections without `wind_a`, each asserting that all of its points lie on one sheet;
* `abs_winding.json`: collections with `metadata.winding_is_absolute`, whose `wind_a` are integer windings in the spiral's branch frame;
* `umbilicus.json`: the scroll axis as control points.

Parsing follows `ScrollPrize/villa/spiral-fitting/point_collection.py`: chain order is point-id order, unannotated points of an annotated collection are skipped, and each collection is classified by content exactly as upstream `classify_pcl` does. A collection whose content contradicts the file it sits in is rejected as an input error rather than reinterpreted.

## Frame and labels

For each point, `r` and `theta` are measured around the umbilicus interpolated at the point's z, with `theta = atan2(y, x)` in `[0, 2pi)` as in upstream `sample_spiral.get_theta`. The integer winding branch cut is at `theta = 0`. For a clockwise-outward scan (`spiral_outward_sense = "CW"`, sense `s = +1`) radius grows with theta, so a sheet followed forward across the cut gains one winding; `"ACW"` gives `s = -1`.

Each collection becomes one node with one unknown integer gauge `g`. A point's integer winding is `W = label + g`, where

* relative collections: `label = wind_a + s * c`,
* same-winding collections: `label = s * c`,
* absolute collections: `label = wind_a`,

and `c` is the signed number of cut crossings accumulated along the chain up to that point (`+1` for a forward crossing on the short arc). Two points `a`, `b` are on the same sheet when

`label_a - label_b + s * cut(a -> b) = g_b - g_a`,

with `cut(a -> b)` the crossing on the short arc from `a` to `b`. Chains with a step wider than 120 degrees are marked seam-unsafe and excluded, since their crossings cannot be counted reliably. Absolute collections are split into angular clusters (gap above 8 degrees), and clusters within 10 degrees of the cut are excluded because the scan-space cut only approximates the fitted spiral's cut.

## Calibration

Absolute anchors from the same sector (`|dz| <= 200`, `|dtheta| <= 8 deg`) with different windings give samples `|dr| / |dw|`. The median is the wrap spacing `D`. Fixed multiples of `D` set the radial tolerances: same-sheet tolerance `r_same = 0.35 D` and crowding floor `0.70 D`. Sector size is 40 voxels in z by 4 degrees. `anchor_order_agreement` reports the fraction of anchor pairs whose radial order matches their winding order. The frozen configuration is serialised and hashed (`config_sha256`).

PHercParis4: `D = 23.30` voxels from 1,053 pairs (10th to 90th percentile 18.0 to 26.9), order agreement 1.0.

## Wrap tilt

Wraps are not vertical. For every point, a local slope `dr/dz` is fitted over points of the same collection that lie on the same sheet (equal branch-corrected label) within 60 voxels in z and 6 degrees. Where the sheet spans less than 10 voxels in z there is no slope. Radial comparisons use the gap at a common z: both points are projected along their own slopes to the midpoint z; a point without a slope stays in place and its partner is projected to it; two slope-less points are compared only within 10 voxels in z.

## Tier 1: single-collection audits

**Order audit.** For every same-sector pair inside one relative or absolute collection, a sign violation is a pair whose projected radial gap exceeds `r_same` in the direction opposite to its label difference; a split is an equal-label pair more than `0.70 D` apart. `shuffled_control` repeats the audit with labels permuted within each relative collection; the gap between the real and shuffled rates is the audit's power on the corpus.

**Same-winding step audit.** For each consecutive step of a same-winding trace, the residual is the radial move minus the median slope of up to three neighbouring steps on each side times the z move. A residual above `0.70 D` is a wrap-scale step against the trace's own trend and is sent to the review queue with both endpoints.

## Tier 2: cross-collection graph

**Links.** Points are binned into sector cells; each cell is compared with itself and its forward neighbours, so every pair within one sector is seen exactly once. Two points from different collections form a link when their projected gap is at most `r_same`. The link's offset is `label_a - label_b + s * cut(a -> b)`.

**Crowding veto.** A cell neighbourhood is skipped when any single collection places two consecutive windings closer than `0.70 D` there, after projection to the neighbourhood's median z and translation to one branch.

**Edges and conflicts.** Links between two collections are grouped by offset. An offset group counts only when it has at least two links (`min_edge_support = 2`); single links are counted as `n_unsupported_links` and dropped. One counted group gives an edge with that offset; two or more give a conflict, which names both collections.

**Cycles.** A spanning forest assigns potentials; each non-tree edge closes one fundamental cycle whose holonomy is `potential_b - potential_a - offset`. A cycle is non-trivial ("effective") unless every node on it has a constant label and every offset on it is zero; trivial cycles cannot fail and are reported separately.

**Repair.** Both repairs share one integer-potential model with big-M rows `|g_b - g_a - offset| <= M * release`. Minimum quarantine minimises, in this order, the number of removed collections, the total support of the union of incident edges (each lost edge counted once), and a canonical choice by node order; minimum edge cut minimises the removed link support and then a canonical choice by edge order. Each stage is solved with HiGHS through `scipy.optimize.milp` with a zero optimality gap, and its optimum is pinned as an equality before the next stage. The canonical stage visits candidates in order and keeps each one out of the removal whenever an optimal repair without it exists, so the mathematical optimum is unique for a given graph; numerical agreement across solver versions must still be tested. Potentials are bounded by `min(sum |offset|, n * max |offset|) + 1`, which contains a solution for every consistent subgraph. A repair is reported as optimal only when every stage succeeds and an independent traversal confirms that the remainder is consistent; tests compare the solver with exhaustive search and check the tie-breaking.

**Certainty of a repair.** The canonical tie-break returns one optimal quarantine, but a reviewer needs to know whether a named collection is implicated by the evidence or merely by that convention. Both objective optima are pinned as equalities, which makes the program's feasible set exactly the set of optimal quarantines; the range of each collection's indicator over that set is then obtained by minimising and maximising it. A collection is `certain` when the range is `[1, 1]`, `excluded` when it is `[0, 0]`, and `possible` otherwise. Both directions are optimisations of a program already known to be feasible rather than feasibility tests of a pinned one, so an infeasible subproblem — which some HiGHS builds report as a generic solve error rather than as infeasibility — never has to be told apart from a failure; a solve that does not reach optimality yields `unknown`. Tests compare the labels with exhaustive enumeration of all optima.

**Export.** Links on edges that survive the minimum edge cut are written as two-point collections without `wind_a`, which upstream loads as same-winding constraints.

## Validation

**Sites.** A boundary inside a collection is a site when evidence independent of that collection connects its two sides: each side needs at least `min_edge_support` links to some partner, and the partners must share a component of the edge graph with the planted collection removed. For a same-winding collection that component must also contain a label-bearing collection, since zero-offset paths alone cannot expose a defect. `site_coverage` is sites over all interior boundaries.

**Planted defects.** At every site, in both directions: a skipped wrap adds `+/-1` to every `wind_a` after the boundary; a sheet switch moves every point after the boundary radially by `D`. The whole pipeline reruns on the mutated corpus. `alarm_rate` is the fraction of plants after which some collection is implicated (conflict or quarantine) that was not implicated on the unmodified corpus; `localization_rate` requires the planted collection itself to be newly implicated; when several collections on one failing cycle are equally good repairs, the canonical tie-break decides which is named, so this rate reflects that convention. `shortlist_rate` is free of it: it requires the planted collection to be newly on a failing cycle or in a conflict, with `mean_shortlist_size` the number of collections a reviewer would open. `certain_localization_rate` is free of it in the stronger sense: it counts trials in which the planted collection is removed by *every* optimal quarantine, and `membership_counts` gives the distribution of those labels over the trials.

**Noise null.** Every relative and same-winding point receives independent Gaussian noise with a standard deviation of 1 voxel on each axis; `rounds_with_alarm` counts rounds that implicate any new collection.

**Sensitivity sweep.** Matching and repair rerun for `r_same` in `{0.25, 0.35, 0.45} D`, sector width in `{3, 4, 6}` degrees and sector height in `{30, 40, 60}` voxels. `robustness` is, per collection, the fraction of the 27 settings that implicate it.

## Protocol for the committed results

1. Inputs were pinned by SHA-256 (`data/paris4/SHA256SUMS`, provenance in `THIRD_PARTY_DATA.md`).
2. The edge-support rule and `r_same` were chosen on the grid in `results/selection_grid.json` (seed 777): among settings with at most 2 of 20 noise rounds raising an alarm, the one with the highest planted-defect alarm rates. The selected setting was then frozen in `windaudit/calibrate.py`.
3. The committed report was produced once with the frozen configuration and seed 20260917 (`python -m windaudit run --inputs data/paris4 --out results/paris4`); `scripts/reproduce.sh` regenerates it and `scripts/compare_reports.py` confirms the reproducible sections agree.

## PHercParis4 in numbers

| Quantity | Value |
|---|---|
| Nodes (non-empty collections, absolute clusters) | 384 |
| Seam-unsafe nodes | 0 |
| Nodes crossing the branch cut | 4 |
| Order pairs tested / violations | 11,447 / 0 |
| Shuffled order-violation rate | 49.4% (48.0 to 51.0%) |
| Trace steps tested / flagged | 3,536 / 2 |
| Links / edges / conflicts | 2,371 / 98 / 0 |
| Single links dropped by the support rule | 7 |
| Links across the branch cut | 1 |
| Cells / vetoed as crowded | 1,109 / 412 |
| Graph nodes / components / largest | 89 / 18 / 20 |
| Cycles / non-trivial / failing | 27 / 3 / 0 |
| Relative ladders linked to another collection | 5 of 254 |
| Same-winding traces linked | 84 of 125 |
| Skipped-wrap sites / plants / alarms | 7 / 14 / 14 |
| Sheet-switch sites / plants / alarms | 48 / 96 / 55 |
| Noise rounds with alarm, in-report null (seed 20260917) | 0 of 10 |
| Noise rounds with alarm, frozen 500-round null (seed 20260921) | 86 of 500 (17%; 95% interval 14–21%) |
| Largest sweep robustness of any collection | 0.19 |

## Beyond the annotations

This document defines the annotation audit, whose only inputs are the four point-collection files. The patch-level work — measuring upstream's integer winding transport in scan coordinates rather than through a fitted checkpoint, building a real patch graph from public tifxyz surfaces, and running the pinned upstream repair solver on it — is defined in [SCANSPACE_TRANSPORT.md](SCANSPACE_TRANSPORT.md), with its results in [WIDE_SCANSPACE_AUDIT.md](WIDE_SCANSPACE_AUDIT.md). Those two paths share the frame convention and the branch-transport rule defined above; they do not share inputs, and neither changes any quantity reported here.

## Scope

The audit certifies the mutual consistency of annotations. Its reach is set by the corpus: a constraint can be checked only where independent annotations overlap it, which the coverage and site figures report directly. The effect of certified or repaired constraints on a downstream spiral fit is outside this release.

## Evidence boundaries

Cycle contradictions are exact only conditional on the inferred same-sheet correspondences and branch model. A contradiction does not prove that the human annotation, rather than an inferred link, is wrong. Passing the graph audit does not certify physical correctness or unsampled regions. Absolute collections are represented in the gauge graph; this is not an independently held-out absolute accuracy measurement.

The selected configuration uses the same corpus for calibration and model selection. Seed changes do not create a held-out corpus. The reported 14/14 alarm result covers only 7 of 1,823 eligible-type boundaries (0.38%); 55/96 sheet-switch alarms cover 48 of 5,288 boundaries (0.91%). Trials in opposite directions at one boundary are not independent samples. The ten-round null run inside every report is a smoke check, not an estimate. The frozen null in `results/frozen_null.json` runs the same perturbation 500 times with the configuration exactly as frozen and a seed that played no part in selecting it: 86 rounds implicate at least one collection the clean run does not, a rate of 17% with a 95% Clopper–Pearson interval of 14–21%. Only 10 of the 384 collections are ever implicated, and three same-winding traces account for 65 of the 86 alarm rounds. The null says nothing directly about the field false-positive rate under a different noise model.

Edge-loss indicator variables encode the OR of the removed-endpoint indicators. This corrects double-counting in the secondary quarantine objective. Canonicalisation accepts only solver status 0 or proven infeasibility; other statuses invalidate optimality. Certainty labels are read from optimisation, not from infeasibility detection, and a solve that does not reach optimality is reported as `unknown` rather than as an answer. Export additionally rechecks the retained graph and fails closed. The sensitivity sweep and robustness values were regenerated after this correction.
