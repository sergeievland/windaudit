# A rejected extension, and why it was rejected

Papyrus wraps do not cross, so it is tempting to turn radial order into winding order and harvest constraints wherever two annotations share a sector rather than only where they share a wrap. On this corpus that would have multiplied the checkable structure roughly tenfold — from 89 collections in the graph to 177, and from 7 checkable skipped-wrap boundaries to 142.

It does not survive contact with the data. Every one of the seven candidate margins failed the noise control, the unmodified corpus is declared broken even at a three-wrap margin — 18 collections would have to be quarantined with nothing planted, and there is an explicit smooth counterexample showing that non-intersection does not imply radial ordering between different angles. The extension therefore ships disabled, with no margin frozen, and the audit's published results are exactly what they were before it existed.

This document records that experiment in full: the pre-registered rule, the grid, what the extension would have bought, and the counterexample. It is kept because a measured negative result is worth more than an unmeasured feature, and because the next person to have this idea should be able to see what it costs.

## Decision

The extension is implemented and tested, but **no margin in the recorded grid
passed the noise gate**. It must remain experimental. No deployment margin is
frozen, and no new CT-confirmed annotation error is claimed. The requested
second-corpus validation is blocked by the public data actually available.
These are failed empirical acceptance gates, not successful validations.

All original Python modules, original tests and published Paris4 report/export
bytes are preserved. New evidence is in separate files. The original report's
17 reproducible sections matched repeated complete legacy runs.

## Margin selection

The grid is 0.5, 0.75, 1, 1.25, 1.5, 2 and 3 times the calibrated median wrap
spacing. The grid was recorded after exploratory baseline inspection and before
null-control evaluation. Twenty identical seed-777 Gaussian jitter draws,
sigma one voxel, were applied to relative and same-winding points; absolute
anchors and calibration stayed fixed, following the legacy protocol.

| Margin | Null rounds with a new alarm | Decision |
|---|---:|---|
| 0.5 | 20/20 | Reject |
| 0.75 | 20/20 | Reject |
| 1 | 18/20 | Reject |
| 1.25 | 20/20 | Reject |
| 1.5 | 20/20 | Reject |
| 2 | 20/20 | Reject |
| 3 | 18/20 | Reject |

The requested strict rule allows zero alarms. Even the old selection allowance
of two alarms in twenty would reject every candidate. `selected` and
`frozen_configuration` are therefore null. Each experimental configuration,
including its margin, nevertheless has a reproducible hash. Planted-based
ranking is not performed for disqualified candidates. These results do not
prove that every conceivable margin must fail.

An alarm follows the legacy convention: a newly quarantined collection or a
newly implicated endpoint of a conflicting equality pair, relative to the
unmodified baseline. The initial exploratory logs used quarantine alone;
`align_control_alarms.py` replayed the identical mutations to include both
conflicting endpoints before finalising the measurements. The final selection
and diagnostic scripts implement this convention directly.

## Before and after

The after column uses **diagnostic margin 1, not a selected or approved margin**.
Structural coverage does not imply that a one-wrap defect will be detected:
inequalities may have slack, and the underlying geometry can be wrong.

| Measurement | Legacy | Order extension, diagnostic m=1 |
|---|---:|---:|
| Collections in graph | 89 | 177 |
| Skipped-wrap boundaries | 7/1,823 (0.384%) | 142/1,823 (7.789%) |
| Sheet-switch boundaries | 48/5,288 (0.908%) | 3,963/5,288 (74.943%) |
| Alarms on the same 14 legacy skipped-wrap trials | 100% | 100% |
| Certain localisation on those 14 trials | 14.286% | 85.714% |
| Alarms on the same 96 legacy sheet-switch trials | 57.292% | 85.417% |
| Certain localisation on those 96 trials | 11.458% | 43.750% |

All legacy reachable sites were tested in both directions. Newly covered sites
were separately sampled: 32 sites per family, seed 20260918, both directions.
They were not used to rescue or tune the rejected margin.

| Newly covered sample | Trials | Alarm rate | Newly certain localisation |
|---|---:|---:|---:|
| Skipped wrap | 64 | 23.438% | 6.250% |
| Sheet switch | 64 | 40.625% | 25.000% |

All 238 diagnostic trials have resolved final results. Three initially hit a
60-second repair budget; after mathematically valid cycle cuts strengthened the
solver, those same cases were retried with the same budget and resolved.
`order_retries.json` preserves the original unknown outcomes and the retries.
A fresh run on a slower machine can still return unknown; rate bounds include
unresolved trials, rather than treating them as success or failure. This is a
sampled new-site evaluation, not an exhaustive evaluation of thousands of new
boundaries. Localisation on newly reached sites is much weaker than on the
old sites, despite the coverage increase.

The pair census also differs from the proposed informal count. Applying seam
safety and projection availability gives 54,401 comparable pairs within one
cell, 882 collection pairs and 168 nodes. Including neighbouring cells as the
legacy matcher does gives 139,308 pairs, 998 collection pairs and 186 nodes.
After support and margin gates the combined graph has 177 nodes. The proposed
54,548/893/187 count is not silently substituted for these measured quantities.

## Certificate and physical limitation

At m=1 the recorded negative cycle contains `same:158` and `same:159` and has
bound sum -1. The JSON lists directions, bounds, original point IDs and voxel
coordinates. Both collections have all-optima membership **possible**, not
certain. The same cycle relations survive all 27 geometric sweep settings.
The exact quarantine optimum is recorded in `order_report.json`; the exact
minimum unit-cost relation-group cut is 96 groups. Edge cut is decomposed over
strongly connected components, which is valid because its groups have fixed
endpoints. This decomposition is not applied to quarantine's coupled secondary
objective. Both formulations and canonical answers are checked against small
exhaustive references.

Neither robustness nor optimisation certainty establishes physical truth.
Consider the smooth polar spiral

`r(t) = 1000 + 20*t/(2*pi) + 100*sin(20*t)`.

Along any fixed ray, successive turns differ by exactly 20 and do not intersect.
Yet between t=1.2 and t=1.25, points on the same turn differ in radius by more
than 20. Both angles fit inside the four-degree tolerance. Their z coordinates
may be identical, so z projection cannot repair the error. The automated test
constructs this counterexample and verifies that the proposed inequality
rejects the true equal gauges. Non-intersection alone therefore does not
justify the claimed cross-angle implication. A useful next model needs
independent angular transport or an explicit, independently justified bound on
local deformation before its order inequalities can be called reliable.

## Second corpus

The public index listed Paris4 and 25 other datasets. Every non-Paris root,
scan and tracks directory index was inspected. None of those listed locations
contained all four required winding/umbilicus JSON inputs. The others exposed
track data, which was not downloaded or relabelled as winding annotations.
`results/second_corpus/discovery.json` contains the URLs, directory contents and
index hashes. No second-corpus run or generalisation claim is fabricated.

## Reproduction and scope

`bash scripts/reproduce.sh` reruns the complete legacy pipeline, compares its
17 reproducible sections, checks both exported files byte-for-byte, and rebuilds
the additive graph to validate its recorded coverage, certificates and control
records. It does **not** repeat all expensive order controls. For those, use
`scripts/reproduce_order.sh`. Full repeat measurement of the 140 order null
runs and 238 diagnostic trials was not performed as a second independent batch.
The repeated legacy runs and deterministic graph reconstruction are separately
recorded in `results/tests.txt`.

The genuine upstream loader function at commit
`2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7` parses all 2,371 exported links in the
new test. Its original AST is executed directly; no parser is reimplemented.
This is not a full import/execution of Villa's GPU and patch-linking dependency
stack. The source, licence and hash are bundled.

Unit-cost evidence groups are used only by the additive order repair, preserving
the old weighted results and original membership implementation. Old SciPy
compatibility and the independent mathematical checks are documented in the
test log.

## Borderline legacy equality links

`results/paris4/borderline_links.json` audits the margin-3 certificate's
same:149, same:150 and same:151 using the unchanged 0.2.1 equality matcher,
without order constraints. Coordinates, every local baseline witness, all
27 settings and separate support checks are recorded. The one-hop neighbour
union adds same:148. All three target traces have constant label zero and
z = 12418.1259765625; no estimated z slope is available for them.

| Equality edge | Support at 0.35D | Support at 0.25D | Present in 27 settings |
| --- | ---: | ---: | ---: |
| 148–149 | 44 | 27 | 27 |
| 148–150 | 0 | 0 | 9 |
| 149–150 | 25 | 23 | 27 |
| 149–151 | 12 | 8 | 24 |
| 150–151 | 45 | 31 | 27 |

The 0.28187D certificate witness is not the sole support for 149–151:
the smallest witness gap is 0.00591D. Radial tightening alone to 0.25D,
even with minimum support five, preserves the edge. Joint tightening to
0.25D and three degrees leaves one witness and removes the edge at all
three z windows (30, 40, 60). At 0.35D and three degrees, four witnesses
remain, so minimum support five also removes it. The baseline witnesses
use seven points per trace and angular separations of 2.457–3.994 degrees.
This is a parameter-sensitive angular-alias matching candidate in the
legacy tool, independently of the new order mathematics; physical error
is not established. Persistence under radial tightening does not prove
that the fault must lie elsewhere in the cycle.

The legacy crowding veto makes zero adjacent-label tests in the target
neighbourhoods, both at baseline and across the sweep: its zero vetoes do
not establish absence of crowding. Descriptive same-angle interpolation
finds 150–151 gaps below the 0.7D threshold (16.307 voxels), with a minimum
of 9.591 voxels; 149–151 stays at least 32.071 voxels apart. Physical wrap
adjacency remains unverified. `results/paris4/borderline_links.png` shows
annotation geometry, not CT imagery, and the differing angular positions
of equality witnesses.

Across the corpus, 964 of 2,371 point-pair witnesses (40.66%) lie in the
0.25–0.35D band, affecting 76 of 98 equality edges; nine edges have only
band witnesses. Tightening to 0.25D leaves 1,398 witnesses and 80 edges.
The 973 lost witnesses comprise the 964 band witnesses plus nine inside
the tighter tolerance whose remaining edge support falls below two.
Thus 18 edges disappear. These are sensitivity measurements, not a count
of proven physical errors. Existing exports and all previous results are
unchanged; the new JSON records verification of this audit.

## Implementation notes

The extension is fully implemented and tested even though it is disabled; these notes describe what is in the code, so that the rejection can be checked rather than taken on trust.

The additive `windaudit.order_constraints` module keeps the 0.2.1 audit,
exports, configuration hash and all published sections unchanged. Its own
configuration includes the radial margin and has a separate hash. The legacy
package version stays at 0.2.1 to preserve its byte-for-byte exports; the order
extension identifies itself as `0.3.0-experimental`.

An arc `u -> v` with integer bound `c` asserts `g_v - g_u <= c`. A same-wrap
relation contributes two opposite arcs. If the projected radius of point `a`
is greater than that of `b` by more than `m * D`, an order relation contributes
`g_B - g_A <= label_diff(a, b) - 1`. `label_diff` includes branch-cut transport.
Bellman-Ford returns a directed cycle whose bounds sum to a negative integer,
with the contributing collection IDs, point IDs and voxel coordinates.

This certificate proves inconsistency of the **inferred constraints**. It does
not establish a CT-confirmed annotation error. Non-intersection alone does not
justify radial ordering between different angles. A smooth, non-intersecting
spiral can change radius by more than one wrap within a four-degree sector.
The explicit counterexample is tested in `tests/test_order_geometry.py`.
Projection to a common z does not remove this angular effect. Membership marked
`certain` means certain within the optimisation model, not certain physical
truth.

Repairs use one unit of cost per relation group, independent of the number of
supporting point pairs. The two arcs of an equality share a release variable.
Quarantine first minimises removed collections and then lost relation groups;
edge cut minimises lost groups. The additive implementation supports the same
big-M integer formulation and an exact cycle-separation formulation for larger
graphs. The edge-cut solver additionally decomposes cyclic strongly connected components; its exactness is tested against the monolithic solver. Both formulations reuse the unchanged lexicographic optimiser and its variable-range
queries over all optimal repairs. Independent Floyd-Warshall and exhaustive
small-graph tests check feasibility, both repair objectives and membership.
No legacy repair implementation or membership test is changed.

`results/selection_grid_order.json` records all seven margins, 20 Gaussian
one-voxel null rounds with seed 777, coverage and the selection decision. The
rule is stricter than 0.2.1: **any** incremental null alarm rejects a margin.
The grid was fixed after exploratory baseline inspection and before the null
and planted evaluations. Rejected margins cannot become defaults. If no margin
passes, `selected` and `frozen_configuration` are null, rather than selecting a
convenient value. Only passing margins qualify for exhaustive planted-based
selection.

To inspect an explicitly unvalidated reference margin:

```bash
python -m windaudit.order_cli --inputs data/paris4 --margin 1 \
  --diagnostic --sweep --out out/paris4/order_report.json
```

`results/paris4/order_report.json` contains the reference certificate, exact
quarantine and edge cut, membership for the review candidates and a 27-setting
sweep. The reported robustness is survival of the **same cycle relations**,
not merely occurrence of some contradiction. The diagnostic command requires
an explicit margin and acknowledgement; it does not modify certified links.

Coverage for the new graph requires a directed path independent of the target
collection, with the frozen support requirement on both boundary sides.
Inequality slack can still hide a one-wrap defect. Counts are structural
coverage, not sensitivity. `order_pair_census.json` distinguishes comparisons
inside a single cell from comparisons including adjacent cells.
`order_diagnostic.json` evaluates both defect directions at every legacy site
and at 32 seeded, newly covered sites per family. Its two cohorts are reported
separately; sampled detection rates must not be presented as an exhaustive
measurement of the extended corpus. Baseline quarantines are subtracted using
the existing incremental-alarm convention.

The genuine upstream `load_point_collection` function at the requested Villa
commit is executed directly from the pinned source AST in an offline test.
This isolates the file loader from optional GPU/patch imports without replacing
the parser. It is a loader compatibility test, not a full Villa integration or
surface-linking run. Source, licence and hash are under `upstream/villa/`.

The second-corpus availability check is recorded under `results/second_corpus/`.
A missing complete four-file corpus is a blocked validation, never a successful
transfer result.

For the expensive order-control reruns:

```bash
python scripts/select_order.py out/selection_grid_order.json
python scripts/diagnostic_order.py
```

The second command regenerates `results/paris4/order_diagnostic.json`.

The diagnostic repair budget is 60 seconds per trial; unresolved outcomes are
reported as unknown with rate bounds. Three initial unknowns were resolved by
retries after valid cycle-cut strengthening; the retry ledger is retained.
