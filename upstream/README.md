# Upstream: `find_inconsistent_windings.py`

This directory holds everything needed to check three defects in ScrollPrize/villa's winding diagnostic, and a candidate patch for two of them.

Base: ScrollPrize/villa `2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7`. The file `spiral-fitting/find_inconsistent_windings.py` on `main` as of 20 September 2026 (`23c3d759`) is byte-identical to the pinned one (SHA-256 `2fe5ee4f5b71981d…`), so everything below applies to current `main`. `villa/` holds the pinned source files used by the harnesses, with their licence and provenance.

## The defects

**1. Non-editable edges lose their equations.** `solve_min_edge_fix` is documented as restricting which edges may *change*. The code instead skips every edge outside `allowed_edge_keys` when building the model, so the equations that pin the other potentials disappear, and a repair can close one inconsistent cycle by breaking a consistent one that was never in the model. On the real PHercParis4 patch graph (`results/wide/`, `results/wide_corrected/`) the model keeps 7 of 645 equations, or 2 of 645 once the attachment gap is transported.

**2. The potential box can exclude the optimum.** Potentials are boxed around a BFS labelling, `big = max|u0 - mean(u0)| + m + 10`. On a graph with large offsets that box is too small to contain the true optimum, and the solver reports an optimum of the wrong problem.

**3. The attachment gap is never transported.** `build_rel_adjacency` computes each edge's branch correction along the annotation at the annotation points' own coordinates (`_tour_unwrap_adjustments` on `p['zyx']`), while every strip starts and ends at the point's attachment on the surface (`on_patch['ij']`). The segment between a point and its attachment is counted by neither. When the branch ray passes through it, the edge is off by one winding and the diagnostic reports a contradiction that is not there. On the real patch graph this happens at `col109` point 1560; run as it is, the diagnostic then recommends changing `col106` 1530 → 1529 from −1 to −2, an edit that breaks two consistent cycles once the gap is counted. See `WIDE_SCANSPACE_AUDIT.md`.

## The candidate patch

`solver-candidate.patch` fixes defects 1 and 2 in `solve_min_edge_fix`:

* every measurable reached edge stays in the model; `allowed_edge_keys` now sets only whether an edge's edit indicator may be non-zero;
* potentials are bounded by the sum of absolute edge offsets plus one, a bound that contains an optimum for every admissible retained graph, with one root fixed per model;
* the MILP is asked for a zero relative gap, so a returned integer optimum is an optimum.

It applies cleanly to both the pinned revision and current `main`. Defect 3 is not in the patch: in upstream's code the fix needs the attachment's position in spiral space at the point where each strip ends, and it is implemented and tested here in the patch-graph harness instead (`WIDE_ATTACHMENT_GAP=1` in `scripts/wide_patch_audit.py`), where it adds the crossing count from each strip's last valid sample to its annotation point.

## Which patch to look at

`solver-candidate.patch` is the patch proposed for upstream review. It changes only `solve_min_edge_fix` in `find_inconsistent_windings.py`. `scanspace-combined.patch` contains the same solver change plus three new standalone files (`scanspace_transport.py`, `scanspace_surface.py`, `audit_scanspace_transport.py`) that run the checkpoint-free scan-space audit inside a villa checkout; it is an optional integration artefact, recorded in `scanspace_provenance.json`, not a replacement for the candidate patch.

## Checking it

```bash
# in a clean checkout of villa, at the pinned commit or at current main
python /path/to/windaudit/upstream/check_solver.py spiral-fitting/find_inconsistent_windings.py
# exit 1: both regressions fail
git apply /path/to/windaudit/upstream/solver-candidate.patch
python /path/to/windaudit/upstream/check_solver.py spiral-fitting/find_inconsistent_windings.py
# exit 0: both regressions pass
```

`check_solver.py` loads the named function from source, so it needs neither a GPU nor villa's patch-loading stack. Its two regressions are small synthetic graphs, one for each solver defect; `before.json` and `after.json` are its recorded outputs. The real-data evidence for defects 1 and 3 comes from the full patch-graph runs in `results/`, reproduced by `scripts/reproduce_wide.sh`.

The patch touches a diagnostic, not `fit_spiral`'s optimisation, so applying it does not change any fit.
