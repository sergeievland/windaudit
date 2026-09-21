# Protocol for a downstream fit comparison

This repository makes no claim about the quality of a fitted surface. Measuring one honestly is expensive, and the design has to be fixed before any number exists, or the comparison can be tuned after the fact. This document fixes that design in advance, so that a later run either follows it or visibly does not.

## Splits and controls

Freeze one development ROI, one validation ROI and one untouched final test ROI before any tuning. Whole absolute-anchor collections or spatial blocks must be withheld; neighbouring points are not independent holdout data. Exclude final-test anchors from calibration, overlap construction, threshold selection, quarantine decisions, fitting and model selection — the annotation audit in this repository calibrates on all supplied absolute anchors, so it cannot be reused unchanged as a held-out evaluation. Inspect the retained patches and fibres for duplicated or nearby evidence that leaks the same test labels, and define exclusion distances from physical scale before seeing any test outcome.

The baseline is unmodified pinned villa with identical training inputs. The treatment applies only exclusions that were confirmed independently. The control applies random exclusions matched on removed valid surface area, point budget, z extent and approximate winding region — matching patch count alone is not enough — with the selection seed frozen and the selected identifiers recorded.

Use paired optimiser seeds 1, 2 and 3 with identical step budgets, ROI, settings and input-sampling policy, recording any change in effective input counts. Interleave baseline, treatment and control within each seed and rotate the order across seeds, rather than spending the budget on three baselines first and many untested ideas afterwards. Explore on development results only, freeze one treatment, then evaluate once. Three seeds are limited evidence, not a precise success probability.

## Outcomes

The primary outcome is held-out geometric and winding error, measured by an evaluator that has not consumed test labels, reporting voxel errors and whole-wrap errors separately with their denominators, the spatial distribution and per-seed paired differences. That evaluator still has to be written against the installed villa revision.

Secondary outcomes are the native satisfaction metrics, runtime, peak VRAM, retained surface area and the native ink metrics. Satisfaction on training inputs is a diagnostic, not a result: a reduced denominator can raise it without improving geometry. `total_fg_pixels` is a repository ink-coverage metric rather than a standalone score — it can rise with duplicate coverage or false-positive ink — so unique physical surface area, line and column coherence and fixed-scale visual examples belong alongside it.

Use the existing upstream `runners/run_single.py` for fit, render and ink scoring once its prerequisites are present: an ink-prediction volume, `vc_render_tifxyz`, flattening support and an ink-model environment. `uv sync` for the fitter alone does not establish those.

## Release gate

Publish the measured result or the null result, with exact commands, commit, configuration, input hashes, every seed, failure logs, checkpoints and matched cross-sections. Keep planted defects separate from native failures, and never tune on the final test set.

[`spiral_lab/`](../../spiral_lab/) holds the pinned setup, data fetch and run launcher this protocol would use. Its orchestration, override handling, failure logging and refusal to overwrite a run directory are covered by the test suite, and its overrides are accepted by the pinned upstream `Config` class; the installation itself and the CUDA fit have not been executed here.

The patch-level work in [SCANSPACE_TRANSPORT.md](../../SCANSPACE_TRANSPORT.md) took a different route to real-data evidence, one that needs no fit at all. This protocol remains the design for the question that route cannot answer: whether a corrected constraint set produces a measurably better fitted surface.
