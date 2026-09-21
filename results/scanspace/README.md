# Pilot run — not the headline results

This directory is the first, deliberately small exercise of the scan-space transport: a band one z-slice wide around z = 12418, 46 patches, 357 attached annotation points, 22 reached patches and 53 measurable equations. That graph is consistent, so it demonstrates operation rather than any contradiction.

The numbers quoted in the README's patch-level results come from the wide run instead:

* `results/wide/` — the wide graph (1,637 patches, 645 equations) built exactly as upstream builds it;
* `results/wide_corrected/` — the same graph with the attachment gap transported, which is the construction the frame-invariance proposition applies to;
* `results/frame_invariance/` — the invariance test on both.

`results/scanspace_initial/` holds the same pilot on its original 42-patch selection.
