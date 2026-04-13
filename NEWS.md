# ggped 2.0.0

This is a fork of [moritzlindner/ggped](https://github.com/moritzlindner/ggped)
with major API changes to support multi-feature pedigree visualization.

## New features

* **Multi-feature segmented symbols**: Each pedigree symbol (square/circle/diamond)
  can be partitioned into up to 4 equal-area segments using an optimal
  ray-from-centroid algorithm. Each segment maps independently to a feature.
* **`scale_feature_continuous()`**: Map numeric feature values to colour gradients,
  with optional diverging scales via `mid`/`midpoint`.
* **`scale_feature_discrete()`**: Map categorical feature values to specific colours
  via a named character vector.
* **`geom_pedigreepoint()`**: New ggplot2 geom that draws segmented pedigree symbols
  with per-segment fill colours and visible segment boundaries.
* **`get_segment_polygons()`**: Utility to retrieve precomputed polygon vertices for
  each segment of a given shape and segment count.
* **`save_ggped()`**: Export helper that correctly handles the custom legend grobs
  attached to ggped plots.
* **Per-feature legends**: Each feature receives its own legend entry (colour bar
  for continuous, key-value for discrete).

## Breaking changes

* The old binary `features` argument (character vector) in `ggdraw.pedigree()` is
  replaced by a named list of `scale_feature_*` objects. Existing code using the
  v1 API will need to be updated.

## Bug fixes

* Fixed de-duplication logic in `geom_pedigreepoint()` that incorrectly removed
  duplicate drawing positions for individuals appearing multiple times in
  multi-mate pedigree layouts.

## Internal

* Precomputed optimal rotation angles for all shape x n combinations (n = 1..4)
  are hardcoded in `segmentation.R` for consistency across individuals.
* Legend construction refactored into internal helpers in `ggdraw.pedigree.R`.
