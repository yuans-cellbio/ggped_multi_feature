# ggped v2 — Multi-Feature Segmented Pedigree Symbols

> **This is a fork of [moritzlindner/ggped](https://github.com/moritzlindner/ggped)**
> by [Moritz Lindner](http://lindnerlab.de), extended with multi-feature
> segmented pedigree symbol support.

A package to draw pedigree charts using ggplot2 and kinship2, with support for
**n ≤ 4 features** per individual, each mapped to an independently coloured
segment of the pedigree symbol.

## What this fork adds

- **Multi-feature segments**: Each symbol (square/circle/diamond) is partitioned
  into n equal-area segments using an optimal ray-from-centroid algorithm.
- **Continuous and discrete colour scales**: Per-feature colour mapping via
  `scale_feature_continuous()` and `scale_feature_discrete()`.
- **Per-feature legends**: Each feature gets its own legend entry (colour bar
  for continuous, key-value for discrete).
- **Visible segment boundaries**: Black outlines around each segment.
- **Consistent segmentation**: Optimal rotation angle θ is precomputed for all
  shape × n combinations and hardcoded, ensuring identical partitioning across
  all individuals of the same sex.
- **Bug fix**: Corrected de-duplication logic that dropped valid drawing positions
  in multi-mate pedigree layouts.

See [NEWS.md](NEWS.md) for the full changelog.

## Shapes

| Sex       | Shape   | Symbol |
|-----------|---------|--------|
| Male      | Square  | □      |
| Female    | Circle  | ○      |
| Undefined | Diamond | ◇      |

## Precomputed Segmentation

| Shape   | n=1 | n=2 (θ)   | n=3 (θ)  | n=4 (θ)   |
|---------|-----|-----------|----------|-----------|
| Circle  | 0°  | 0°        | 0°       | 0°        |
| Square  | 0°  | 45°       | 15°      | 45°       |
| Diamond | 0°  | 0°        | 0°       | 0°        |

## Installation

```r
# Install from this fork:
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("[YOUR GITHUB USERNAME]/ggped")

# Or install the original upstream package:
# remotes::install_github("moritzlindner/ggped")
```

## Quick Start

```r
library(ggped)
library(kinship2)

# Build pedigree
ped_df <- data.frame(
  id       = c("F1", "M1", "S1", "D1"),
  fatherid = c(NA, NA, "F1", "F1"),
  motherid = c(NA, NA, "M1", "M1"),
  sex      = c(1, 2, 1, 2)
)
ped <- with(ped_df, pedigree(id, dadid = fatherid, momid = motherid, sex = sex))
df <- dfalign.pedigree(ped)

# Add feature data
df$disease_score <- c(0.2, 0.1, 0.85, 0.45)
df$genotype      <- c("wt/ko", "wt/wt", "ko/ko", "wt/ko")

# Plot with 2 features
p <- ggdraw.pedigree(df, features = list(
  disease_score = scale_feature_continuous(
    low = "white", high = "red", limits = c(0, 1), name = "Disease Score"
  ),
  genotype = scale_feature_discrete(
    values = c("wt/wt" = "white", "wt/ko" = "orange", "ko/ko" = "red"),
    name = "Genotype"
  )
))

# Save
save_ggped(p, "my_pedigree.png", width = 10, height = 8)
```

## API Reference

### `ggdraw.pedigree(dat, features, ...)`
Main plotting function. `features` is a named list where names are column names
in `dat` and values are `scale_feature_*` objects.

### `scale_feature_continuous(low, high, mid, midpoint, limits, na.value, name)`
Maps numeric values to a colour gradient. Supports diverging scales via `mid`/`midpoint`.

### `scale_feature_discrete(values, na.value, name)`
Maps categorical values to specific colours. `values` is a named character vector.

### `dfalign.pedigree(ped, ...)`
Converts a kinship2 pedigree object to a ggplot2-compatible data frame.

### `get_segment_polygons(shape, n_segments, scale)`
Returns polygon vertices for each segment of a shape. Useful for custom rendering.

### `save_ggped(plot, filename, width, height, dpi)`
Saves a ggped plot (handles custom legend grobs correctly).

## Dependencies

- R ≥ 4.0.0
- ggplot2, grid, kinship2, tidyr, cli, scales, gtable

## License

GPL-2. See [LICENSE](LICENSE) for details.

## Acknowledgements

This fork is based on the original [ggped](https://github.com/moritzlindner/ggped)
package by [Moritz Lindner](http://lindnerlab.de). The original package provides
the core pedigree layout engine (`dfalign.pedigree`) and ggplot2 integration.
