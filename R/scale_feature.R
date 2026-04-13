#' Feature Scale Specifications
#'
#' Constructor functions for defining how feature values map to fill colours
#' in pedigree segments.
#'
#' @name scale_feature
NULL

#' Continuous Feature Scale
#'
#' Maps numeric feature values to a colour gradient.
#'
#' @param low Colour for the minimum value (default: "white")
#' @param high Colour for the maximum value (default: "red")
#' @param mid Optional colour for the midpoint (creates a diverging scale)
#' @param midpoint Numeric midpoint value when \code{mid} is specified
#' @param limits Numeric vector of length 2 giving (min, max). If NULL,
#'   computed from data.
#' @param na.value Colour for NA values (default: "grey50")
#' @param name Display name for the legend. If NULL, uses the feature column name.
#' @return An S3 object of class \code{scale_feature_continuous}
#' @export
#' @examples
#' scale_feature_continuous(low = "white", high = "red")
#' scale_feature_continuous(low = "blue", high = "red", mid = "white", midpoint = 0)
scale_feature_continuous <- function(low = "white",
                                     high = "red",
                                     mid = NULL,
                                     midpoint = NULL,
                                     limits = NULL,
                                     na.value = "grey50",
                                     name = NULL) {
  structure(
    list(
      type = "continuous",
      low = low,
      high = high,
      mid = mid,
      midpoint = midpoint,
      limits = limits,
      na.value = na.value,
      name = name
    ),
    class = "scale_feature_continuous"
  )
}

#' Discrete Feature Scale
#'
#' Maps categorical feature values to specific colours.
#'
#' @param values Named character vector mapping factor levels to colours,
#'   e.g. \code{c("wt/wt" = "white", "wt/ko" = "orange", "ko/ko" = "red")}.
#' @param na.value Colour for NA values (default: "grey50")
#' @param name Display name for the legend. If NULL, uses the feature column name.
#' @return An S3 object of class \code{scale_feature_discrete}
#' @export
#' @examples
#' scale_feature_discrete(values = c("wt/wt" = "white", "wt/ko" = "orange", "ko/ko" = "red"))
scale_feature_discrete <- function(values,
                                    na.value = "grey50",
                                    name = NULL) {
  stopifnot(is.character(values), !is.null(names(values)))
  structure(
    list(
      type = "discrete",
      values = values,
      na.value = na.value,
      name = name
    ),
    class = "scale_feature_discrete"
  )
}

#' Resolve feature values to hex colours
#'
#' @param scale A \code{scale_feature_continuous} or \code{scale_feature_discrete} object
#' @param x Vector of feature values
#' @return Character vector of hex colour strings
#' @export
resolve_colours <- function(scale, x) {
  UseMethod("resolve_colours")
}

#' @export
resolve_colours.scale_feature_continuous <- function(scale, x) {
  x_num <- as.numeric(x)
  lims <- scale$limits
  if (is.null(lims)) {
    lims <- range(x_num, na.rm = TRUE)
  }
  if (lims[1] == lims[2]) {
    lims <- lims + c(-0.5, 0.5)
  }

  # Normalise to [0, 1]
  x_norm <- (x_num - lims[1]) / (lims[2] - lims[1])
  x_norm <- pmin(pmax(x_norm, 0), 1)

  # Build colour ramp
  if (!is.null(scale$mid)) {
    mp <- if (!is.null(scale$midpoint)) {
      (scale$midpoint - lims[1]) / (lims[2] - lims[1])
    } else {
      0.5
    }
    mp <- pmin(pmax(mp, 0), 1)
    ramp <- grDevices::colorRampPalette(c(scale$low, scale$mid, scale$high))
    pal <- ramp(256)
  } else {
    ramp <- grDevices::colorRampPalette(c(scale$low, scale$high))
    pal <- ramp(256)
  }

  idx <- pmax(1L, pmin(256L, as.integer(x_norm * 255) + 1L))
  out <- pal[idx]
  out[is.na(x_num)] <- scale$na.value
  out
}

#' @export
resolve_colours.scale_feature_discrete <- function(scale, x) {
  x_chr <- as.character(x)
  out <- scale$values[x_chr]
  out[is.na(out)] <- scale$na.value
  unname(out)
}
