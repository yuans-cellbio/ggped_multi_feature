#' geom_pedigreepoint
#'
#' The pedigreepoint geom draws pedigree symbols (square for male, circle for
#' female, diamond for undefined) split into n segments, each independently
#' filled according to a feature's colour scale.
#'
#' @inheritParams ggplot2::geom_point
#' @name geom_pedigreepoint
#' @section Aesthetics:
#' \code{geom_pedigreepoint()} understands the following aesthetics:
#' \itemize{
#'  \item \strong{x} Numeric x coordinate.
#'  \item \strong{y} Numeric y coordinate.
#'  \item \strong{sex} Factor. 1/M = Male (square), 2/F = Female (circle), 0/U = Unknown (diamond).
#'  \item \strong{isdead} Factor or Logical. 0/FALSE = Alive, 1/TRUE = Dead.
#'  \item \strong{individual_id} Unique individual identifier.
#'  \item adopted Logical. Brackets around symbol.
#'  \item pregnancy Logical. P in symbol.
#'  \item colour Stroke colour. Default: "black".
#'  \item alpha Symbol alpha. Default: 1.
#'  \item stroke Stroke thickness.
#' }
#' @importFrom ggplot2 ggproto Geom aes
#' @importFrom grid unit polygonGrob polylineGrob segmentsGrob pointsGrob
#'   textGrob gpar gList unit.c
#' @export
NULL

GeomPedigreePoint <- ggproto(
"GeomPedigreePoint",
ggplot2::Geom,
required_aes = c("x", "y", "sex", "isdead"),
optional_aes = c(
  "individual_id",
  "adopted",
  "pregnancy",
  "colour",
  "alpha",
  "stroke"
),
default_aes = aes(
  adopted   = FALSE,
  pregnancy = FALSE,
  colour    = "black",
  alpha     = 1,
  stroke    = 1,
  isdead    = 0
),

# ---- Extra params carried through from the layer ----
# segment_fills: named list, keyed by individual_id, each a character vector
#   of n hex colours (one per segment, in segment order)
# n_features: integer, number of features (= number of segments)

setup_data = function(data, params) {
  if (!is.factor(data$sex)) {
    data$sex <- as.factor(data$sex)
  }
  if (!is.factor(data$isdead) && !is.logical(data$isdead)) {
    data$isdead <- as.factor(data$isdead)
  }
  data
},

draw_panel = function(data, panel_scales, coord,
                      segment_fills = list(),
                      n_features    = 1L,
                      size          = 5,
                      na.colour     = "#808080") {

  coords <- coord$transform(data, panel_scales)

  n_seg <- as.integer(n_features)
  grob_list <- grid::gList()

  for (i in seq_len(nrow(coords))) {
    row <- coords[i, ]
    cx <- grid::unit(row$x, "npc")
    cy <- grid::unit(row$y, "npc")

    shape <- sex_to_shape(row$sex)
    seg_polys <- get_segment_polygons(shape, n_seg, scale = 1)

    # Retrieve fill colours for this individual
    ind_id <- if ("individual_id" %in% names(row)) as.character(row$individual_id) else paste(row$x, row$y)
    fills <- segment_fills[[ind_id]]
    if (is.null(fills)) {
      fills <- rep("#FFFFFF", n_seg)
    }

    # Draw each segment as a polygon
    for (s in seq_len(n_seg)) {
      poly_pts <- seg_polys[[s]]
      n_pts <- nrow(poly_pts)

      px <- cx + grid::unit(poly_pts[, 1] * size, "points")
      py <- cy + grid::unit(poly_pts[, 2] * size, "points")

      seg_grob <- grid::polygonGrob(
        x  = px,
        y  = py,
        gp = grid::gpar(
          fill = fills[s],
          col  = row$colour,
          alpha = row$alpha,
          lwd  = row$stroke
        )
      )
      grob_list <- grid::gList(grob_list, seg_grob)
    }

    # Overlay: deceased
    if (!is.null(row$isdead) && as.character(row$isdead) %in% c("1", "TRUE")) {
      # Diagonal slash through symbol
      slash_size <- size * 1.3
      grob_list <- grid::gList(grob_list,
        grid::segmentsGrob(
          x0 = cx - grid::unit(slash_size, "points"),
          y0 = cy - grid::unit(slash_size, "points"),
          x1 = cx + grid::unit(slash_size, "points"),
          y1 = cy + grid::unit(slash_size, "points"),
          gp = grid::gpar(col = row$colour, lwd = row$stroke * 1.5)
        )
      )
    }

    # Overlay: adopted (brackets)
    if (!is.null(row$adopted) && isTRUE(as.logical(row$adopted))) {
      grob_list <- grid::gList(grob_list,
        grid::textGrob(
          "[  ]",
          x  = cx,
          y  = cy,
          vjust = 0.4,
          gp = grid::gpar(
            cex   = 1.5,
            alpha = row$alpha,
            col   = row$colour
          )
        )
      )
    }

    # Overlay: pregnancy
    if (!is.null(row$pregnancy) && isTRUE(as.logical(row$pregnancy))) {
      grob_list <- grid::gList(grob_list,
        grid::textGrob(
          "P",
          x  = cx,
          y  = cy,
          vjust = 0.5,
          gp = grid::gpar(
            cex   = 0.8,
            alpha = row$alpha,
            col   = row$colour
          )
        )
      )
    }
  }

  grob_list
},

draw_key = function(data, params, size) {
  grid::pointsGrob(
    0.5, 0.5,
    pch = 22,
    gp  = grid::gpar(
      fill  = "#CCCCCC",
      col   = data$colour,
      alpha = data$alpha
    )
  )
}
)

#' @describeIn geom_pedigreepoint Create a pedigree point layer
#' @param segment_fills Named list of fill colour vectors per individual.
#' @param n_features Integer, number of features (segments).
#' @param na.colour Colour for NA features.
#' @param size Symbol size in grid points.
#' @importFrom ggplot2 layer
#' @export
geom_pedigreepoint <- function(data = NULL,
                               mapping = NULL,
                               stat = "identity",
                               position = "identity",
                               ...,
                               na.rm = FALSE,
                               na.colour = "#808080",
                               size = 5,
                               segment_fills = list(),
                               n_features = 1L,
                               show.legend = NA,
                               inherit.aes = TRUE) {
  ggplot2::layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomPedigreePoint,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      na.rm          = na.rm,
      na.colour      = na.colour,
      size           = size,
      segment_fills  = segment_fills,
      n_features     = n_features,
      ...
    )
  )
}
