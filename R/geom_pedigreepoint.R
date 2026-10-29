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
  if (!nrow(coords)) return(grid::nullGrob())
  n_seg <- as.integer(n_features)
  shapes <- vapply(as.character(coords$sex), sex_to_shape, character(1))
  templates <- setNames(lapply(unique(shapes), get_segment_polygons,
    n_segments = n_seg, scale = 1), unique(shapes))
  polygons <- unlist(unname(templates[shapes]), recursive = FALSE)
  lengths <- vapply(polygons, nrow, integer(1))
  vertices <- do.call(rbind, polygons) * size
  owner <- rep(seq_len(nrow(coords)), each = n_seg)
  fills <- unlist(lapply(seq_len(nrow(coords)), function(i) {
    id <- if ("individual_id" %in% names(coords)) as.character(coords$individual_id[i]) else paste(coords$x[i], coords$y[i])
    cols <- segment_fills[[id]]
    if (is.null(cols)) rep("#FFFFFF", n_seg) else cols
  }), use.names = FALSE)
  if (length(fills) != length(lengths)) cli::cli_abort("Each segment_fills entry must contain n_features colours.")
  # id.lengths keeps one polygon and colour per individual/feature while grid
  # handles all vertices in one grob. Template reuse is local to this draw.
  grobs <- list(grid::polygonGrob(
    x = grid::unit(rep(coords$x[owner], lengths), "npc") + grid::unit(vertices[, 1], "points"),
    y = grid::unit(rep(coords$y[owner], lengths), "npc") + grid::unit(vertices[, 2], "points"),
    id.lengths = lengths,
    gp = grid::gpar(fill = fills, col = coords$colour[owner], alpha = coords$alpha[owner], lwd = coords$stroke[owner])))
  dead <- which(as.character(coords$isdead) %in% c("1", "TRUE"))
  if (length(dead)) grobs[[length(grobs) + 1L]] <- grid::segmentsGrob(
    x0 = grid::unit(coords$x[dead], "npc") - grid::unit(size * 1.3, "points"),
    y0 = grid::unit(coords$y[dead], "npc") - grid::unit(size * 1.3, "points"),
    x1 = grid::unit(coords$x[dead], "npc") + grid::unit(size * 1.3, "points"),
    y1 = grid::unit(coords$y[dead], "npc") + grid::unit(size * 1.3, "points"),
    gp = grid::gpar(col = coords$colour[dead], lwd = coords$stroke[dead] * 1.5))
  adopted <- which(as.logical(coords$adopted))
  if (length(adopted)) grobs[[length(grobs) + 1L]] <- grid::textGrob("[  ]",
    x = grid::unit(coords$x[adopted], "npc"), y = grid::unit(coords$y[adopted], "npc"), vjust = 0.4,
    gp = grid::gpar(cex = 1.5, alpha = coords$alpha[adopted], col = coords$colour[adopted]))
  pregnancy <- which(as.logical(coords$pregnancy))
  if (length(pregnancy)) grobs[[length(grobs) + 1L]] <- grid::textGrob("P",
    x = grid::unit(coords$x[pregnancy], "npc"), y = grid::unit(coords$y[pregnancy], "npc"), vjust = 0.5,
    gp = grid::gpar(cex = 0.8, alpha = coords$alpha[pregnancy], col = coords$colour[pregnancy]))
  grid::gTree(children = do.call(grid::gList, grobs))
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
  .check_positive(size, "size")
  if (length(n_features) != 1L || !is.numeric(n_features) || is.na(n_features) || !n_features %in% 1:4)
    cli::cli_abort("'n_features' must be an integer from 1 to 4.")
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
