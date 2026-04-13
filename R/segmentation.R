#' Segmentation Engine for Pedigree Point Symbols
#'
#' Provides precomputed optimal rotation angles and functions to generate
#' segment polygons for pedigree symbols (circle, square, diamond) split
#' into n equal-area regions.
#'
#' @name segmentation
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Precomputed optimal theta lookup table
# ---------------------------------------------------------------------------
# Each entry: list(theta, vertex_hits)
# Computed by minimizing area variance of n ray-partitioned segments,
# with tie-breaking by max vertex hits then smallest theta.

# Precomputed optimal theta lookup table
# All shapes use the SAME theta for cross-shape visual consistency.
# Chosen to avoid overlap with the 45-degree deceased slash mark:
#   n=2: θ=π/2 (vertical split → left, right)
#   n=3: θ=π/2 (rays at 90°, 210°, 330° → top, bottom-left, bottom-right)
#   n=4: θ=0   (cross → right, top, left, bottom quadrants)

.THETA_LOOKUP <- list(
  circle  = list(`1` = 0,        `2` = pi / 2,   `3` = pi / 2,     `4` = 0),
  square  = list(`1` = 0,        `2` = pi / 2,   `3` = pi / 2,     `4` = 0),
  diamond = list(`1` = 0,        `2` = pi / 2,   `3` = pi / 2,     `4` = 0)
)

# Sex → shape mapping (matching original ggped)
.SEX_SHAPE_MAP <- c("1" = "square", "M" = "square", "male" = "square",
                     "2" = "circle", "F" = "circle", "female" = "circle",
                     "0" = "diamond", "U" = "diamond", "unknown" = "diamond")

#' Get optimal theta for a shape and number of segments
#' @param shape One of "circle", "square", "diamond"
#' @param n Integer 1-4
#' @return Numeric angle in radians
#' @keywords internal
get_optimal_theta <- function(shape, n) {
  stopifnot(shape %in% c("circle", "square", "diamond"))
  stopifnot(n >= 1 && n <= 4)
  .THETA_LOOKUP[[shape]][[as.character(n)]]
}

# ---------------------------------------------------------------------------
# Shape vertex definitions
# ---------------------------------------------------------------------------

.square_vertices <- function() {
  # CCW from top-right
  matrix(c(1, 1, -1, 1, -1, -1, 1, -1), ncol = 2, byrow = TRUE)
}

.diamond_vertices <- function() {
  # CCW from right

  matrix(c(1, 0, 0, 1, -1, 0, 0, -1), ncol = 2, byrow = TRUE)
}

# ---------------------------------------------------------------------------
# Ray-polygon intersection
# ---------------------------------------------------------------------------

#' Find intersection of ray from origin in direction d with convex polygon edge
#' @keywords internal
.ray_polygon_intersect <- function(origin, direction, vertices) {
  n <- nrow(vertices)
  best_t <- Inf
  best_pt <- NULL

  for (i in seq_len(n)) {
    j <- ifelse(i == n, 1L, i + 1L)
    p1 <- vertices[i, ]
    p2 <- vertices[j, ]
    edge <- p2 - p1
    d <- direction

    denom <- d[1] * (-edge[2]) - d[2] * (-edge[1])
    if (abs(denom) < 1e-14) next

    diff <- p1 - origin
    t_val <- (diff[1] * (-edge[2]) - diff[2] * (-edge[1])) / denom
    s_val <- (d[1] * diff[2] - d[2] * diff[1]) / denom

    if (t_val > 1e-10 && s_val >= -1e-10 && s_val <= 1 + 1e-10) {
      if (t_val < best_t) {
        best_t <- t_val
        best_pt <- origin + t_val * d
      }
    }
  }
  best_pt
}

# ---------------------------------------------------------------------------
# Segment polygon generation
# ---------------------------------------------------------------------------

#' Generate segment polygons for a pedigree symbol
#'
#' @param shape One of "circle", "square", "diamond"
#' @param n_segments Integer 1-4
#' @param scale Numeric size scaling factor (default 1)
#' @param n_circle_pts Number of points for circle arc approximation
#' @return List of n_segments matrices, each with columns (x, y) defining
#'   a closed polygon for that segment.
#' @export
get_segment_polygons <- function(shape, n_segments, scale = 1, n_circle_pts = 80) {
  stopifnot(shape %in% c("circle", "square", "diamond"))
  stopifnot(n_segments >= 1 && n_segments <= 4)

  theta <- get_optimal_theta(shape, n_segments)

  if (shape == "circle") {
    .get_circle_segments(n_segments, theta, scale, n_circle_pts)
  } else {
    verts <- if (shape == "square") .square_vertices() else .diamond_vertices()
    verts <- verts * scale
    .get_polygon_segments(verts, n_segments, theta)
  }
}

.get_circle_segments <- function(n, theta, radius, n_pts) {
  if (n == 1) {
    angles <- seq(0, 2 * pi, length.out = n_pts + 1)
    pts <- cbind(radius * cos(angles), radius * sin(angles))
    return(list(pts))
  }

  center <- c(0, 0)
  segments <- vector("list", n)
  for (k in seq_len(n)) {
    a_start <- theta + 2 * pi * (k - 1) / n
    a_end   <- theta + 2 * pi * k / n
    arc_angles <- seq(a_start, a_end, length.out = n_pts)
    arc_pts <- cbind(radius * cos(arc_angles), radius * sin(arc_angles))
    # Close: center → arc → center
    segments[[k]] <- rbind(center, arc_pts, center)
  }
  segments
}

.get_polygon_segments <- function(vertices, n, theta) {
  if (n == 1) {
    return(list(rbind(vertices, vertices[1, , drop = FALSE])))
  }

  centroid <- colMeans(vertices)
  n_verts <- nrow(vertices)

  # Compute ray intersection points
  ray_data <- lapply(seq_len(n), function(k) {
    angle <- theta + 2 * pi * (k - 1) / n
    d <- c(cos(angle), sin(angle))
    pt <- .ray_polygon_intersect(centroid, d, vertices)
    # Find which edge this point is on
    edge_idx <- .find_edge(pt, vertices)
    list(pt = pt, edge = edge_idx, angle = angle %% (2 * pi))
  })

  # Sort by angle
  angles <- sapply(ray_data, function(x) x$angle)
  ord <- order(angles)
  ray_data <- ray_data[ord]

  # Build segment polygons
  segments <- vector("list", n)
  for (i in seq_len(n)) {
    j <- ifelse(i == n, 1L, i + 1L)
    pt_start  <- ray_data[[i]]$pt
    edge_start <- ray_data[[i]]$edge
    pt_end    <- ray_data[[j]]$pt
    edge_end  <- ray_data[[j]]$edge

    seg_pts <- rbind(centroid, pt_start)

    # Walk boundary CCW from edge_start to edge_end
    e <- edge_start
    safety <- 0
    while (safety < n_verts + 2) {
      next_v <- ifelse(e == n_verts, 1L, e + 1L)
      if (e == edge_end) {
        seg_pts <- rbind(seg_pts, pt_end)
        break
      } else {
        seg_pts <- rbind(seg_pts, vertices[next_v, ])
        e <- next_v
      }
      safety <- safety + 1
    }
    # Close back to centroid
    seg_pts <- rbind(seg_pts, centroid)
    segments[[i]] <- seg_pts
  }
  segments
}

#' Determine which edge of the polygon a point lies on
#' @keywords internal
.find_edge <- function(pt, vertices) {
  n <- nrow(vertices)
  for (i in seq_len(n)) {
    j <- ifelse(i == n, 1L, i + 1L)
    p1 <- vertices[i, ]
    p2 <- vertices[j, ]
    # Check collinearity and between-ness
    cross <- (pt[1] - p1[1]) * (p2[2] - p1[2]) - (pt[2] - p1[2]) * (p2[1] - p1[1])
    if (abs(cross) < 1e-8) {
      # Check if pt is between p1 and p2
      dot <- sum((pt - p1) * (p2 - p1))
      len2 <- sum((p2 - p1)^2)
      if (dot >= -1e-8 && dot <= len2 + 1e-8) {
        return(i)
      }
    }
  }
  1L  # fallback
}

#' Map sex factor level to shape name
#' @param sex_level Character or factor level
#' @return Shape name: "circle", "square", or "diamond"
#' @keywords internal
sex_to_shape <- function(sex_level) {
  sl <- as.character(sex_level)
  if (sl %in% names(.SEX_SHAPE_MAP)) {
    return(.SEX_SHAPE_MAP[[sl]])
  }
  # ggped encodes: 1=male→22 (square pch), 2=female→1 (circle pch), else→0 (diamond pch)
  if (sl == "22") return("square")
  if (sl == "0")  return("diamond")
  "circle"  # default female
}
