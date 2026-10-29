# All collision calculations use inches on the current graphics device.
.label_options <- function(control) {
  defaults <- list(size = 8, family = "", fontface = "plain", gap = 0.04,
    wrap = 20L, stagger = 4L, shrink = FALSE, min.size = 6)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || any(!names(control) %in% names(defaults)))))
    cli::cli_abort("'label.control' must be a named list of supported label controls.")
  opts <- utils::modifyList(defaults, control)
  for (nm in c("size", "gap", "wrap", "stagger", "min.size")) .check_positive(opts[[nm]], paste0("label.control$", nm))
  .check_flag(opts$shrink, "label.control$shrink")
  if (opts$min.size > opts$size) cli::cli_abort("'min.size' cannot exceed label size.")
  opts
}

.text_bounds <- function(label, size = 8, family = "", fontface = "plain") {
  g <- grid::textGrob(label, gp = grid::gpar(fontsize = size, fontfamily = family, fontface = fontface, lineheight = 1.2))
  c(width = grid::convertWidth(grid::grobWidth(g), "inches", valueOnly = TRUE),
    height = grid::convertHeight(grid::grobHeight(g), "inches", valueOnly = TRUE))
}
.text_bounds_many <- function(labels, size = 8, family = "", fontface = "plain") {
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = size,
    fontfamily = family, fontface = fontface, lineheight = 1.2)))
  on.exit(grid::popViewport())
  cbind(width = grid::convertWidth(grid::stringWidth(labels), "inches", valueOnly = TRUE),
    height = grid::convertHeight(grid::stringHeight(labels), "inches", valueOnly = TRUE))
}
.wrap_label <- function(x, width) {
  paste(vapply(strsplit(x, "\n", fixed = TRUE)[[1]], function(line) {
    paste(strwrap(line, width = width), collapse = "\n")
  }, character(1)), collapse = "\n")
}
.box <- function(x, y, w, h) c(left = x - w / 2, right = x + w / 2, bottom = y - h / 2, top = y + h / 2)
.overlap_boxes <- function(box, boxes, pad = 0) {
  if (!nrow(boxes)) return(logical())
  box[1] < boxes[, 2] + pad & box[2] > boxes[, 1] - pad & box[3] < boxes[, 4] + pad & box[4] > boxes[, 3] - pad
}
.line_hits_box <- function(box, s, pad = 0.015) {
  if (!nrow(s)) return(logical())
  box <- box + c(-pad, pad, -pad, pad)
  lo <- rep(0, nrow(s)); hi <- rep(1, nrow(s)); possible <- rep(TRUE, nrow(s))
  for (axis in c("x", "y")) {
    start <- s[[axis]]; delta <- s[[paste0(axis, "end")]] - start
    lower <- box[if (axis == "x") 1 else 3]; upper <- box[if (axis == "x") 2 else 4]
    flat <- abs(delta) < 1e-12
    possible[flat & (start < lower | start > upper)] <- FALSE
    delta[flat] <- 1
    a <- (lower - start) / delta; b <- (upper - start) / delta
    lo[!flat] <- pmax(lo[!flat], pmin(a, b)[!flat]); hi[!flat] <- pmin(hi[!flat], pmax(a, b)[!flat])
  }
  possible & lo <= hi
}

.label_requests <- function(x) {
  d <- x$dat
  labels <- data.frame(label_id = character(), occurrence_id = character(), person_id = character(),
    kind = character(), text = character(), x = numeric(), y = numeric())
  if (x$options$plot.names) labels <- data.frame(label_id = paste0("name:", d$occurrence_id),
    occurrence_id = d$occurrence_id, person_id = d$person_id, kind = "name",
    text = ifelse(is.na(d$Name), "NA", as.character(d$Name)), x = d$x, y = d$y)
  k <- x$relationships$partnerships
  k <- k[k$kinship > 0, , drop = FALSE]
  if (x$options$plot.kinship.label && nrow(k)) labels <- rbind(labels, data.frame(
    label_id = paste0("kinship:", k$relationship_id), occurrence_id = NA_character_, person_id = NA_character_,
    kind = "kinship", text = paste0("Kinship: ", format(k$kinship)), x = k$x, y = -k$y))
  labels
}

.validate_overrides <- function(overrides, requests) {
  if (is.null(overrides)) return(NULL)
  if (!is.data.frame(overrides) || !"label_id" %in% names(overrides))
    cli::cli_abort("'label.overrides' must be a data frame with 'label_id' and optional x, y, dx, dy (dx/dy in inches).")
  if (anyNA(overrides$label_id) || anyDuplicated(overrides$label_id) || !all(overrides$label_id %in% requests$label_id))
    cli::cli_abort("Label override IDs must be unique and match requested label IDs.")
  for (nm in intersect(c("x", "y", "dx", "dy"), names(overrides))) {
    if (!is.numeric(overrides[[nm]]) || any(!is.na(overrides[[nm]]) & !is.finite(overrides[[nm]])))
      cli::cli_abort("Label override coordinates must be finite numeric values or NA.")
  }
  overrides
}

.place_labels <- function(x, panel_width, panel_height, ranges) {
  opts <- x$options$label.control
  req <- .label_requests(x); overrides <- .validate_overrides(x$options$label.overrides, req)
  xr <- ranges$x.range; yr <- ranges$y.range
  sx <- function(v) (v - xr[1]) / diff(xr) * panel_width
  sy <- function(v) (v - yr[1]) / diff(yr) * panel_height
  invx <- function(v) v / panel_width * diff(xr) + xr[1]
  invy <- function(v) v / panel_height * diff(yr) + yr[1]
  d <- x$dat; px <- sx(d$x); py <- sy(d$y)
  half <- rep(x$options$size / 72.27, nrow(d))
  half[as.character(d$status) %in% c("1", "TRUE")] <- half[as.character(d$status) %in% c("1", "TRUE")] * 1.3
  hw <- hh <- half
  if (any(d$adopted)) hw[d$adopted] <- pmax(hw[d$adopted], .text_bounds("[  ]", 18)[1] / 2)
  symbols <- cbind(px - hw, px + hw, py - hh, py + hh)
  lines <- x$connections
  lines$x <- sx(lines$x); lines$xend <- sx(lines$xend); lines$y <- sy(lines$y); lines$yend <- sy(lines$yend)
  used <- matrix(numeric(), ncol = 4); placed <- list()
  req <- req[order(req$kind != "name", -req$y, req$x, req$label_id), , drop = FALSE]
  generations <- sort(unique(d$y))
  for (i in seq_len(nrow(req))) {
    r <- req[i, ]; text <- .wrap_label(r$text, opts$wrap)
    sizes <- if (opts$shrink) unique(c(seq(opts$size, opts$min.size, by = -0.5), opts$min.size)) else opts$size
    manual <- if (is.null(overrides)) NULL else overrides[overrides$label_id == r$label_id, , drop = FALSE]
    best <- NULL; best_score <- Inf
    for (fs in sizes) {
      wh <- .text_bounds(text, fs, opts$family, opts$fontface)
      anchor_x <- sx(r$x); anchor_y <- sy(r$y)
      direction <- if (r$kind == "name") -1 else 1
      own <- match(r$occurrence_id, d$occurrence_id)
      extent <- if (is.na(own)) x$options$size / 72.27 else hh[own]
      base_y <- anchor_y + direction * (extent + opts$gap + wh[2] / 2)
      gen <- match(r$y, generations)
      lower <- if (gen > 1) sy(mean(generations[c(gen - 1, gen)])) else 0
      upper <- if (gen < length(generations)) sy(mean(generations[c(gen, gen + 1)])) else panel_height
      candidates <- expand.grid(stagger = 0:as.integer(opts$stagger), shift = c(0, -0.25, 0.25, -0.5, 0.5))
      if (!is.null(manual) && nrow(manual)) candidates <- candidates[1, , drop = FALSE]
      for (j in seq_len(nrow(candidates))) {
        cx <- anchor_x + candidates$shift[j] * wh[1]
        cy <- base_y + direction * candidates$stagger[j] * (wh[2] + opts$gap)
        # A small edge shift is allowed, while preserving proximity to the owner.
        if (wh[1] < panel_width) cx <- max(wh[1] / 2, min(panel_width - wh[1] / 2, cx))
        if (!is.null(manual) && nrow(manual)) {
          if ("x" %in% names(manual) && !is.na(manual$x)) cx <- sx(manual$x)
          if ("y" %in% names(manual) && !is.na(manual$y)) cy <- sy(manual$y)
          if ("dx" %in% names(manual) && !is.na(manual$dx)) cx <- cx + manual$dx
          if ("dy" %in% names(manual) && !is.na(manual$dy)) cy <- cy + manual$dy
        }
        box <- .box(cx, cy, wh[1], wh[2])
        hits <- c(symbol = sum(.overlap_boxes(box, symbols, opts$gap / 2)),
          label = sum(.overlap_boxes(box, used, opts$gap / 2)), line = sum(.line_hits_box(box, lines)),
          canvas = as.integer(box[1] < -1e-8 || box[2] > panel_width + 1e-8 || box[3] < 0 || box[4] > panel_height),
          band = as.integer(box[3] < lower || box[4] > upper))
        score <- sum(hits * c(1000, 100, 10, 10000, 2))
        if (score < best_score) {
          best_score <- score; best <- list(x = cx, y = cy, box = box, size = fs, hits = hits)
        }
        if (score == 0) break
      }
      if (best_score == 0 || (!is.null(manual) && nrow(manual))) break
    }
    used <- rbind(used, best$box)
    placed[[i]] <- data.frame(label_id = r$label_id, occurrence_id = r$occurrence_id, kind = r$kind,
      text = text, x = invx(best$x), y = invy(best$y), size = best$size,
      left = best$box[1], right = best$box[2], bottom = best$box[3], top = best$box[4],
      symbol = unname(best$hits[1]), label = unname(best$hits[2]), line = unname(best$hits[3]),
      canvas = unname(best$hits[4]), band = unname(best$hits[5]), manual = !is.null(manual) && nrow(manual) > 0)
  }
  if (!length(placed)) return(data.frame(label_id = character(), text = character(), x = numeric(), y = numeric(), size = numeric(), symbol = integer(), label = integer(), line = integer(), canvas = integer(), band = integer()))
  do.call(rbind, placed)
}

.label_diagnostics <- function(labels) {
  labels[rowSums(labels[, c("symbol", "label", "line", "canvas", "band"), drop = FALSE]) > 0, , drop = FALSE]
}

.panel_size <- function(gt, width, height) {
  panel <- gt$layout[gt$layout$name == "panel", , drop = FALSE]
  if (nrow(panel) != 1L) cli::cli_abort("Automatic labels require a single pedigree panel.")
  outside_w <- setdiff(seq_along(gt$widths), panel$l:panel$r)
  outside_h <- setdiff(seq_along(gt$heights), panel$t:panel$b)
  fixed_w <- if (length(outside_w)) grid::convertWidth(sum(gt$widths[outside_w]), "inches", valueOnly = TRUE) else 0
  fixed_h <- if (length(outside_h)) grid::convertHeight(sum(gt$heights[outside_h]), "inches", valueOnly = TRUE) else 0
  c(width = max(0.05, width - fixed_w), height = max(0.05, height - fixed_h))
}

.plot_with_labels <- function(x, width, height) {
  if (x$options$label.mode != "auto") return(list(grob = ggplot2::ggplotGrob(x$plot), labels = NULL))
  built <- ggplot2::ggplot_build(x$plot)
  gt <- ggplot2::ggplot_gtable(built)
  panel <- .panel_size(gt, width, height)
  labels <- .place_labels(x, panel[1], panel[2], built$layout$panel_params[[1]])
  p <- x$plot
  # Freeze the measured ranges so adding labels cannot change their conversion.
  original_coord <- p$coordinates
  p$coordinates <- ggplot2::ggproto(NULL, original_coord,
    limits = list(x = built$layout$panel_params[[1]]$x.range,
                  y = built$layout$panel_params[[1]]$y.range), expand = FALSE)
  if (nrow(labels)) p <- p + ggplot2::geom_text(data = labels,
    ggplot2::aes(x = x, y = y, label = text), inherit.aes = FALSE,
    size = labels$size / (72.27 / 25.4), family = x$options$label.control$family,
    fontface = x$options$label.control$fontface, lineheight = 1.2)
  list(grob = ggplot2::ggplotGrob(p), labels = labels, panel = panel)
}
