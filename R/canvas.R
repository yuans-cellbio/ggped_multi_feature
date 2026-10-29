.size_options <- function(control) {
  defaults <- list(min.width = 4, min.height = 3, max.width = 40, max.height = 40,
    x.spacing = 0.75, y.spacing = 1, margin = 0.15, max.iter = 4L)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || any(!names(control) %in% names(defaults)))))
    cli::cli_abort("'size.control' must be a named list of supported sizing controls.")
  opts <- utils::modifyList(defaults, control)
  for (nm in names(opts)) .check_positive(opts[[nm]], paste0("size.control$", nm))
  if (opts$max.width < opts$min.width || opts$max.height < opts$min.height) cli::cli_abort("Maximum canvas dimensions cannot be smaller than minimum dimensions.")
  if (opts$max.iter > 10 || opts$max.iter != as.integer(opts$max.iter)) cli::cli_abort("'max.iter' must be an integer from 1 to 10.")
  opts
}

.measurement_device <- function(code) {
  opened <- grDevices::dev.cur() == 1L
  if (opened) grDevices::pdf(NULL, width = 10, height = 8)
  if (opened) on.exit(grDevices::dev.off())
  force(code)
}

.recommended_dimensions <- function(x) .measurement_device({
  opts <- x$options$size.control; lab <- x$options$label.control
  req <- .label_requests(x)
  bounds <- if (nrow(req)) .text_bounds_many(vapply(req$text, .wrap_label, character(1), width = lab$wrap), lab$size, lab$family, lab$fontface) else matrix(c(0, 0), nrow = 1)
  symbol <- 2 * x$options$size / 72.27
  dx <- max(opts$x.spacing, symbol + 0.15, max(bounds[, 1]) + 0.12)
  dy <- max(opts$y.spacing, symbol + max(bounds[, 2]) * 2 + 0.4)
  w <- max(opts$min.width, (diff(range(x$dat$x)) + 1) * dx + 2 * opts$margin)
  legends <- .legend_layout(x, w)
  h <- max(opts$min.height, (diff(range(x$dat$y)) + 1) * dy + legends$height + 2 * opts$margin)
  c(width = min(opts$max.width, w), height = min(opts$max.height, h))
})

.pack_legends <- function(grobs, width, gap = 0.12) {
  if (!length(grobs)) return(list(grob = grid::nullGrob(), width = 0, height = 0, overflow = FALSE))
  ws <- vapply(grobs, function(g) grid::convertWidth(sum(g$widths), "inches", valueOnly = TRUE), numeric(1))
  hs <- vapply(grobs, function(g) grid::convertHeight(sum(g$heights), "inches", valueOnly = TRUE), numeric(1))
  rows <- list(); row <- integer(); occupied <- 0
  for (i in seq_along(grobs)) {
    if (length(row) && occupied + gap + ws[i] > width) {
      rows[[length(rows) + 1L]] <- row; row <- integer(); occupied <- 0
    }
    occupied <- occupied + if (length(row)) gap + ws[i] else ws[i]; row <- c(row, i)
  }
  rows[[length(rows) + 1L]] <- row
  heights <- vapply(rows, function(r) max(hs[r]) + gap, numeric(1))
  table <- gtable::gtable(widths = grid::unit(width, "inches"), heights = grid::unit(heights, "inches"))
  for (j in seq_along(rows)) {
    r <- rows[[j]]; widths <- as.vector(rbind(ws[r], rep(gap, length(r))))
    widths <- head(widths, -1); spare <- max(0, (width - sum(widths)) / 2)
    row_table <- gtable::gtable(widths = grid::unit(c(spare, widths, spare), "inches"), heights = grid::unit(max(hs[r]), "inches"))
    for (k in seq_along(r)) row_table <- gtable::gtable_add_grob(row_table, grobs[[r[k]]], t = 1, l = 2 * k, clip = "off")
    table <- gtable::gtable_add_grob(table, row_table, t = j, l = 1, clip = "off")
  }
  list(grob = table, width = max(ws), height = sum(heights), overflow = any(ws > width + 1e-6))
}

.legend_layout <- function(x, width) {
  margin <- x$options$size.control$margin
  available <- max(0.25, width - 2 * margin)
  .pack_legends(.build_all_legend_grobs(x$features, x$dat, max_width = min(4, available)), available)
}

.prepare_canvas <- function(x, width, height) {
  legends <- .legend_layout(x, width)
  margin <- x$options$size.control$margin
  plot_height <- max(0.05, height - legends$height - 2 * margin)
  plot_width <- max(0.05, width - 2 * margin)
  prepared <- .plot_with_labels(x, plot_width, plot_height)
  labels <- prepared$labels
  collisions <- if (is.null(labels)) data.frame() else .label_diagnostics(labels)
  list(plot = prepared$grob, labels = labels, collisions = collisions, legends = legends,
    plot_height = plot_height, plot_width = plot_width, margin = margin,
    overflow = legends$overflow || legends$height + 2 * margin >= height)
}

.choose_dimensions <- function(x, width = NULL, height = NULL, auto.size = FALSE) .measurement_device({
  fixed_w <- !is.null(width); fixed_h <- !is.null(height)
  if (fixed_w) .check_positive(width, "width")
  if (fixed_h) .check_positive(height, "height")
  recommended <- if (auto.size && (!fixed_w || !fixed_h)) .recommended_dimensions(x) else c(10, 8)
  if (!fixed_w) width <- if (auto.size) recommended[1] else 10
  if (!fixed_h) height <- if (auto.size) recommended[2] else 8
  opts <- x$options$size.control
  iterations <- 0L
  if (auto.size) for (i in seq_len(opts$max.iter)) {
    iterations <- i
    prepared <- .prepare_canvas(x, width, height)
    if (!nrow(prepared$collisions) && !prepared$overflow) break
    next_w <- if (fixed_w) width else min(opts$max.width, width * 1.2)
    next_h <- if (fixed_h) height else min(opts$max.height, height * 1.2)
    if (identical(next_w, width) && identical(next_h, height)) break
    width <- next_w; height <- next_h
  }
  c(width = unname(width), height = unname(height), iterations = iterations)
})

.open_ggped_device <- function(filename, ext, width, height, dpi) {
  if (ext == "pdf") grDevices::pdf(filename, width = width, height = height)
  else if (ext == "svg") grDevices::svg(filename, width = width, height = height)
  else if (ext %in% c("png", "tif", "tiff", "jpg", "jpeg")) {
    device <- switch(ext, png = grDevices::png, tif = grDevices::tiff, tiff = grDevices::tiff,
      jpg = grDevices::jpeg, jpeg = grDevices::jpeg)
    device(filename, width = width, height = height, units = "in", res = dpi)
  } else cli::cli_abort("Unsupported export extension '{ext}'. Use PNG, PDF, SVG, TIFF, or JPEG.")
}
