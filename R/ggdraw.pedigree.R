#' Plot Pedigree Using ggplot2 with Multi-Feature Segments
#'
#' Generates a pedigree plot using ggplot2, where each individual's symbol
#' is divided into n segments (one per feature), each independently filled
#' according to a continuous or discrete colour scale.
#'
#' @param dat A \code{kinship2::pedigree} object or data frame from
#'   \code{\link{dfalign.pedigree}}.
#' @param features Named list of feature scale specifications. Names must
#'   correspond to columns in \code{dat}. Values are objects from
#'   \code{\link{scale_feature_continuous}} or
#'   \code{\link{scale_feature_discrete}}. Maximum 4 features.
#' @param plot.names Logical. Plot individual names?
#' @param plot.kinship.label Logical. Show kinship coefficient for
#'   consanguineous matings?
#' @param column.names List mapping internal column names to actual column
#'   names in \code{dat}. Defaults match \code{dfalign.pedigree} output.
#' @param col.tree Colour for pedigree tree lines.
#' @param col.double Colour for lines connecting repeated subjects.
#' @param chrtype Chromosome type: "autosome" or "X".
#' @param packed Logical. Compress layout?
#' @param align Logical. Align children under parents?
#' @param width Minimum width for packed layout.
#' @param allow.repeated Allow individuals plotted multiple times?
#' @param size Symbol size in grid points (default 7).
#' @param debug Show relationship IDs, drawing occurrences, and anchors.
#' @param label.mode Use legacy offsets (default) or measured automatic placement.
#' @param label.control Named list: size and min.size in points; family, fontface;
#'   gap in inches; wrap in characters; stagger (number of extra rows);
#'   shrink (FALSE by default). Unresolved constraints are reported, and all
#'   requested labels remain in the drawing.
#' @param label.overrides Data frame with label_id and optional x, y, dx, dy.
#'   x/y are drawing coordinates (negative generations); dx/dy are inch offsets.
#'   Label IDs are available in plot$diagnostics$labels after drawing.
#' @param auto.size Opt in to automatic export sizing. Interactive printing always
#'   uses the current graphics device dimensions.
#' @param size.control Named list: min.width, min.height, max.width, max.height,
#'   x.spacing, y.spacing, margin (inches), and max.iter (1 to 10).
#' @param ... Further arguments passed to \code{\link{geom_pedigreepoint}}.
#' @return A ggped_plot containing a ggplot in `plot`, relationship and connection
#'   tables, recommended_dimensions (inches), options, and a diagnostics
#'   environment updated on each print or save. The environment contains labels,
#'   collisions, device size, and export dimensions when applicable.
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom cli cli_abort cli_warn
#' @export
#' @examples
#' \dontrun{
#' library(kinship2)
#' data(minnbreast)
#' bpeds <- with(minnbreast,
#'   pedigree(id, fatherid, motherid, sex,
#'            affected = proband, famid = famid))
#' bped.id8 <- bpeds['8']
#' df <- dfalign.pedigree(bped.id8)
#'
#' # Single binary feature (as discrete)
#' ggdraw.pedigree(df, features = list(
#'   affected = scale_feature_discrete(
#'     values = c("TRUE" = "black", "FALSE" = "white"))
#' ))
#' }
ggdraw.pedigree <- function(dat = NULL,
                            features = list(),
                            plot.names = TRUE,
                            plot.kinship.label = TRUE,
                            column.names = list(
                              x = "x", y = "y", Name = "Name",
                              mate.id = "mate.id", twin.id = "twin.id",
                              mate.centerpoint = "mate.centerpoint",
                              family.centerpoint = "family.centerpoint",
                              twin.centerpoint = "twin.centerpoint",
                              sex = "sex", adopted = "adopted",
                              pregnancy = "pregnancy"
                            ),
                            col.tree = "#000000",
                            col.double = "#808080",
                            chrtype = "autosome",
                            packed = TRUE,
                            align = TRUE,
                            width = 10,
                            allow.repeated = FALSE,
                            size = 7,
                            debug = FALSE,
                            label.mode = c("legacy", "auto"),
                            label.control = list(),
                            label.overrides = NULL,
                            auto.size = FALSE,
                            size.control = list(),
                            ...) {
  # ---- Input validation ----
  if (is.null(dat)) cli::cli_abort("Input data 'dat' cannot be NULL.")
  if (length(features) == 0) cli::cli_abort("'features' must contain at least one feature specification.")
  if (length(features) > 4) cli::cli_abort("Maximum 4 features supported (got {length(features)}).")
  if (is.null(names(features)) || anyNA(names(features)) || any(!nzchar(names(features))) || anyDuplicated(names(features)))
    cli::cli_abort("'features' must have unique, nonempty names.")
  if (!all(vapply(features, function(s) inherits(s, c("scale_feature_continuous", "scale_feature_discrete")), logical(1))))
    cli::cli_abort("Every feature must use scale_feature_continuous() or scale_feature_discrete().")
  .check_positive(size, "size")
  .check_flag(plot.names, "plot.names"); .check_flag(plot.kinship.label, "plot.kinship.label")
  .check_flag(allow.repeated, "allow.repeated"); .check_flag(debug, "debug")
  label.mode <- match.arg(label.mode)
  label.control <- .label_options(label.control)
  .check_flag(auto.size, "auto.size")
  size.control <- .size_options(size.control)

  # Convert pedigree object if needed
  if (inherits(dat, "pedigree")) {
    dat <- dfalign.pedigree(ped = dat, chrtype = chrtype,
                            packed = packed, align = align, width = width)
  } else if (!is.data.frame(dat)) {
    cli::cli_abort("'dat' must be a data.frame or kinship2 pedigree object.")
  }

  if (nrow(dat) < 2) cli::cli_abort("Pedigree must have at least two individuals.")

  # Check feature columns exist
  for (fn in names(features)) {
    if (!(fn %in% colnames(dat))) {
      cli::cli_abort("Feature column '{fn}' not found in data. Available: {paste(colnames(dat), collapse=', ')}")
    }
  }

  # ---- Prepare fill colours ----
  n_feat <- length(features)
  feat_names <- names(features)

  # Ensure unique individual identifier
  if (!"ID" %in% colnames(dat)) {
    dat$ID <- seq_len(nrow(dat))
  }
  dat$individual_id <- as.character(dat$ID)

  # Resolve colours for each feature × individual
  colour_matrix <- matrix(NA_character_, nrow = nrow(dat), ncol = n_feat)
  colnames(colour_matrix) <- feat_names

  for (j in seq_len(n_feat)) {
    fn <- feat_names[j]
    colour_matrix[, j] <- resolve_colours(features[[fn]], dat[[fn]])
  }

  # Build segment_fills lookup: individual_id → c(col1, col2, ...)
  segment_fills <- setNames(
    lapply(seq_len(nrow(dat)), function(i) colour_matrix[i, ]),
    dat$individual_id
  )

  # ---- Column renaming ----
  for (nm in names(column.names)) {
    target <- column.names[[nm]]
    if (target %in% colnames(dat) && nm != target) {
      colnames(dat)[colnames(dat) == target] <- nm
    }
    if (!(nm %in% colnames(dat)) && nm %in% c("adopted", "pregnancy")) {
      dat[[nm]] <- FALSE
    }
  }
  for (nm in c("adopted", "pregnancy")) if (!nm %in% names(dat)) dat[[nm]] <- FALSE
  if (!"status" %in% names(dat)) dat$status <- 0
  if (anyNA(dat$adopted) || anyNA(dat$pregnancy)) cli::cli_abort("'adopted' and 'pregnancy' cannot contain missing values.")

  model <- .relationship_model(dat)
  dat <- model$occurrences
  if (anyDuplicated(dat$person_id) && !allow.repeated) {
    cli::cli_abort("This layout repeats individuals. Set 'allow.repeated = TRUE'.")
  }
  dat$individual_id <- dat$occurrence_id
  names(segment_fills) <- dat$individual_id
  segments <- .connection_segments(model)
  segments$y <- -segments$y
  segments$yend <- -segments$yend
  segments$colour <- ifelse(segments$kind == "repeated", col.double, col.tree)
  dat$y <- -dat$y
  shape.size <- size
  text.size <- 8
  .pt <- 72.27 / 25.4
  plt <- ggplot(dat, aes(x = x, y = y)) +
    geom_segment(data = segments,
      aes(x = x, y = y, xend = xend, yend = yend,
          group = relationship_id, linetype = I(linetype), colour = I(colour)),
      inherit.aes = FALSE)
  if (plot.kinship.label && label.mode == "legacy") {
    k <- model$partnerships[model$partnerships$kinship > 0, , drop = FALSE]
    if (nrow(k)) plt <- plt + geom_text(data = k,
      aes(x = x, y = -y + 0.08, label = paste0("Kinship: ", kinship)),
      inherit.aes = FALSE, size = text.size / .pt, vjust = 0)
  }
  # ---- Individual symbols ----
  plt <- plt + geom_pedigreepoint(
    mapping = aes(
      sex           = as.factor(sex),
      isdead        = as.factor(status),
      individual_id = individual_id,
      adopted       = adopted,
      pregnancy     = pregnancy
    ),
    segment_fills = segment_fills,
    n_features    = n_feat,
    size          = shape.size,
    ...
  )

  # ---- Name labels ----
  voffset <- 0
  if ("Name" %in% colnames(dat) && plot.names && label.mode == "legacy") {
    plt <- plt + geom_text(
      data = dat,
      aes(label = Name),
      vjust = -(shape.size * 1.5 / min(dat$y)),
      hjust = "outward",
      size = text.size / .pt
    )
    voffset <- -(shape.size * 1.5 / min(dat$y)) + text.size / .pt * 1.01
  }

  # ---- Formatting ----
  plt <- plt +
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(expand = expansion(add = 0.5))
  if (debug) {
    plt <- plt + geom_text(aes(label = occurrence_id), colour = "#005a9c",
      nudge_y = 0.14, size = 2) +
      geom_point(data = model$partnerships, aes(x = x, y = -y),
        inherit.aes = FALSE, colour = "red", shape = 4) +
      geom_text(data = model$partnerships,
        aes(x = x, y = -y + 0.25, label = relationship_id),
        inherit.aes = FALSE, colour = "red", size = 2)
  }

  # ---- Build result object with embedded legend specs ----
  result <- structure(
    list(
      plot     = plt,
      features = features,
      dat      = dat,
      relationships = model,
      connections = segments,
      options = list(size = size, plot.names = plot.names, plot.kinship.label = plot.kinship.label,
        label.mode = label.mode, label.control = label.control, label.overrides = label.overrides,
        auto.size = auto.size, size.control = size.control),
      diagnostics = new.env(parent = emptyenv())
    ),
    class = "ggped_plot"
  )
  .validate_overrides(label.overrides, .label_requests(result))
  result$recommended_dimensions <- .recommended_dimensions(result)
  result
}

# =========================================================================
# Legend building helpers
# =========================================================================

#' Build all legend grobs for a ggped plot
#' @keywords internal
#' @noRd
.build_all_legend_grobs <- function(features, dat, max_width = Inf) {
  all_legends <- list()

  # Shape legend
  shape_df <- data.frame(
    x     = c(0, 0, 0),
    y     = c(0, 0, 0),
    shape = factor(c("Male", "Female", "Unknown"),
                    levels = c("Male", "Female", "Unknown"))
  )
  shape_rows <- if (max_width < 3) 3 else 1
  p_shape <- ggplot(shape_df, aes(x = x, y = y, shape = shape)) +
    geom_point(size = 3, fill = "white", colour = "black") +
    scale_shape_manual(
      values = c(Male = 22, Female = 21, Unknown = 23),
      name   = "Sex", guide = guide_legend(nrow = shape_rows, title.position = "top")
    ) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9, face = "bold"),
          legend.text  = element_text(size = 8))
  g <- .extract_legend_safe(p_shape)
  if (!is.null(g)) all_legends[["shape:Sex"]] <- .pad_legend(g)

  # Feature legends
  feat_names <- names(features)
  for (fn in feat_names) {
    scale <- features[[fn]]
    display_name <- if (!is.null(scale$name)) scale$name else fn
    if (is.finite(max_width)) display_name <- .wrap_label(display_name, max(4, floor(max_width * 11)))

    if (inherits(scale, "scale_feature_continuous")) {
      vals <- as.numeric(dat[[fn]])
      lims <- scale$limits
      if (is.null(lims)) lims <- .feature_limits(vals)
      if (lims[1] == lims[2]) lims <- lims + c(-0.5, 0.5)
      colours <- if (!is.null(scale$mid)) {
        c(scale$low, scale$mid, scale$high)
      } else {
        c(scale$low, scale$high)
      }
      dummy_df <- data.frame(x = 0, y = 0, val = mean(lims, na.rm = TRUE))
      p_leg <- ggplot(dummy_df, aes(x = x, y = y, fill = val)) +
        geom_tile() +
        scale_fill_gradientn(
          colours = colours, limits = lims, name = display_name,
          guide = guide_colorbar(
            barwidth = grid::unit(min(1.2, max_width * 0.7), "inches"),
            barheight = grid::unit(0.4, "cm"),
            title.position = "top"
          )
        ) +
        theme_void() +
        theme(legend.position = "bottom",
              legend.title = element_text(size = 9, face = "bold"),
              legend.text  = element_text(size = 8))

    } else {
      lvls <- names(scale$values)
      dummy_df <- data.frame(
        x = seq_along(lvls), y = 1,
        val = factor(lvls, levels = lvls)
      )
      p_leg <- ggplot(dummy_df, aes(x = x, y = y, fill = val)) +
        geom_tile(colour = "black", linewidth = 0.3) +
        scale_fill_manual(
          values = scale$values, name = display_name,
          labels = vapply(lvls, .wrap_label, character(1), width = if (is.finite(max_width)) max(4, floor(max_width * 9)) else 40),
          guide = guide_legend(title.position = "top", ncol = 1)
        ) +
        theme_void() +
        theme(legend.position = "bottom",
              legend.title = element_text(size = 9, face = "bold"),
              legend.text  = element_text(size = 8))
    }

    g <- .extract_legend_safe(p_leg)
    if (!is.null(g)) all_legends[[paste0("feature:", fn)]] <- .pad_legend(g)
  }

  all_legends
}

#' Safely extract legend grob from a ggplot
#' @keywords internal
#' @noRd
.extract_legend_safe <- function(p) {
    g <- ggplot2::ggplotGrob(p)
    idx <- which(grepl("^guide-box", g$layout$name) & !vapply(g$grobs, inherits, logical(1), "zeroGrob"))
    if (length(idx) == 0) return(NULL)
    g$grobs[[idx[1]]]
}

.pad_legend <- function(g) {
  # Colourbar endpoint tick labels can extend beyond the bar's gtable width.
  # Reserve their physical overhang before packing adjacent legends.
  g <- gtable::gtable_add_cols(g, grid::unit(0.16, "inches"), pos = 0)
  gtable::gtable_add_cols(g, grid::unit(0.16, "inches"), pos = -1)
}

# =========================================================================
# Drawing function (used by both print and save)
# =========================================================================

#' Draw a ggped_plot onto the current graphics device
#' @param x A ggped_plot object
#' @keywords internal
#' @noRd
.draw_ggped <- function(x) {
  device <- grDevices::dev.size("in")
  prepared <- .prepare_canvas(x, device[1], device[2])
  x$diagnostics$labels <- prepared$labels
  x$diagnostics$collisions <- prepared$collisions
  x$diagnostics$device <- device
  x$diagnostics$legend_overflow <- prepared$overflow
  issues <- character()
  if (nrow(prepared$collisions)) issues <- c(issues, paste0(nrow(prepared$collisions), " label(s) have unresolved placement constraints. Inspect plot$diagnostics$collisions."))
  if (prepared$overflow) issues <- c(issues, "Legends exceed the available canvas; increase width or height.")
  if (length(issues)) cli::cli_warn(paste(issues, collapse = " "))
  main <- grid::viewport(layout = grid::grid.layout(4, 3,
    widths = grid::unit(c(prepared$margin, prepared$plot_width, prepared$margin), "inches"),
    heights = grid::unit(c(prepared$margin, prepared$plot_height, prepared$legends$height, prepared$margin), "inches")))
  grid::pushViewport(main)
  on.exit(grid::popViewport())
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid::grid.draw(prepared$plot)
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  grid::grid.draw(prepared$legends$grob)
  grid::popViewport()
  invisible(x)
}
# =========================================================================
# Print and save methods
# =========================================================================

#' Print method for ggped plots
#' @param x A ggped_plot object
#' @param ... Ignored
#' @export
print.ggped_plot <- function(x, ...) {
  grid::grid.newpage()
  .draw_ggped(x)
}

#' Plot method for ggped plots (RStudio compatibility)
#' @param x A ggped_plot object
#' @param ... Ignored
#' @export
plot.ggped_plot <- function(x, ...) {
  print.ggped_plot(x, ...)
}

#' Save a ggped plot to file
#'
#' Supports png, pdf, svg, tiff, and jpeg output based on file extension.
#'
#' @param plot A ggped_plot object (or a regular ggplot)
#' @param filename Output file path
#' @param width Width in inches (default 10)
#' @param height Height in inches (default 8)
#' @param dpi Resolution for raster formats (default 150)
#' @param auto.size Override the plot's automatic sizing option. NULL inherits it.
#'   Explicit width and height always take precedence. Omitted dimensions use
#'   automatic recommendations when enabled, or 10 by 8 inches otherwise.
#' @return Invisibly, the output filename. See plot$diagnostics for placement
#'   constraints and the dimensions actually used.
#' @export
save_ggped <- function(plot, filename, width = NULL, height = NULL, dpi = 150,
                       auto.size = NULL) {
  .check_positive(dpi, "dpi")
  if (!is.null(width)) .check_positive(width, "width")
  if (!is.null(height)) .check_positive(height, "height")
  if (!inherits(plot, "ggped_plot")) {
    ggplot2::ggsave(filename, plot, width = if (is.null(width)) 10 else width,
      height = if (is.null(height)) 8 else height, dpi = dpi)
    return(invisible(filename))
  }
  if (is.null(auto.size)) auto.size <- plot$options$auto.size
  .check_flag(auto.size, "auto.size")
  ext <- tolower(tools::file_ext(filename))
  if (ext == "") { filename <- paste0(filename, ".png"); ext <- "png" }
  # Measure on the requested device type, including raster-specific font metrics.
  probe <- tempfile(fileext = paste0(".", ext))
  .open_ggped_device(probe, ext, 10, 8, dpi)
  probe_device <- grDevices::dev.cur()
  dims <- tryCatch(.choose_dimensions(plot, width, height, auto.size),
    finally = { grDevices::dev.off(probe_device); unlink(probe) })
  plot$diagnostics$export_dimensions <- dims[c("width", "height")]
  plot$diagnostics$sizing_iterations <- unname(dims["iterations"])
  .open_ggped_device(filename, ext, dims["width"], dims["height"], dpi)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device))
  grid::grid.newpage()
  .draw_ggped(plot)
  invisible(filename)
}
