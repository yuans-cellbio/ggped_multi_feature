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
#' @param ... Further arguments passed to \code{\link{geom_pedigreepoint}}.
#' @return A ggplot object.
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
                            ...) {
  # ---- Input validation ----
  if (is.null(dat)) cli::cli_abort("Input data 'dat' cannot be NULL.")
  if (length(features) == 0) cli::cli_abort("'features' must contain at least one feature specification.")
  if (length(features) > 4) cli::cli_abort("Maximum 4 features supported (got {length(features)}).")

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

  # Flip y for top-down rendering
  dat$y <- -dat$y

  # ---- Robust per-mating descent target computation ----
  # The original dfalign.pedigree sets family.centerpoint per INDIVIDUAL,
  # which breaks when a parent has multiple mates: the last family processed
  # overwrites the centerpoint, so descent lines point to the wrong litter.
  #
  # Fix: compute a descent target per MATING by finding the specific children
  # whose dad.id AND mom.id match the two parents in each mating pair.
  dat$descent.target.x <- dat$family.centerpoint  # start with original

  if ("dad.id" %in% colnames(dat) && "mom.id" %in% colnames(dat)) {
    mate_ids <- unique(dat$mate.id[!is.na(dat$mate.id) & !is.na(dat$mate.centerpoint)])
    for (mid in mate_ids) {
      mate_rows <- which(dat$mate.id == mid)
      if (length(mate_rows) < 2) next

      # The two parents in this mating
      parent_internal_ids <- dat$ID[mate_rows]

      # Find children whose dad.id AND mom.id are BOTH in this parent pair
      child_rows <- which(
        dat$dad.id %in% parent_internal_ids &
        dat$mom.id %in% parent_internal_ids
      )

      if (length(child_rows) == 0) {
        # Fallback: find children where EITHER parent matches
        child_rows <- which(
          dat$dad.id %in% parent_internal_ids |
          dat$mom.id %in% parent_internal_ids
        )
        # Further filter: children should be one generation below
        if (length(child_rows) > 0) {
          parent_y <- dat$y[mate_rows[1]]
          child_rows <- child_rows[dat$y[child_rows] < parent_y - 0.1]
        }
      }

      if (length(child_rows) > 0) {
        target_x <- mean(dat$x[child_rows], na.rm = TRUE)
        # Set on the row that has mate.centerpoint (first spouse)
        cp_rows <- mate_rows[!is.na(dat$mate.centerpoint[mate_rows])]
        if (length(cp_rows) > 0) {
          dat$descent.target.x[cp_rows[1]] <- target_x
        }
      }
    }
  }

  shape.size <- size
  text.size <- 8
  .pt <- 72.27 / 25.4  # ggplot2 internal constant (points per mm)

  # ---- Tree structure ----
  plt <- ggplot(dat, aes(x = x, y = y)) +
    # Mating lines
    geom_line(
      data = dat[!is.na(dat$mate.id), ],
      aes(group = floor(mate.id)),
      colour = col.tree
    )

  # Double line for consanguineous matings
  if (any(dat$kinship > 0 & !is.na(dat$mate.id))) {
    consang <- dat[dat$mate.id %in% dat$mate.id[dat$kinship > 0] &
                     !is.na(dat$mate.id), ]
    plt <- plt + geom_line(
      data = consang,
      aes(group = mate.id, y = y + 0.02),
      colour = col.tree
    )
    if (plot.kinship.label) {
      consang_label <- consang[!is.na(consang$mate.centerpoint), ]
      if (nrow(consang_label) > 0) {
        plt <- plt + geom_text(
          data = consang_label,
          aes(label = paste0("Kinship:\n", as.character(kinship)),
              x = mate.centerpoint),
          vjust = shape.size / 2 / min(dat$y),
          size = text.size / .pt
        )
      }
    }
  }

  # Repeated subjects
  if (length(dat$Name) != length(unique(dat$Name))) {
    if (!allow.repeated) {
      cli::cli_abort("Drawing subjects multiple times is not supported. Set 'allow.repeated = TRUE'.")
    }
    plt <- plt + geom_line(
      data = dat,
      aes(group = ID, y = y - 0.02),
      linetype = 2, colour = col.double
    )
  }

  # Descending tree segments — use per-mating descent.target.x
  desc_data <- dat[!is.na(dat$mate.id) & !is.na(dat$mate.centerpoint), ]
  if (nrow(desc_data) > 0) {
    plt <- plt +
      geom_segment(
        data = desc_data,
        aes(group = mate.id, x = mate.centerpoint,
            xend = mate.centerpoint, yend = y - 0.5),
        colour = col.tree
      ) +
      geom_segment(
        data = desc_data[!is.na(desc_data$descent.target.x), ],
        aes(group = mate.id, x = mate.centerpoint,
            xend = descent.target.x, y = y - 0.5, yend = y - 0.75),
        colour = col.tree
      )
  }

  # Children connections - twins
  if (any(!is.na(dat$twin.id))) {
    dat$xfambarcoord <- mapply(function(x_val, twin_cp, fam_cp) {
      if (!is.na(twin_cp)) {
        if (x_val >= fam_cp) max(x_val, twin_cp) else min(x_val, twin_cp)
      } else {
        x_val
      }
    }, dat$x, dat$twin.centerpoint, dat$family.centerpoint)

    twin_fam_data <- dat[dat$family %in% dat$family[!is.na(dat$twin.id)], ]
    if (nrow(twin_fam_data) > 0) {
      plt <- plt + geom_line(
        data = twin_fam_data,
        aes(group = family, y = y + 0.25, x = xfambarcoord),
        colour = col.tree
      )
    }

    twin_data <- dat[dat$family %in% dat$family[!is.na(dat$twin.id)] &
                       !is.na(dat$twin.id), ]
    if (nrow(twin_data) > 0) {
      plt <- plt + geom_segment(
        data = twin_data,
        aes(group = family, xend = twin.centerpoint, x = x, yend = y + 0.25),
        colour = col.tree
      )
    }

    twin_type_data <- dat[dat$family %in% dat$family[!is.na(dat$twin.id)] &
                            !is.na(dat$twin.type), ]
    if (nrow(twin_type_data) > 0) {
      plt <- plt + geom_segment(
        data = twin_type_data,
        aes(group = twin.id, linetype = twin.type,
            y = y + 0.125, yend = y + 0.125,
            x = twin.centerpoint + 0.25, xend = twin.centerpoint - 0.25),
        colour = col.tree
      )
    }
  }

  # Non-twin children
  nontwin <- dat[!is.na(dat$family) & is.na(dat$twin.id), ]
  if (nrow(nontwin) > 0) {
    plt <- plt +
      geom_line(data = nontwin, aes(group = family, y = y + 0.25, x = x),
                colour = col.tree) +
      geom_segment(data = nontwin,
                   aes(group = family, xend = x, yend = y + 0.25),
                   colour = col.tree)
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
  if ("Name" %in% colnames(dat) && plot.names) {
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
    theme(legend.position = "none")  # we build legends ourselves

  # ---- Build result object with embedded legend specs ----
  result <- structure(
    list(
      plot     = plt,
      features = features,
      dat      = dat
    ),
    class = "ggped_plot"
  )
  result
}

# =========================================================================
# Legend building helpers
# =========================================================================

#' Build all legend grobs for a ggped plot
#' @keywords internal
.build_all_legend_grobs <- function(features, dat) {
  all_legends <- list()

  # Shape legend
  shape_df <- data.frame(
    x     = c(0, 0, 0),
    y     = c(0, 0, 0),
    shape = factor(c("Male", "Female", "Unknown"),
                    levels = c("Male", "Female", "Unknown"))
  )
  p_shape <- ggplot(shape_df, aes(x = x, y = y, shape = shape)) +
    geom_point(size = 3, fill = "white", colour = "black") +
    scale_shape_manual(
      values = c(Male = 22, Female = 21, Unknown = 23),
      name   = "Sex"
    ) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9, face = "bold"),
          legend.text  = element_text(size = 8))
  g <- .extract_legend_safe(p_shape)
  if (!is.null(g)) all_legends[["Sex"]] <- g

  # Feature legends
  feat_names <- names(features)
  for (fn in feat_names) {
    scale <- features[[fn]]
    display_name <- if (!is.null(scale$name)) scale$name else fn

    if (inherits(scale, "scale_feature_continuous")) {
      vals <- as.numeric(dat[[fn]])
      lims <- scale$limits
      if (is.null(lims)) lims <- range(vals, na.rm = TRUE)
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
            barwidth = grid::unit(3, "cm"),
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
          guide = guide_legend(title.position = "top")
        ) +
        theme_void() +
        theme(legend.position = "bottom",
              legend.title = element_text(size = 9, face = "bold"),
              legend.text  = element_text(size = 8))
    }

    g <- .extract_legend_safe(p_leg)
    if (!is.null(g)) all_legends[[fn]] <- g
  }

  all_legends
}

#' Safely extract legend grob from a ggplot
#' @keywords internal
.extract_legend_safe <- function(p) {
  tryCatch({
    g <- ggplot2::ggplotGrob(p)
    idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    if (length(idx) == 0) return(NULL)
    g$grobs[[idx[1]]]
  }, error = function(e) NULL)
}

#' Arrange legend grobs into rows with automatic wrapping
#' @param legend_grobs List of grob objects
#' @param max_per_row Maximum legends per row (default 3)
#' @return A gtable grob
#' @keywords internal
.arrange_legends_multirow <- function(legend_grobs, max_per_row = 3) {
  n <- length(legend_grobs)
  if (n == 0) return(grid::nullGrob())

  n_rows <- ceiling(n / max_per_row)
  row_height <- grid::unit(1.5, "cm")

  # Build a gtable with n_rows rows
  gt <- gtable::gtable(
    widths  = rep(grid::unit(1, "null"), max_per_row),
    heights = rep(row_height, n_rows)
  )

  for (i in seq_len(n)) {
    r <- ceiling(i / max_per_row)
    c <- ((i - 1) %% max_per_row) + 1
    gt <- gtable::gtable_add_grob(gt, legend_grobs[[i]], t = r, l = c)
  }

  gt
}

# =========================================================================
# Drawing function (used by both print and save)
# =========================================================================

#' Draw a ggped_plot onto the current graphics device
#' @param x A ggped_plot object
#' @keywords internal
.draw_ggped <- function(x) {
  legend_grobs <- .build_all_legend_grobs(x$features, x$dat)
  n_legends <- length(legend_grobs)

  if (n_legends == 0) {
    # No legends — just draw the plot
    grid::grid.draw(ggplot2::ggplotGrob(x$plot))
    return(invisible(NULL))
  }

  legend_table <- .arrange_legends_multirow(legend_grobs, max_per_row = 3)
  n_legend_rows <- ceiling(n_legends / 3)
  legend_height <- grid::unit(1.5 * n_legend_rows, "cm")

  # Create a viewport layout: plot on top, legends on bottom
  main_vp <- grid::viewport(
    layout = grid::grid.layout(
      nrow = 2, ncol = 1,
      heights = grid::unit.c(grid::unit(1, "null"), legend_height)
    )
  )

  grid::pushViewport(main_vp)

  # Draw main plot
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(ggplot2::ggplotGrob(x$plot))
  grid::popViewport()

  # Draw legends
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.draw(legend_table)
  grid::popViewport()

  grid::popViewport()
  invisible(NULL)
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
#' @export
save_ggped <- function(plot, filename, width = 10, height = 8, dpi = 150) {
  if (!inherits(plot, "ggped_plot")) {
    # Fall back to ggsave for standard ggplot objects
    ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi)
    return(invisible(filename))
  }

  ext <- tolower(tools::file_ext(filename))
  if (ext %in% c("png", "")) {
    if (ext == "") filename <- paste0(filename, ".png")
    grDevices::png(filename, width = width, height = height,
                   units = "in", res = dpi)
  } else if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height)
  } else if (ext == "svg") {
    grDevices::svg(filename, width = width, height = height)
  } else if (ext %in% c("tif", "tiff")) {
    grDevices::tiff(filename, width = width, height = height,
                    units = "in", res = dpi)
  } else if (ext %in% c("jpg", "jpeg")) {
    grDevices::jpeg(filename, width = width, height = height,
                    units = "in", res = dpi)
  } else {
    grDevices::png(filename, width = width, height = height,
                   units = "in", res = dpi)
  }

  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  .draw_ggped(plot)
  invisible(filename)
}
