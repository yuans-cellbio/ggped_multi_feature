#' Align a pedigree and retain its explicit relationships
#'
#' Uses kinship2 for initial coordinates. A person can have several drawing
#' occurrences. Parentage is retained separately from per-row family columns.
#' @param ped A single kinship2 pedigree.
#' @param chrtype Chromosome type passed to kinship2::kinship.
#' @param packed,width,align,hints Arguments passed to kinship2::align.pedigree.
#' @return A data frame with legacy drawing columns, stable `person_id` and
#'   `occurrence_id` columns, and a `ggped_relationships` attribute.
#' @export
#' @examples
#' ped <- kinship2::pedigree(1:3, c(0, 0, 1), c(0, 0, 2), c(1, 2, 1))
#' dfalign.pedigree(ped)
dfalign.pedigree <- function(ped, chrtype = "autosome", packed = TRUE,
                             width = 10, align = TRUE, hints = ped$hints) {
  if (!inherits(ped, "pedigree")) cli::cli_abort("'ped' must be a single kinship2 pedigree.")
  chrtype <- tolower(chrtype)
  if (length(chrtype) != 1L || is.na(chrtype) || !chrtype %in% c("autosome", "x", "y", "mt"))
    cli::cli_abort("'chrtype' must be autosome, x, y, or mt.")
  .check_flag(packed, "packed"); .check_flag(align, "align"); .check_positive(width, "width")
  if (!is.null(hints) && !is.list(hints)) cli::cli_abort("'hints' must be a list or NULL.")
  if (!length(ped$id)) cli::cli_abort("'ped' must contain individuals.")
  s <- kinship2::align.pedigree(ped, packed = packed, width = width, align = align, hints = hints)
  cells <- which(s$nid > 0, arr.ind = TRUE)
  cells <- cells[order(cells[, 1], cells[, 2]), , drop = FALSE]
  ids <- s$nid[cells]; n <- length(ids)
  d <- data.frame(ID = ids, Name = as.character(ped$id[ids]),
    family = ifelse(s$fam[cells] > 0, paste(cells[, 1], s$fam[cells], sep = ":"), NA),
    mate.id = NA_real_, mate.type = factor(s$spouse[cells]), twin.id = NA_real_, twin.type = NA_real_,
    dad.id = ped$findex[ids], mom.id = ped$mindex[ids], sex = ped$sex[ids],
    status = if (is.null(ped$status)) rep(0, n) else ped$status[ids], x = s$pos[cells], y = cells[, 1],
    mate.centerpoint = NA_real_, family.centerpoint = NA_real_, twin.centerpoint = NA_real_, kinship = 0,
    person_id = as.character(ped$id[ids]), stringsAsFactors = FALSE)
  if (is.null(ped$affected)) d$affected <- FALSE
  else if (is.matrix(ped$affected)) {
    a <- ped$affected[ids, , drop = FALSE]
    if (is.null(colnames(a))) colnames(a) <- paste0("affected", seq_len(ncol(a)))
    for (nm in colnames(a)) d[[nm]] <- as.logical(a[, nm])
  } else d$affected <- as.logical(ped$affected[ids])
  people <- data.frame(person_id = as.character(ped$id), ID = seq_along(ped$id),
    dad = ifelse(ped$findex > 0, as.character(ped$id[pmax(1L, ped$findex)]), NA_character_),
    mom = ifelse(ped$mindex > 0, as.character(ped$id[pmax(1L, ped$mindex)]), NA_character_))
  pairs <- unique(people[!is.na(people$dad), c("dad", "mom")]); names(pairs) <- c("parent1", "parent2")
  twins <- data.frame(person1 = character(), person2 = character(), code = integer())
  if (!is.null(ped$relation)) {
    r <- ped$relation; code <- as.integer(r$code)
    if (any(code == 4L)) pairs <- unique(rbind(pairs, data.frame(parent1 = as.character(ped$id[r$indx1[code == 4L]]), parent2 = as.character(ped$id[r$indx2[code == 4L]]))))
    if (any(code < 4L)) twins <- data.frame(person1 = as.character(ped$id[r$indx1[code < 4L]]), person2 = as.character(ped$id[r$indx2[code < 4L]]), code = code[code < 4L])
  }
  if (!is.null(hints$spouse)) pairs <- unique(rbind(pairs, data.frame(parent1 = as.character(ped$id[hints$spouse[, 1]]), parent2 = as.character(ped$id[hints$spouse[, 2]]))))
  pairs$relationship_id <- .pair_id(pairs$parent1, pairs$parent2)
  pairs <- pairs[!duplicated(pairs$relationship_id), , drop = FALSE]
  ks <- kinship2::kinship(ped, chrtype)
  pairs$kinship <- vapply(seq_len(nrow(pairs)), function(i) ks[match(pairs$parent1[i], ped$id), match(pairs$parent2[i], ped$id)], numeric(1))
  d <- .adjacent_partners(d, pairs)
  d$occurrence_id <- .occurrence_ids(d$person_id)
  attr(d, "ggped_relationships") <- list(people = people, partnerships = pairs, twins = twins)
  model <- .relationship_model(d)
  for (i in seq_len(nrow(model$partnerships))) {
    p <- model$partnerships[i, ]; rows <- match(c(p$occurrence1, p$occurrence2), d$occurrence_id)
    d$mate.id[rows] <- i; d$mate.centerpoint[rows[1]] <- p$x; d$kinship[rows[1]] <- p$kinship
  }
  for (i in seq_len(nrow(model$families))) {
    f <- model$families[i, ]; cr <- match(model$children$occurrence_id[model$children$family_id == f$family_id], d$occurrence_id)
    d$family[cr] <- f$family_id
    pr <- match(c(f$occurrence1, f$occurrence2), d$occurrence_id); d$family.centerpoint[pr] <- f$target_x
  }
  for (i in seq_len(nrow(model$twin_members))) {
    t <- model$twin_members[i, ]; row <- match(t$occurrence_id, d$occurrence_id)
    d$twin.id[row] <- match(t$twin_id, unique(model$twin_members$twin_id)); d$twin.type[row] <- t$code; d$twin.centerpoint[row] <- t$apex_x
  }
  d$family <- factor(d$family); d$status <- factor(d$status); d$twin.type <- factor(d$twin.type)
  omitted <- setdiff(as.character(ped$id), d$person_id)
  if (length(omitted)) cli::cli_alert_info("No connected layout for subject(s): {paste(omitted, collapse = ', ')}.")
  d
}
