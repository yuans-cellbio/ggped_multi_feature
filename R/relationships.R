.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) cli::cli_abort("'{name}' must be TRUE or FALSE.")
}
.check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) cli::cli_abort("'{name}' must be a positive finite number.")
}
.pair_id <- function(a, b) {
  if (!length(a) || !length(b)) return(character())
  a <- as.character(a); b <- as.character(b)
  lo <- ifelse(a < b, a, b); hi <- ifelse(a < b, b, a)
  paste0("pair:", nchar(lo), ":", lo, ":", nchar(hi), ":", hi)
}
.occurrence_ids <- function(ids) paste0("occ:", nchar(ids), ":", ids, ":", ave(seq_along(ids), ids, FUN = seq_along))

.adjacent_partners <- function(d, pairs) {
  for (i in seq_len(nrow(pairs))) {
    a <- which(d$person_id == pairs$parent1[i]); b <- which(d$person_id == pairs$parent2[i])
    if (!length(a) || !length(b)) cli::cli_abort("The initial layout omitted a partnership endpoint.")
    candidates <- expand.grid(a = a, b = b)
    same <- candidates[d$y[candidates$a] == d$y[candidates$b], , drop = FALSE]
    if (!nrow(same)) {
      row <- if (min(d$y[a]) < min(d$y[b])) a[1] else b[1]
      other <- if (row %in% a) b[which.max(d$y[b])] else a[which.max(d$y[a])]
    } else {
      distances <- abs(d$x[same$a] - d$x[same$b]); j <- which.min(distances)
      between <- which(d$y == d$y[same$a[j]] & d$x > min(d$x[c(same$a[j], same$b[j])]) + 1e-7 & d$x < max(d$x[c(same$a[j], same$b[j])]) - 1e-7)
      if (!length(between)) next
      degree <- function(id) sum(pairs$parent1 == id | pairs$parent2 == id)
      row <- if (degree(pairs$parent1[i]) >= degree(pairs$parent2[i])) same$a[j] else same$b[j]
      other <- if (row == same$a[j]) same$b[j] else same$a[j]
    }
    pivot <- d$x[other]
    move <- d$x >= pivot - 0.5 - 1e-7; d$x[move] <- d$x[move] + 1
    copy <- d[row, , drop = FALSE]; copy$x <- d$x[other] - 1; copy$y <- d$y[other]; copy$family <- NA_character_
    d <- rbind(d, copy)
  }
  d[order(d$y, d$x), , drop = FALSE]
}

.legacy_relationships <- function(d) {
  if (!"ID" %in% names(d)) {
    if (anyDuplicated(d$Name)) cli::cli_abort("Repeated legacy rows require stable 'ID' values.")
    d$ID <- seq_len(nrow(d))
  }
  if (!"person_id" %in% names(d)) d$person_id <- as.character(d$ID)
  if (!all(c("dad.id", "mom.id") %in% names(d))) {
    if ("family" %in% names(d) && any(!is.na(d$family) & as.character(d$family) != "0")) cli::cli_abort("Ambiguous legacy data: child families require both 'dad.id' and 'mom.id'. Recreate with dfalign.pedigree().")
    d$dad.id <- d$mom.id <- 0
  }
  ids <- as.character(d$ID); lookup <- setNames(d$person_id[!duplicated(ids)], ids[!duplicated(ids)])
  parent <- function(x) { x <- as.character(x); x[is.na(x) | x == "0" | x == ""] <- NA_character_; unname(lookup[x]) }
  missing_parent <- function(x) !is.na(x) & !as.character(x) %in% c("0", "") & !as.character(x) %in% names(lookup)
  if (any(missing_parent(d$dad.id) | missing_parent(d$mom.id))) cli::cli_abort("Parent IDs must refer to existing 'ID' values.")
  people <- unique(data.frame(person_id = d$person_id, ID = d$ID, dad = parent(d$dad.id), mom = parent(d$mom.id)))
  if (anyDuplicated(people$person_id)) cli::cli_abort("Repeated people have conflicting parentage.")
  if (any(xor(is.na(people$dad), is.na(people$mom)))) cli::cli_abort("Each person must have two parents or neither parent.")
  pairs <- unique(people[!is.na(people$dad), c("dad", "mom")]); names(pairs) <- c("parent1", "parent2")
  if ("mate.id" %in% names(d)) for (mid in unique(stats::na.omit(d$mate.id))) {
    rows <- which(d$mate.id == mid); members <- unique(d$person_id[rows])
    if (length(members) != 2L) cli::cli_abort("Ambiguous legacy mating '{mid}': expected exactly two people. Recreate with dfalign.pedigree().")
    key <- .pair_id(members[1], members[2]); known <- .pair_id(pairs$parent1, pairs$parent2)
    if (!key %in% known) {
      if (any(members %in% c(pairs$parent1, pairs$parent2))) cli::cli_abort("Ambiguous legacy mating '{mid}' conflicts with supplied parent pairs. Recreate with dfalign.pedigree().")
      pairs <- rbind(pairs, data.frame(parent1 = members[1], parent2 = members[2]))
    }
  }
  pairs$relationship_id <- .pair_id(pairs$parent1, pairs$parent2)
  pairs$kinship <- vapply(seq_len(nrow(pairs)), function(i) {
    rows <- d$person_id %in% c(pairs$parent1[i], pairs$parent2[i])
    if (!"kinship" %in% names(d)) return(0)
    vals <- d$kinship[rows]; if (all(is.na(vals))) 0 else max(vals, na.rm = TRUE)
  }, numeric(1))
  twins <- data.frame(person1 = character(), person2 = character(), code = integer())
  if ("twin.id" %in% names(d)) for (tid in unique(stats::na.omit(d$twin.id))) {
    rows <- which(d$twin.id == tid); rows <- rows[order(d$x[rows])]
    if (length(rows) < 2L) cli::cli_abort("Legacy twin groups require at least two members.")
    codes <- unique(stats::na.omit(as.integer(as.character(d$twin.type[rows]))))
    if (length(codes) != 1L) cli::cli_abort("Ambiguous legacy twin type. Recreate with dfalign.pedigree().")
    twins <- rbind(twins, data.frame(person1 = head(d$person_id[rows], -1), person2 = tail(d$person_id[rows], -1), code = codes))
  }
  list(people = people, partnerships = pairs, twins = twins)
}

.relationship_model <- function(d) {
  required <- c("x", "y", "Name", "sex")
  if (!all(required %in% names(d))) cli::cli_abort("Drawing data require columns: {paste(required, collapse = ', ')}.")
  if (!nrow(d) || !is.numeric(d$x) || !is.numeric(d$y) || any(!is.finite(d$x) | !is.finite(d$y))) cli::cli_abort("Drawing coordinates must be finite numeric values.")
  source <- attr(d, "ggped_relationships")
  if (is.null(source)) source <- .legacy_relationships(d)
  if (!"person_id" %in% names(d)) d$person_id <- as.character(d$ID)
  if (!"occurrence_id" %in% names(d)) d$occurrence_id <- .occurrence_ids(d$person_id)
  if (anyNA(d$person_id) || anyDuplicated(d$occurrence_id)) cli::cli_abort("Person IDs must be present and occurrence IDs unique.")
  if (!all(d$person_id %in% source$people$person_id)) cli::cli_abort("Unknown person IDs in drawing data; recreate with dfalign.pedigree().")
  if (all(c("ID", "dad.id", "mom.id") %in% names(d))) {
    ref <- source$people[match(d$person_id, source$people$person_id), ]
    lookup <- setNames(as.character(source$people$ID), source$people$person_id)
    expected_dad <- unname(lookup[ref$dad]); expected_mom <- unname(lookup[ref$mom])
    normalize <- function(v) { v <- as.character(v); v[is.na(v) | v == "0" | v == ""] <- NA_character_; v }
    equal <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    if (any(as.character(d$ID) != as.character(ref$ID)) ||
        !all(equal(normalize(d$dad.id), expected_dad)) || !all(equal(normalize(d$mom.id), expected_mom)))
      cli::cli_abort("Parentage or IDs conflict with retained relationships. Recreate with dfalign.pedigree() after editing the pedigree.")
  }
  parts <- source$partnerships
  parts$occurrence1 <- parts$occurrence2 <- character(nrow(parts)); parts$x <- parts$y <- numeric(nrow(parts))
  by_person <- split(seq_len(nrow(d)), d$person_id)
  for (i in seq_len(nrow(parts))) {
    a <- by_person[[parts$parent1[i]]]; b <- by_person[[parts$parent2[i]]]
    if (!length(a) || !length(b)) cli::cli_abort("A partnership endpoint is absent from the drawing.")
    cand <- expand.grid(a = a, b = b)
    score <- 1000 * abs(d$y[cand$a] - d$y[cand$b]) + abs(d$x[cand$a] - d$x[cand$b])
    j <- which.min(score); rows <- c(cand$a[j], cand$b[j])
    parts$occurrence1[i] <- d$occurrence_id[rows[1]]; parts$occurrence2[i] <- d$occurrence_id[rows[2]]
    parts$x[i] <- mean(d$x[rows]); parts$y[i] <- max(d$y[rows])
  }
  people <- source$people
  child_index <- split(which(!is.na(people$dad)), .pair_id(people$dad, people$mom)[!is.na(people$dad)])
  children <- data.frame(person_id = character(), occurrence_id = character(), family_id = character())
  families <- data.frame(family_id = character(), relationship_id = character(), occurrence1 = character(), occurrence2 = character(), x = numeric(), y = numeric(), target_x = numeric(), child_y = numeric())
  child_records <- family_records <- list()
  for (i in seq_len(nrow(parts))) {
    p <- parts[i, ]; child_ids <- people$person_id[child_index[[p$relationship_id]]]
    rows <- which(d$person_id %in% child_ids)
    rows <- rows[order(d$y[rows], d$x[rows])]; rows <- rows[!duplicated(d$person_id[rows])]
    for (gen in unique(d$y[rows])) {
      cr <- rows[d$y[rows] == gen]
      if (gen <= p$y) cli::cli_abort("Children must be below their parent partnership.")
      fid <- paste0(p$relationship_id, "@", gen)
      family_records[[length(family_records) + 1L]] <- data.frame(family_id = fid, relationship_id = p$relationship_id, occurrence1 = p$occurrence1, occurrence2 = p$occurrence2, x = p$x, y = p$y, target_x = mean(range(d$x[cr])), child_y = gen)
      child_records[[length(child_records) + 1L]] <- data.frame(person_id = d$person_id[cr], occurrence_id = d$occurrence_id[cr], family_id = fid)
    }
  }
  if (length(family_records)) families <- do.call(rbind, family_records)
  if (length(child_records)) children <- do.call(rbind, child_records)
  twin_members <- data.frame(twin_id = character(), occurrence_id = character(), family_id = character(), code = integer(), apex_x = numeric())
  twins <- source$twins; twins$twin_id <- character(nrow(twins))
  for (i in seq_len(nrow(twins))) {
    ids <- c(twins$person1[i], twins$person2[i])
    repeat {
      linked <- twins$person1 %in% ids | twins$person2 %in% ids
      new <- unique(c(ids, twins$person1[linked], twins$person2[linked]))
      if (length(new) == length(ids)) break
      ids <- new
    }
    sorted <- sort(ids)
    tid <- paste0("twins:", paste0(nchar(sorted), ":", sorted, collapse = ":")); twins$twin_id[i] <- tid
    if (tid %in% twin_members$twin_id) next
    cr <- children[children$person_id %in% ids, , drop = FALSE]
    if (nrow(cr) != length(ids) || length(unique(cr$family_id)) != 1L) cli::cli_abort("Twin groups must have the same two parents and drawing generation.")
    if (any(!twins$code[linked] %in% 1:3)) cli::cli_abort("Twin codes must be 1, 2, or 3.")
    x <- d$x[match(cr$occurrence_id, d$occurrence_id)]
    twin_members <- rbind(twin_members, data.frame(twin_id = tid, occurrence_id = cr$occurrence_id, family_id = cr$family_id, code = twins$code[i], apex_x = mean(range(x))))
  }
  list(people = people, occurrences = d, partnerships = parts, families = families, children = children, twins = twins, twin_members = twin_members)
}

.connection_segments <- function(model) {
  d <- model$occurrences; segments <- list()
  add <- function(id, kind, x, y, xend, yend, from = NA_character_, to = NA_character_, linetype = "solid") {
    if (abs(x - xend) + abs(y - yend) < 1e-9) return(invisible(NULL))
    segments[[length(segments) + 1L]] <<- list(id, kind, from, to, x, y, xend, yend, linetype)
  }
  for (i in seq_len(nrow(model$partnerships))) {
    p <- model$partnerships[i, ]; rows <- match(c(p$occurrence1, p$occurrence2), d$occurrence_id); a <- rows[1]; b <- rows[2]
    add(p$relationship_id, "partnership", d$x[a], d$y[a], d$x[b], d$y[b], p$occurrence1, p$occurrence2)
    if (p$kinship > 0) add(p$relationship_id, "consanguinity", d$x[a], d$y[a] - 0.02, d$x[b], d$y[b] - 0.02, p$occurrence1, p$occurrence2)
  }
  for (i in seq_len(nrow(model$families))) {
    f <- model$families[i, ]; bar_y <- f$child_y - 0.25
    cr <- model$children[model$children$family_id == f$family_id, , drop = FALSE]; rows <- match(cr$occurrence_id, d$occurrence_id)
    tm <- model$twin_members[model$twin_members$family_id == f$family_id, , drop = FALSE]; anchors <- d$x[rows]
    if (nrow(tm)) anchors[match(tm$occurrence_id, cr$occurrence_id)] <- tm$apex_x
    target <- mean(range(anchors)); mid_y <- f$y + (bar_y - f$y) * 2 / 3
    add(f$family_id, "descent", f$x, f$y, f$x, mid_y)
    add(f$family_id, "descent", f$x, mid_y, target, bar_y)
    add(f$family_id, "sibship", min(anchors), bar_y, max(anchors), bar_y)
    for (j in seq_along(rows)) add(f$family_id, if (cr$occurrence_id[j] %in% tm$occurrence_id) "twin_branch" else "child", anchors[j], bar_y, d$x[rows[j]], d$y[rows[j]], to = cr$occurrence_id[j])
    tr <- model$twins[model$twins$twin_id %in% tm$twin_id, , drop = FALSE]
    for (j in seq_len(nrow(tr))) {
      if (tr$code[j] == 2L) next
      a <- rows[match(tr$person1[j], cr$person_id)]; b <- rows[match(tr$person2[j], cr$person_id)]; apex <- tm$apex_x[match(tr$twin_id[j], tm$twin_id)]
      add(tr$twin_id[j], "zygosity", (apex + d$x[a]) / 2, (bar_y + d$y[a]) / 2, (apex + d$x[b]) / 2, (bar_y + d$y[b]) / 2, linetype = if (tr$code[j] == 1L) "solid" else "dashed")
    }
  }
  for (id in unique(d$person_id[duplicated(d$person_id)])) {
    rows <- which(d$person_id == id); rows <- rows[order(d$y[rows], d$x[rows])]
    for (j in seq_len(length(rows) - 1L)) {
      a <- rows[j]; b <- rows[j + 1L]; offset <- min(d$y[c(a, b)]) - 0.12; rid <- paste0("repeat:", id, ":", j)
      add(rid, "repeated", d$x[a], d$y[a], d$x[a], offset, d$occurrence_id[a], linetype = "dashed")
      add(rid, "repeated", d$x[a], offset, d$x[b], offset, d$occurrence_id[a], d$occurrence_id[b], "dashed")
      add(rid, "repeated", d$x[b], offset, d$x[b], d$y[b], to = d$occurrence_id[b], linetype = "dashed")
    }
  }
  if (!length(segments)) return(data.frame(relationship_id = character(), kind = character(), from_occurrence = character(), to_occurrence = character(), x = numeric(), y = numeric(), xend = numeric(), yend = numeric(), linetype = character()))
  # Construct columns once, avoiding a data.frame allocation for every segment.
  col <- function(i, prototype) vapply(segments, `[[`, prototype, i)
  s <- data.frame(relationship_id = col(1, ""), kind = col(2, ""),
    from_occurrence = col(3, ""), to_occurrence = col(4, ""),
    x = col(5, 0), y = col(6, 0), xend = col(7, 0), yend = col(8, 0), linetype = col(9, ""))
  s[!duplicated(s), , drop = FALSE]
}
