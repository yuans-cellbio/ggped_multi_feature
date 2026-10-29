#' @importFrom stats ave setNames
#' @importFrom utils head tail
NULL

# ggplot2 evaluates these names in layer data rather than the package namespace.
utils::globalVariables(c("x", "y", "xend", "yend", "shape", "val", "text",
  "relationship_id", "linetype", "colour", "kinship", "sex", "status",
  "individual_id", "adopted", "pregnancy", "Name", "occurrence_id"))
