# ===========================================================================
# ggped v2 Test Example
# Demonstrates multi-feature segmented pedigree symbols
# ===========================================================================

# ---- Load dependencies FIRST ----
library(kinship2)
library(ggplot2)
library(grid)
library(gtable)

# ---- Load ggped ----
# Option A: source from local directory
# Set pkg_dir to the ggped package root (parent of R/)
pkg_dir <- getwd()
# If running from inst/ subfolder:
if (basename(pkg_dir) == "inst") pkg_dir <- dirname(pkg_dir)

r_dir <- file.path(pkg_dir, "R")
if (dir.exists(r_dir)) {
  cat("Sourcing R files from:", r_dir, "\n")
  for (f in list.files(r_dir, full.names = TRUE, pattern = "\\.R$")) {
    source(f)
  }
} else {
  # Option B: installed package
  library(ggped)
}


# ===========================================================================
# Example 1: minnbreast data — single discrete feature
# ===========================================================================

cat("\n=== Example 1: Single feature (affected status) ===\n")

data(minnbreast)
bpeds <- with(minnbreast,
  pedigree(id, fatherid, motherid, sex,
           affected = proband, famid = famid))
bped.id8 <- bpeds['8']
df <- dfalign.pedigree(bped.id8)

# Convert affected to character for discrete scale
df$affected <- as.character(df$affected)

p1 <- ggdraw.pedigree(
  dat = df,
  features = list(
    affected = scale_feature_discrete(
      values = c("TRUE" = "black", "FALSE" = "white"),
      na.value = "grey80",
      name = "Proband"
    )
  ),
  size = 7
)

save_ggped(p1, "test_example1_single_feature.png", width = 12, height = 8)
cat("Saved: test_example1_single_feature.png\n")


# ===========================================================================
# Example 2: Synthetic pedigree — 3 features (1 continuous + 2 discrete)
# ===========================================================================

cat("\n=== Example 2: Three features (continuous + discrete) ===\n")

# Build a simple 2-generation pedigree
ped_df <- data.frame(
  id       = c("F1", "M1", "S1", "D1", "S2"),
  fatherid = c(NA, NA, "F1", "F1", "F1"),
  motherid = c(NA, NA, "M1", "M1", "M1"),
  sex      = c(1, 2, 1, 2, 1),
  stringsAsFactors = FALSE
)

ped2 <- with(ped_df,
  pedigree(id = id, dadid = fatherid, momid = motherid, sex = sex))

df2 <- dfalign.pedigree(ped2)

# Add synthetic feature columns
set.seed(42)
df2$disease_score <- c(0.2, 0.1, 0.85, 0.45, 0.6)[seq_len(nrow(df2))]
df2$genotype      <- c("wt/ko", "wt/wt", "ko/ko", "wt/ko", "wt/wt")[seq_len(nrow(df2))]
df2$treatment     <- c("control", "treated", "treated", "control", "treated")[seq_len(nrow(df2))]

p2 <- ggdraw.pedigree(
  dat = df2,
  features = list(
    disease_score = scale_feature_continuous(
      low = "white", high = "red",
      limits = c(0, 1),
      name = "Disease Score"
    ),
    genotype = scale_feature_discrete(
      values = c("wt/wt" = "#f7f7f7", "wt/ko" = "#fc8d59", "ko/ko" = "#d73027"),
      name = "Genotype"
    ),
    treatment = scale_feature_discrete(
      values = c("control" = "#d9d9d9", "treated" = "#3182bd"),
      name = "Treatment"
    )
  ),
  size = 9
)

save_ggped(p2, "test_example2_three_features.png", width = 10, height = 8)
cat("Saved: test_example2_three_features.png\n")


# ===========================================================================
# Example 3: 4 features
# ===========================================================================

cat("\n=== Example 3: Four features ===\n")

df3 <- df2

# Add a 4th feature
df3$pathway_activity <- c(0.3, 0.7, 0.9, 0.5, 0.2)[seq_len(nrow(df3))]

p3 <- ggdraw.pedigree(
  dat = df3,
  features = list(
    disease_score = scale_feature_continuous(
      low = "white", high = "#e41a1c",
      limits = c(0, 1),
      name = "Disease Score"
    ),
    genotype = scale_feature_discrete(
      values = c("wt/wt" = "#f7f7f7", "wt/ko" = "#fc8d59", "ko/ko" = "#d73027"),
      name = "Genotype"
    ),
    treatment = scale_feature_discrete(
      values = c("control" = "#d9d9d9", "treated" = "#3182bd"),
      name = "Treatment"
    ),
    pathway_activity = scale_feature_continuous(
      low = "#f7fcf5", high = "#238b45",
      limits = c(0, 1),
      name = "Pathway Activity"
    )
  ),
  size = 9
)

save_ggped(p3, "test_example3_four_features.png", width = 10, height = 8)
cat("Saved: test_example3_four_features.png\n")


# ===========================================================================
# Example 4: Diverging continuous scale (2 features)
# ===========================================================================

cat("\n=== Example 4: Diverging colour scale ===\n")

df4 <- df2
df4$fold_change <- c(-1.5, 0.2, 2.3, -0.5, 1.1)[seq_len(nrow(df4))]

p4 <- ggdraw.pedigree(
  dat = df4,
  features = list(
    fold_change = scale_feature_continuous(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-3, 3),
      name = "Log2 Fold Change"
    ),
    genotype = scale_feature_discrete(
      values = c("wt/wt" = "#f7f7f7", "wt/ko" = "#fc8d59", "ko/ko" = "#d73027"),
      name = "Genotype"
    )
  ),
  size = 9
)

save_ggped(p4, "test_example4_diverging.png", width = 10, height = 8)
cat("Saved: test_example4_diverging.png\n")


# ===========================================================================
# Segmentation geometry validation
# ===========================================================================

cat("\n=== Segmentation geometry check ===\n")
for (shape in c("circle", "square", "diamond")) {
  for (n in 1:4) {
    segs <- get_segment_polygons(shape, n)
    cat(sprintf("  %s n=%d: %d segments, pts per segment: %s\n",
                shape, n, length(segs),
                paste(sapply(segs, nrow), collapse = ", ")))
  }
}

cat("\nAll tests complete.\n")
