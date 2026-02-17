library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ggplot2)

# ----------------------------
# 1) Subset samples
# ----------------------------
keep <- meta$ExpType == "atRAL5hr"
se_atRAL5hr <- data_filt[, keep]

meta_sub <- as.data.frame(colData(se_atRAL5hr)) %>%
  rownames_to_column("sample")

# ----------------------------
# 2) Extract matrix (PRE-imputation)
# ----------------------------
mat <- assay(se_atRAL5hr)

# Keep proteins with at least 2 observed values
mat <- mat[rowSums(is.finite(mat) & !is.na(mat)) >= 2, , drop = FALSE]

# Remove proteins that still contain NA (prcomp cannot handle NA)
mat <- mat[rowSums(is.na(mat) | !is.finite(mat)) == 0, , drop = FALSE]

# Remove zero-variance proteins
rsd <- apply(mat, 1, sd)
mat <- mat[is.finite(rsd) & rsd > 0, , drop = FALSE]

cat("Proteins used for PCA:", nrow(mat), "\n")

# ----------------------------
# 3) PCA
# ----------------------------
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta_sub, by = "sample") %>%
  mutate(condition = factor(condition))

# ----------------------------
# 4) Plot with ellipses
# ----------------------------
ggplot(scores, aes(PC1, PC2, color = Treatment, fill = Treatment)) +
  stat_ellipse(
    geom  = "polygon",
    type  = "t",
    level = 0.68,
    alpha = 0.2,
    color = NA
  ) +
  stat_ellipse(
    geom  = "path",
    type  = "t",
    linewidth = 1
  ) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA (Pre-imputation): RDH12 atRAL5hr",
    subtitle = paste0("Proteins used: ", nrow(mat)),
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5))
