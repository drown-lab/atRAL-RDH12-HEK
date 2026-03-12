library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ggplot2)


meta_sub <- as.data.frame(colData(data_imp)) %>%
  rownames_to_column("sample")

# ----------------------------
# 2) Extract matrix (PRE-imputation)
# ----------------------------
mat <- assay(data_imp)

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
ggplot(scores, aes(PC1, PC2, color = condition, fill = condition)) +
  stat_ellipse(
    geom  = "polygon",
    type  = "t",
    level = 0.68,
    alpha = 0.2,
    color = NA
  ) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA (Pre-imputation): Control and RDH12 Recovery",
    subtitle = paste0("Proteins used: ", nrow(mat)),
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5))


#PCA plot with ellipse


## 1) Clean scores (as you already did)
scores_clean <- scores |>
  as.data.frame() |>
  tibble::as_tibble() |>
  mutate(
    PC1 = as.numeric(PC1),
    PC2 = as.numeric(PC2),
    condition = as.factor(condition)
  ) |>
  filter(!is.na(PC1), !is.na(PC2))

## 2) Helper to generate ellipse points from mean + covariance
make_ellipse <- function(df, level = 0.68, n = 100, ridge = 1e-8) {
  # Need at least 3 points
  if (nrow(df) < 3) return(NULL)
  
  mu <- colMeans(df[, c("PC1", "PC2")])
  
  S <- stats::cov(df[, c("PC1", "PC2")])
  # Stabilize covariance for n=3 (prevents singular matrix issues)
  S <- S + diag(ridge, 2)
  
  # Radius for given confidence level (chi-square with 2 df)
  r2 <- stats::qchisq(level, df = 2)
  
  # Ellipse parameterization using eigen decomposition
  eig <- eigen(S)
  A <- eig$vectors %*% diag(sqrt(eig$values * r2)) %*% t(eig$vectors)
  
  t <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(t), sin(t))
  
  ell <- sweep(circle %*% t(A), 2, mu, "+")
  
  tibble(PC1 = ell[, 1], PC2 = ell[, 2])
}

## 3) Build ellipse dataframe per condition
ellipse_df <- scores_clean |>
  group_by(condition) |>
  group_modify(~ {
    e <- make_ellipse(.x, level = 0.68, n = 150, ridge = 1e-6)
    if (is.null(e)) return(tibble())
    e
  }) |>
  ungroup()

## 4) Plot: shaded ellipses + points
ggplot(scores_clean, aes(PC1, PC2, color = condition, fill = condition)) +
  geom_polygon(
    data = ellipse_df,
    aes(group = condition),
    alpha = 0.18,
    color = NA
  ) +
  geom_path(
    data = ellipse_df,
    aes(group = condition),
    linewidth = 0
  ) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA (After Imputation): Control and RDH12 Recovery",
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)"),
    color = "Condition",
    fill  = "Condition"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_paletteer_d("ggsci::category10_d3") +
  scale_fill_paletteer_d("ggsci::category10_d3")

ggsave(filename = file.path(path_figures,"PCA_Recovery.png"), width = 7, height = 5)

