
#requires table made from 1_initial_processing_MRMdata_script that produces Recovery_lipidsall_v9

experimental_design<- read.csv("Lipidomics/metadata/Sample_description guide_combined.csv")
colnames(experimental_design)



data_abundance <- Lipidsall_long_v7 %>%
  select(picked_candidate, mrm, sample, absint) %>%
  pivot_wider(
    names_from = sample,
    values_from = absint,
    names_prefix = "LFQ."
  )

data_unique <- make_unique(data_abundance, "picked_candidate", "mrm", delim = ";")

experimental_design2 <- experimental_design %>%
  transmute(
    label = paste0("LFQ.", sample),
    condition = paste(Genetype, Treatment,Experiment ,sep = "_"),
    replicate = Rep
  )

LFQ_columns <- grep("LFQ\\.", colnames(data_unique))

data_se <- make_se(data_unique, LFQ_columns, experimental_design2)
data_se

setdiff(experimental_design2$label, colnames(data_unique))
setdiff(colnames(data_unique)[LFQ_columns], experimental_design2$label)


plot_frequency(data_se)

summary(data_unique[, LFQ_columns])
plot_normalization(data_se)

plot_imputation(data_se)
meanSdPlot(data_se)



data_normvsn<- normalize_vsn(data_se)
meanSdPlot(data_normvsn)
plot_normalization(data_se,data_normvsn)
plot_imputation(data_se, data_normvsn)

plot_missval(data_normvsn)
data_imp <- impute(data_normvsn, fun = "QRILC")

plot_imputation(data_se, data_normvsn, data_imp)
plot_normalization(data_se,data_normvsn, data_imp)
meanSdPlot(data_imp)


table(colData(data_normvsn)$condition)
colData(data_normvsn)$condition
#get correlation between mean and std dev
mat <- assay(data_se)

means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)

cor(means, sds)

mat <- assay(data_normvsn)

means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)

cor(means, sds)

mat <- assay(data_imp)

means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)

cor(means, sds)


test<-get_df_long(data_imp)
testnorm<-get_df_long(data_normvsn)


#====================
library(sva)

mat <- assay(data_imp)
meta <- as.data.frame(colData(data_imp))

batch <- meta$Experiment
mod <- model.matrix(~ condition, data = meta)

mat_combat <- ComBat(
  dat = as.matrix(mat),
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

library(dplyr)
library(tidyr)
library(DEP)
library(vsn)
library(SummarizedExperiment)
library(sva)

# requires table made from 1_initial_processing_MRMdata_script
# that produces Lipidsall_long_v7

experimental_design <- read.csv("Lipidomics/metadata/Sample_description guide_combined.csv")

data_abundance <- Lipidsall_long_v7 %>%
  select(picked_candidate, mrm, sample, absint) %>%
  pivot_wider(
    names_from = sample,
    values_from = absint,
    names_prefix = "LFQ."
  )

data_unique <- make_unique(data_abundance, "picked_candidate", "mrm", delim = ";")

# Keep Experiment separate from biological condition
experimental_design2 <- experimental_design %>%
  transmute(
    label = paste0("LFQ.", sample),
    condition = paste(Genetype, Treatment, sep = "_"),
    replicate = Rep,
    Genetype = Genetype,
    Treatment = Treatment,
    Experiment = Experiment
  )

LFQ_columns <- grep("LFQ\\.", colnames(data_unique))

data_se <- make_se(data_unique, LFQ_columns, experimental_design2)

# QC checks
setdiff(experimental_design2$label, colnames(data_unique))
setdiff(colnames(data_unique)[LFQ_columns], experimental_design2$label)

plot_frequency(data_se)
summary(data_unique[, LFQ_columns])
plot_normalization(data_se)
plot_imputation(data_se)
vsn::meanSdPlot(assay(data_se))

# VSN normalization
data_normvsn <- normalize_vsn(data_se)

vsn::meanSdPlot(assay(data_normvsn))
plot_normalization(data_se, data_normvsn)
plot_imputation(data_se, data_normvsn)

# Missingness / imputation
plot_missval(data_normvsn)
data_imp <- impute(data_normvsn, fun = "QRILC")

plot_imputation(data_se, data_normvsn, data_imp)
plot_normalization(data_se, data_normvsn, data_imp)
vsn::meanSdPlot(assay(data_imp))

# Inspect conditions
table(colData(data_normvsn)$condition)
colData(data_normvsn)$condition

# Correlation between row means and SDs
mat <- assay(data_se)
means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)
cor(means, sds)

mat <- assay(data_normvsn)
means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)
cor(means, sds)

mat <- assay(data_imp)
means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)
cor(means, sds)

test <- get_df_long(data_imp)
testnorm <- get_df_long(data_normvsn)

# ====================
# OPTIONAL: ComBat for exploratory PCA only
# Do not use as your main inference workflow if Acute/Recovery is biologically distinct

library(dplyr)
library(SummarizedExperiment)
library(sva)
library(ggplot2)

# start from the imputed object
meta <- as.data.frame(colData(data_imp))
mat  <- assay(data_imp)

# keep only non-Control samples
keep <- meta$Genetype != "Control"

meta_sub <- meta[keep, , drop = FALSE]
mat_sub  <- mat[, keep, drop = FALSE]

# make sure row order in metadata matches matrix columns
stopifnot(ncol(mat_sub) == nrow(meta_sub))

# optional but strongly recommended:
# use sample names from the matrix as rownames of metadata
rownames(meta_sub) <- colnames(mat_sub)

# check balance across batch and treatment
table(meta_sub$Experiment, meta_sub$Treatment)

# factors
meta_sub$Experiment <- factor(meta_sub$Experiment)
meta_sub$Treatment  <- factor(meta_sub$Treatment)

# ComBat
batch <- meta_sub$Experiment
mod <- model.matrix(~ Treatment, data = meta_sub)

mat_combat <- ComBat(
  dat = as.matrix(mat_sub),
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

# PCA
pca <- prcomp(t(mat_combat), center = TRUE, scale. = TRUE)

# scores table
scores <- as.data.frame(pca$x)
scores$sample <- rownames(scores)

meta_sub2 <- meta_sub %>%
  tibble::rownames_to_column("sample")

scores <- scores %>%
  left_join(meta_sub2, by = "sample")

ggplot(scores, aes(PC1, PC2, color = Treatment, shape = Experiment)) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)")
  )


##make PCA with ellipses
make_ellipse <- function(df, level = 0.68, n = 150, ridge = 1e-6) {
  if (nrow(df) < 3) return(NULL)
  
  mu <- colMeans(df[, c("PC1", "PC2")])
  
  S <- stats::cov(df[, c("PC1", "PC2")])
  S <- S + diag(ridge, 2)
  
  r2 <- stats::qchisq(level, df = 2)
  eig <- eigen(S)
  
  A <- eig$vectors %*% diag(sqrt(eig$values * r2)) %*% t(eig$vectors)
  
  tt <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(tt), sin(tt))
  
  ell <- sweep(circle %*% t(A), 2, mu, "+")
  
  tibble(
    PC1 = ell[, 1],
    PC2 = ell[, 2]
  )
}

scores <- scores %>%
  mutate(group = interaction(Treatment, Experiment, drop = TRUE))

ellipse_df <- scores %>%
  group_by(Treatment, Experiment, group) %>%
  group_modify(~ {
    e <- make_ellipse(.x, level = 0.68, n = 150, ridge = 1e-6)
    if (is.null(e)) return(tibble())
    e
  }) %>%
  ungroup()

ggplot(scores, aes(PC1, PC2, color = Treatment, shape = Experiment)) +
  geom_polygon(
    data = ellipse_df,
    aes(group = group, fill = Treatment),
    alpha = 0.18,
    color = NA
  ) +
  geom_path(
    data = ellipse_df,
    aes(group = group, color = Treatment),
    linewidth = 0.6
  ) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA",
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)"),
    color = "Treatment",
    fill  = "Treatment",
    shape = "Experiment"
  ) +
  theme(plot.title = element_text(hjust = 0.5))



# keep same samples used in ComBat
data_combat <- data_imp[, keep]

# replace assay with corrected matrix
assay(data_combat) <- mat_combat

plot_imputation(data_se, data_normvsn, data_imp,data_combat)
plot_normalization(data_se, data_normvsn, data_imp,data_combat)

test<-get_df_wide(data_combat)
