library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
#load file

df<-Hek_norm

df <- df %>%
  left_join(
    experimental_design,
    by = c("label" = "label")
  )


stopifnot(exists("df"))

# Detection rule: counts 0 as NOT detected (common for proteomics exports)
detected_fun <- function(x) !is.na(x) & is.finite(x) & x > 0
print(unique(df$condition))

df$condition <- gsub("\\+", ".", df$condition)
print(unique(df$condition))

# ----------------------------
# 1) Long table: one row per protein × sample × channel
# ----------------------------
det_long <- df %>%
  mutate(
    detected = detected_fun(PG.MaxLFQ),
    log2_intensity = log2(PG.MaxLFQ),  # will be -Inf for 0; we'll handle below
    Protein.Group = Genes
  )


# Replace non-finite log2 values (0, NA, Inf) with NA for abundance summaries
det_long <- det_long %>%
  mutate(log2_intensity = ifelse(is.finite(log2_intensity), log2_intensity, NA_real_))


# ----------------------------
# 3) "Global abundance" per protein × condition 
#    Use mean log2 intensity across replicates (only where detected/log2 finite)
# ----------------------------
abundance_df <- det_long %>%
  group_by(, Protein.Group, condition) %>%
  summarise(
    mean_log2_intensity = mean(log2_intensity, na.rm = TRUE),
    .groups = "drop"
  )


# ----------------------------
# 3a) Detection class per protein × condition  (replicate-based)
# ----------------------------
#came from Imputation model

plot_df <- protein_class_by_cond %>%
  left_join(abundance_df, by = c( "Protein.Group", "condition")) %>%
  filter(!is.na(mean_log2_intensity))


# ----------------------------
# 4) Density plot: distribution of mean log2 intensity by class
# ----------------------------
p_den <- ggplot(plot_df, aes(x = mean_log2_intensity, fill = missing_class_3rep)) +
  geom_density(alpha = 0.4, na.rm = TRUE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Mean log2(Intensity) density by detection class",
    x = "Mean log2(Intensity) across replicates",
    y = "Density",
    fill = "Detection class"
  )

print(p_den)


p_den <- ggplot(plot_df, aes(x = mean_log2_intensity, fill = missing_class_3rep )) +
  geom_density(alpha = 0.4, na.rm = TRUE) +
  facet_wrap(~ condition) +
  theme_bw(base_size = 12) +
  labs(
    title = "Mean log2(Intensity) density by detection class",
    x = "Mean log2(Intensity) across replicates",
    y = "Density",
    fill = "Detection class"
  )

print(p_den)

ggsave(filename = file.path(path_figures,"Missingness_Abund_Distribution.png"), width = 8, height = 7)

