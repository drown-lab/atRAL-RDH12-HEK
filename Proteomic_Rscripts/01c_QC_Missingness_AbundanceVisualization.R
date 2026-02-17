library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(tibble)

# =========================
# Choose SE (pre-imputation)
# =========================
se <- data_filt

mat  <- assay(se)
meta <- as.data.frame(colData(se)) %>%
  rownames_to_column("sample")

stopifnot("condition" %in% names(meta))

# =========================================================
# 1) Detection class per protein × condition (0/1/2/3 of reps)
#    IMPORTANT: if zeros represent missing in your data, convert 0 -> NA here
# =========================================================
mat2 <- mat
mat2[is.finite(mat2) & mat2 == 0] <- NA_real_  # remove if zeros are real values

det_by_cond <- as.data.frame(!is.na(mat2)) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample", values_to = "detected") %>%
  left_join(meta %>% select(sample, condition), by = "sample") %>%
  group_by(protein, condition) %>%
  summarise(
    n_detect = sum(detected),
    n_reps   = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    detect_class = case_when(
      n_detect == 0 ~ "MNAR_0of3",
      n_detect == 1 ~ "MAR_1of3",
      n_detect == 2 ~ "MAR_2of3",
      n_detect >= 3 ~ "Present_3of3",
      TRUE ~ NA_character_
    ),
    detect_class = factor(
      detect_class,
      levels = c("MNAR_0of3", "MAR_1of3", "MAR_2of3", "Present_3of3")
    )
  )

# =========================================================
# 2) Mean log2 abundance per protein × condition (across reps)
# =========================================================
mean_abund <- as.data.frame(mat2) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample", values_to = "log2_abundance") %>%
  left_join(meta %>% select(sample, condition), by = "sample") %>%
  group_by(protein, condition) %>%
  summarise(
    mean_log2_abundance = mean(log2_abundance, na.rm = TRUE),
    .groups = "drop"
  )

# Join mean abundance onto detection class table
plot_df <- det_by_cond %>%
  left_join(mean_abund, by = c("protein", "condition")) %>%
  filter(!is.na(detect_class), is.finite(mean_log2_abundance))

# =========================================================
# 3) Density plot: mean log2 abundance by MAR class
# =========================================================
p_den <- ggplot(plot_df, aes(x = mean_log2_abundance, fill = detect_class)) +
  geom_density(alpha = 0.4, na.rm = TRUE) +
  facet_wrap(~ condition, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Mean log2 abundance density by detection class (pre-imputation)",
    x = "Mean log2 abundance across replicates",
    y = "Density",
    fill = "Detection class"
  )

print(p_den)
