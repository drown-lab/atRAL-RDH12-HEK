library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
#load file

df<-read.csv("Cell_Perturbed/input_files/SILAC_UTandMG132_data.csv")

#df <- df2_MG132

# Detection rule: counts 0 as NOT detected (common for proteomics exports)
detected_fun <- function(x) !is.na(x) & is.finite(x) & x > 0

# ----------------------------
# 1) Long table: one row per protein × sample × channel
# ----------------------------
det_long <- df %>%
  select(Protein.Group, ExperimentName, condition, replicate, L, H) %>%
  pivot_longer(cols = c(L, H), names_to = "channel", values_to = "intensity") %>%
  mutate(
    detected = detected_fun(intensity),
    log2_intensity = log2(intensity)  # will be -Inf for 0; we'll handle below
  )

# Replace non-finite log2 values (0, NA, Inf) with NA for abundance summaries
det_long <- det_long %>%
  mutate(log2_intensity = ifelse(is.finite(log2_intensity), log2_intensity, NA_real_))

# ----------------------------
# 2) Detection class per protein × condition × channel (replicate-based)
# ----------------------------
det_by_cond <- det_long %>%
  group_by(channel, Protein.Group, condition) %>%
  summarise(
    n_detect = sum(detected, na.rm = TRUE),
    n_reps   = n(),  # number of rows (replicates) for that protein in that condition/channel
    .groups  = "drop"
  ) %>%
  mutate(
    detect_class = case_when(
      n_reps == 3 & n_detect == 0 ~ "MNAR_0",
      n_reps == 3 & n_detect == 1 ~ "MAR_1of3",
      n_reps == 3 & n_detect == 2 ~ "MAR_2of3",
      n_reps == 3 & n_detect == 3 ~ "Present_3",
      TRUE ~ NA_character_  # drop conditions that don't have exactly 3 reps
    ),
    detect_class = factor(
      detect_class,
      levels = c("MNAR_0", "MAR_1of3", "MAR_2of3", "Present_3")
    )
  ) %>%
  filter(!is.na(detect_class))

# ----------------------------
# 3) "Global abundance" per protein × condition × channel
#    Use mean log2 intensity across replicates (only where detected/log2 finite)
# ----------------------------
abundance_df <- det_long %>%
  group_by(channel, Protein.Group, condition) %>%
  summarise(
    mean_log2_intensity = mean(log2_intensity, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- det_by_cond %>%
  left_join(abundance_df, by = c("channel", "Protein.Group", "condition")) %>%
  filter(!is.na(mean_log2_intensity))

# ----------------------------
# 4) Boxplot: distribution of mean log2 intensity by missingness class
# ----------------------------
p_box <- ggplot(plot_df, aes(x = detect_class, y = mean_log2_intensity, fill = detect_class)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ channel, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Global intensity distribution by detection class",
    x = "Detection class (within condition, 3 replicates)",
    y = "Mean log2(Intensity) across replicates"
  )

print(p_box)

# ----------------------------
# 5) Density plot: distribution of mean log2 intensity by class
# ----------------------------
p_den <- ggplot(plot_df, aes(x = mean_log2_intensity, fill = detect_class)) +
  geom_density(alpha = 0.4, na.rm = TRUE) +
  facet_wrap(~ channel, scales = "free") +
  theme_bw(base_size = 12) +
  labs(
    title = "Mean log2(Intensity) density by detection class",
    x = "Mean log2(Intensity) across replicates",
    y = "Density",
    fill = "Detection class"
  )

print(p_den)

p_den <- ggplot(plot_df, aes(x = mean_log2_intensity, fill = detect_class)) +
  geom_density(alpha = 0.4, na.rm = TRUE) +
  facet_wrap(~ channel+condition, ncol=4) +
  theme_bw(base_size = 12) +
  labs(
    title = "Mean log2(Intensity) density by detection class",
    x = "Mean log2(Intensity) across replicates",
    y = "Density",
    fill = "Detection class"
  )

print(p_den)
