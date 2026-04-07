library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

# ---------------------------------------------------------
# 1. INPUTS
# ---------------------------------------------------------

df_long <- read.csv("Lipidomics/output_txts/Acute_DEA_long.csv", check.names = FALSE)

# remove unnamed first column if present
df_long <- df_long[, colnames(df_long) != ""]

# lipid to plot
lipid_to_plot <- "DG 34:1 NL 18:1"

# optional: filter to a specific MRM
mrm_to_plot <- NULL

# plot settings
plot_points <- TRUE
use_sd <- TRUE   # TRUE = SD, FALSE = SEM

condition_order <- c("RDH12_Veh", "RDH12_100", "RDH12_200")

condition_labels <- c(
  "RDH12_Veh" = "Veh",
  "RDH12_100" = "100uM",
  "RDH12_200" = "200uM"
)

condition_fill <- c(
  "RDH12_Veh" = "#BDBDBD",
  "RDH12_100" = "#E08A5F",
  "RDH12_200" = "#E6CDB9"
)

out_file_bar_raw  <- paste0("Lipidomics/Figures/Acute/Lipid_barplot_raw_",  gsub("[^A-Za-z0-9]+", "_", lipid_to_plot), ".svg")
out_file_box_raw  <- paste0("Lipidomics/Figures/Acute/Lipid_boxplot_raw_",  gsub("[^A-Za-z0-9]+", "_", lipid_to_plot), ".svg")
out_file_bar_log2 <- paste0("Lipidomics/Figures/Acute/Lipid_barplot_log2_", gsub("[^A-Za-z0-9]+", "_", lipid_to_plot), ".svg")
out_file_box_log2 <- paste0("Lipidomics/Figures/Acute/Lipid_boxplot_log2_", gsub("[^A-Za-z0-9]+", "_", lipid_to_plot), ".svg")

# ---------------------------------------------------------
# 2. CHECK INPUT TABLE
# ---------------------------------------------------------

print(colnames(df_long))
print(unique(df_long$condition))

# ---------------------------------------------------------
# 3. FILTER TO LIPID OF INTEREST
# ---------------------------------------------------------

lipid_df <- df_long %>%
  mutate(
    condition = factor(condition, levels = condition_order)
  ) %>%
  filter(
    name == lipid_to_plot,
    if (!is.null(mrm_to_plot)) mrm == mrm_to_plot else TRUE
  )

if (nrow(lipid_df) == 0) {
  stop("No rows found for the requested lipid. Check 'lipid_to_plot' and optional 'mrm_to_plot'.")
}

# inspect matching rows
lipid_df %>%
  distinct(name, mrm, ID, picked_candidate) 
 # print(n = 50)

# ---------------------------------------------------------
# 4. SUMMARIZE BY CONDITION (RAW INTENSITY)
# ---------------------------------------------------------

summary_df <- lipid_df %>%
  filter(!is.na(condition), !is.na(intensity)) %>%
  group_by(condition) %>%
  summarise(
    mean_intensity = mean(intensity, na.rm = TRUE),
    sd_intensity   = sd(intensity, na.rm = TRUE),
    n              = sum(!is.na(intensity)),
    sem_intensity  = sd_intensity / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    error_value = if (use_sd) sd_intensity else sem_intensity
  )

print(summary_df)

# ---------------------------------------------------------
# 5. BAR PLOT: RAW INTENSITY
# ---------------------------------------------------------

p_bar_raw <- ggplot(summary_df, aes(x = condition, y = mean_intensity, fill = condition)) +
  geom_point(color = "red", width = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_intensity - error_value,
      ymax = mean_intensity + error_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  {
    if (plot_points) {
      geom_point(
        data = lipid_df,
        aes(x = condition, y = intensity),
        inherit.aes = FALSE,
        size = 2.5,
        position = position_jitter(width = 0.08, height = 0)
      )
    }
  } +
  scale_x_discrete(labels = condition_labels) +
  scale_fill_manual(values = condition_fill, guide = "none") +
  labs(
    title = paste(lipid_to_plot, "abundance across acute treatment"),
    subtitle = if (!is.null(mrm_to_plot)) paste("MRM:", mrm_to_plot) else NULL,
    x = "Condition",
    y = "Intensity"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar_raw)

# ---------------------------------------------------------
# 6. BOXPLOT + POINTS: RAW INTENSITY
# ---------------------------------------------------------

p_box_raw <- ggplot(lipid_df, aes(x = condition, y = intensity, fill = condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black") +
  geom_point(
    size = 2.5,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  scale_x_discrete(labels = condition_labels) +
  scale_fill_manual(values = condition_fill, guide = "none") +
  labs(
    title = paste(lipid_to_plot, "distribution across acute treatment"),
    subtitle = if (!is.null(mrm_to_plot)) paste("MRM:", mrm_to_plot) else NULL,
    x = "Condition",
    y = "Intensity"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_box_raw)

# ---------------------------------------------------------
# 7. LOG2 INTENSITY VERSION
# ---------------------------------------------------------

lipid_df_log2 <- lipid_df %>%
  mutate(log2_intensity = log2(intensity + 1))

summary_df_log2 <- lipid_df_log2 %>%
  filter(!is.na(condition), !is.na(log2_intensity)) %>%
  group_by(condition) %>%
  summarise(
    mean_log2 = mean(log2_intensity, na.rm = TRUE),
    sd_log2   = sd(log2_intensity, na.rm = TRUE),
    n         = sum(!is.na(log2_intensity)),
    sem_log2  = sd_log2 / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    error_value = if (use_sd) sd_log2 else sem_log2
  )

print(summary_df_log2)

# ---------------------------------------------------------
# 8. BAR PLOT: LOG2 INTENSITY
# ---------------------------------------------------------

p_bar_log2 <- ggplot(summary_df_log2, aes(x = condition, y = mean_log2, fill = condition)) +
  geom_point(color = "red", width = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_log2 - error_value,
      ymax = mean_log2 + error_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  {
    if (plot_points) {
      geom_point(
        data = lipid_df_log2,
        aes(x = condition, y = log2_intensity),
        inherit.aes = FALSE,
        size = 2.5,
        position = position_jitter(width = 0.08, height = 0)
      )
    }
  } +
  scale_x_discrete(labels = condition_labels) +
  scale_fill_manual(values = condition_fill, guide = "none") +
  labs(
    title = paste(lipid_to_plot, "log2 abundance across acute treatment"),
    subtitle = if (!is.null(mrm_to_plot)) paste("MRM:", mrm_to_plot) else NULL,
    x = "Condition",
    y = "log2(Intensity)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar_log2)

# ---------------------------------------------------------
# 9. BOXPLOT + POINTS: LOG2 INTENSITY
# ---------------------------------------------------------

p_box_log2 <- ggplot(lipid_df_log2, aes(x = condition, y = log2_intensity, fill = condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black") +
  geom_point(
    size = 2.5,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  scale_x_discrete(labels = condition_labels) +
  scale_fill_manual(values = condition_fill, guide = "none") +
  labs(
    title = paste(lipid_to_plot, "log2 distribution across acute treatment"),
    subtitle = if (!is.null(mrm_to_plot)) paste("MRM:", mrm_to_plot) else NULL,
    x = "Condition",
    y = "log2(Intensity + 1)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_box_log2)

# ---------------------------------------------------------
# 10. OPTIONAL: PRINT DEA STATS FOR THIS LIPID
# ---------------------------------------------------------

lipid_stats <- lipid_df %>%
  select(
    name, mrm, ID,
    starts_with("RDH12_100_vs_RDH12_Veh"),
    starts_with("RDH12_200_vs_RDH12_100"),
    starts_with("RDH12_200_vs_RDH12_Veh")
  ) %>%
  distinct()

#print(lipid_stats, n = 50)

# ---------------------------------------------------------
# 11. SAVE
# ---------------------------------------------------------

ggsave(out_file_bar_raw,  p_bar_raw,  width = 8, height = 6)
ggsave(out_file_box_raw,  p_box_raw,  width = 8, height = 6)
ggsave(out_file_bar_log2, p_bar_log2, width = 8, height = 6)
ggsave(out_file_box_log2, p_box_log2, width = 8, height = 6)