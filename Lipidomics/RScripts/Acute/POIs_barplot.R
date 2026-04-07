library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(janitor)

# ---------------------------------------------------------
# 1. INPUTS
# ---------------------------------------------------------

df_long <- read.csv("Lipidomics/output_txts/Acute_filteredlipidtable.csv", check.names = FALSE)
meta <- read.csv("Lipidomics/metadata/Sample_description guide_5hratral.csv", check.names = FALSE)

# clean metadata column names
meta <- meta %>%
  janitor::clean_names()

# check metadata names
print(colnames(meta))

# if your metadata has 'treament' misspelled, fix it
if ("treament" %in% colnames(meta)) {
  meta <- meta %>% rename(treatment = treament)
}

# if sample column is not lower-case after import, fix that too
if (!"sample" %in% colnames(meta)) {
  stop("Metadata file must contain a sample column.")
}
if (!"treatment" %in% colnames(meta)) {
  stop("Metadata file must contain a treatment column.")
}

colnames(df_long)
df_long <- df_long[, colnames(df_long) != ""]
df_long <- df_long |> select(-blank)

# lipid to plot
lipid_to_plot <- "PE 34:1"
mrm_to_plot <- NULL

plot_points <- TRUE
use_sd <- TRUE

treatment_order <- c("Veh", "100", "200")
treatment_labels <- c(
  "Veh" = "Veh",
  "100" = "100uM",
  "200" = "200uM"
)

treatment_fill <- c(
  "Veh" = "#BDBDBD",
  "100" = "#E08A5F",
  "200" = "#E6CDB9"
)

out_file_bar <- paste0("Lipid_barplot_", gsub("[^A-Za-z0-9]+", "_", lipid_to_plot), ".svg")
out_file_box <- paste0("Lipid_boxplot_", gsub("[^A-Za-z0-9]+", "_", lipid_to_plot), ".svg")

# ---------------------------------------------------------
# 2. RESHAPE LIPID TABLE TO LONG FORMAT
# ---------------------------------------------------------

sample_cols <- grep("^s\\d+$", colnames(df_long), value = TRUE)

lipid_long <- df_long %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "intensity"
  ) %>%
  left_join(meta, by = "sample") %>%
  mutate(
    Treatment = factor(as.character(treatment), levels = treatment_order)
  )

# quick check
print(colnames(lipid_long))
print(unique(lipid_long$Treatment))

# ---------------------------------------------------------
# 3. FILTER TO LIPID OF INTEREST
# ---------------------------------------------------------

lipid_df <- lipid_long %>%
  filter(
    lipid_name == lipid_to_plot,
    if (!is.null(mrm_to_plot)) mrm == mrm_to_plot else TRUE
  )

if (nrow(lipid_df) == 0) {
  stop("No rows found for the requested lipid. Check 'lipid_to_plot' and optional 'mrm_to_plot'.")
}

lipid_df %>%
  distinct(lipid_name, mrm, lipid_class, lipid_class1, precursor, product) %>%
  print(n = 50)

# ---------------------------------------------------------
# 4. SUMMARIZE BY TREATMENT
# ---------------------------------------------------------

summary_df <- lipid_df %>%
  filter(!is.na(Treatment), !is.na(intensity)) %>%
  group_by(Treatment) %>%
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
# 5. BAR PLOT: mean +/- SD/SEM with replicate points
# ---------------------------------------------------------

p_bar <- ggplot(summary_df, aes(x = Treatment, y = log2(mean_intensity), fill = Treatment)) +
  geom_col(color = "black", width = 0.7) +
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
        aes(x = Treatment, y = intensity),
        inherit.aes = FALSE,
        size = 2.5,
        position = position_jitter(width = 0.08, height = 0)
      )
    }
  } +
  scale_x_discrete(labels = treatment_labels) +
  scale_fill_manual(values = treatment_fill, guide = "none") +
  labs(
    title = paste(lipid_to_plot, "abundance across acute treatment"),
    subtitle = if (!is.null(mrm_to_plot)) paste("MRM:", mrm_to_plot) else NULL,
    x = "Treatment",
    y = "Intensity"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar)

# ---------------------------------------------------------
# 6. BOXPLOT + POINTS
# ---------------------------------------------------------

p_box <- ggplot(lipid_df, aes(x = Treatment, y = intensity, fill = Treatment)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black") +
  geom_point(
    size = 2.5,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  scale_x_discrete(labels = treatment_labels) +
  scale_fill_manual(values = treatment_fill, guide = "none") +
  labs(
    title = paste(lipid_to_plot, "distribution across acute treatment"),
    subtitle = if (!is.null(mrm_to_plot)) paste("MRM:", mrm_to_plot) else NULL,
    x = "Treatment",
    y = "Intensity"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_box)

# ---------------------------------------------------------
# 7. OPTIONAL: LOG2 INTENSITY VERSION
# ---------------------------------------------------------

lipid_df_log2 <- lipid_df %>%
  mutate(log2_intensity = log2(intensity + 1))

summary_df_log2 <- lipid_df_log2 %>%
  filter(!is.na(Treatment), !is.na(log2_intensity)) %>%
  group_by(Treatment) %>%
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

p_bar_log2 <- ggplot(summary_df_log2, aes(x = Treatment, y = mean_log2, fill = Treatment,)) +
  geom_point(color = "red", width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_log2 - error_value, ymax = mean_log2 + error_value),
    width = 0.2,
    linewidth = 0.6
  ) +
  {
    if (plot_points) {
      geom_point(
        data = lipid_df_log2,
        aes(x = Treatment, y = log2_intensity),
        inherit.aes = FALSE,
        size = 2.5,
        position = position_jitter(width = 0.08, height = 0)
      )
    }
  } +
  scale_x_discrete(labels = treatment_labels) +
  scale_fill_manual(values = treatment_fill, guide = "none") +
  labs(
    title = paste(lipid_to_plot, "log2 abundance across acute treatment"),
    subtitle = if (!is.null(mrm_to_plot)) paste("MRM:", mrm_to_plot) else NULL,
    x = "Treatment",
    y = "log2(Intensity)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_bar_log2)

# ---------------------------------------------------------
# 8. SAVE
# ---------------------------------------------------------

ggsave(out_file_bar, p_bar_log2, width = 8, height = 6)
ggsave(out_file_box, p_box, width = 8, height = 6)