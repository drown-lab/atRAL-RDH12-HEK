# =========================================================
# Recovery lipidomics: distinct precursor counts per lipid class
# - Applies treatment-specific signal-to-blank detection filtering
# - Generates compact stacked composition bars split by genotype
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(paletteer)
  library(tidyr)
})

# ---------------------------------------------------------
# 1. INPUTS
# ---------------------------------------------------------
recovery_file <- "Lipidomics/output_txts/Recovery_filteredlipidtable.csv"
out_dir <- "Lipidomics/Figures/Recovery/Distinct_precursors_by_class/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

blank_cutoff <- 1.3

group_sample_cols <- list(
  "Control_Veh" = c("s1", "s2", "s3"),
  "Control_100" = c("s4", "s5", "s6"),
  "Control_200" = c("s7", "s8", "s9"),
  "RDH12_Veh" = c("s10", "s11", "s12"),
  "RDH12_100" = c("s13", "s14", "s15"),
  "RDH12_200" = c("s16", "s17", "s18")
)

out_long_csv <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_long.csv")
out_wide_csv <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_wide.csv")
out_main_pdf <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_main_text.pdf")
out_main_png <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_main_text.png")
out_composition_pdf <- file.path(out_dir, "Recovery_lipid_class_precursor_composition_stacked_bar.pdf")
out_composition_png <- file.path(out_dir, "Recovery_lipid_class_precursor_composition_stacked_bar.png")

# ---------------------------------------------------------
# 2. HELPERS
# ---------------------------------------------------------
prepare_detected_precursors <- function(df, group, sample_cols, blank_cutoff = 1.3) {
  missing_cols <- setdiff(c(sample_cols, "blank", "precursor", "lipid_class"), colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns for ", group, ":\n", paste(missing_cols, collapse = "\n"))
  }
  
  df |>
    rowwise() |>
    mutate(
      max_group_signal = max(c_across(all_of(sample_cols)), na.rm = TRUE),
      signal_to_blank = max_group_signal / blank
    ) |>
    ungroup() |>
    filter(
      is.finite(signal_to_blank),
      signal_to_blank >= blank_cutoff,
      !is.na(precursor),
      !is.na(lipid_class)
    ) |>
    transmute(
      Group = group,
      Genotype = if_else(grepl("^Control", group), "WT", "RDH12"),
      Treatment = case_when(
        grepl("_Veh$", group) ~ "Veh",
        grepl("_100$", group) ~ "100 uM",
        grepl("_200$", group) ~ "200 uM",
        TRUE ~ group
      ),
      LipidClass = as.character(lipid_class),
      precursor = precursor
    ) |>
    distinct(Group, Genotype, Treatment, LipidClass, precursor)
}

# ---------------------------------------------------------
# 3. BUILD TABLES
# ---------------------------------------------------------
lipid_df <- read.csv(recovery_file, check.names = TRUE, stringsAsFactors = FALSE)

detected_precursors <- bind_rows(
  lapply(names(group_sample_cols), function(group) {
    prepare_detected_precursors(lipid_df, group, group_sample_cols[[group]], blank_cutoff)
  })
)

count_long <- detected_precursors |>
  count(Group, Genotype, Treatment, LipidClass, name = "Count") |>
  complete(
    Group = names(group_sample_cols),
    LipidClass = sort(unique(detected_precursors$LipidClass)),
    fill = list(Count = 0)
  ) |>
  mutate(
    Genotype = if_else(grepl("^Control", Group), "WT", "RDH12"),
    Treatment = case_when(
      grepl("_Veh$", Group) ~ "Veh",
      grepl("_100$", Group) ~ "100 uM",
      grepl("_200$", Group) ~ "200 uM",
      TRUE ~ Group
    ),
    Genotype = factor(Genotype, levels = c("WT", "RDH12")),
    Treatment = factor(Treatment, levels = c("Veh", "100 uM", "200 uM"))
  ) |>
  group_by(Group) |>
  mutate(
    TotalPrecursors = sum(Count),
    Proportion = if_else(TotalPrecursors > 0, Count / TotalPrecursors, 0)
  ) |>
  ungroup()

class_order <- count_long |>
  group_by(LipidClass) |>
  summarise(MaxCount = max(Count), .groups = "drop") |>
  arrange(desc(MaxCount), LipidClass) |>
  pull(LipidClass)

count_long <- count_long |>
  mutate(LipidClass = factor(LipidClass, levels = rev(class_order)))

count_wide <- count_long |>
  select(LipidClass, Group, Count) |>
  mutate(LipidClass = as.character(LipidClass)) |>
  pivot_wider(names_from = Group, values_from = Count) |>
  arrange(desc(Control_Veh), LipidClass)

count_wide <- bind_rows(
  count_wide,
  count_wide |>
    summarise(
      LipidClass = "Sum",
      Control_Veh = sum(Control_Veh),
      Control_100 = sum(Control_100),
      Control_200 = sum(Control_200),
      RDH12_Veh = sum(RDH12_Veh),
      RDH12_100 = sum(RDH12_100),
      RDH12_200 = sum(RDH12_200)
    )
)

write.csv(count_long, out_long_csv, row.names = FALSE)
write.csv(count_wide, out_wide_csv, row.names = FALSE)

# ---------------------------------------------------------
# 4. FIGURES
# ---------------------------------------------------------
treatment_colors <- c("Veh" = "#4D4D4D", "100 uM" = "#E69F00", "200 uM" = "#0072B2")

p_main <- ggplot(count_long, aes(x = Count, y = LipidClass, color = Treatment)) +
  geom_line(aes(group = interaction(Genotype, LipidClass)), color = "grey82", linewidth = 0.4) +
  geom_point(size = 2.4, alpha = 0.95) +
  facet_wrap(~ Genotype, nrow = 1) +
  scale_color_manual(values = treatment_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(color = "grey20"),
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey25"),
    plot.margin = margin(6, 12, 6, 6)
  ) +
  labs(
    title = "Detected recovery lipid precursor classes",
    subtitle = paste0("Distinct precursors after treatment-specific detection filtering (", blank_cutoff, "x blank)"),
    x = "Number of distinct precursors",
    y = NULL
  )

label_df <- count_long |>
  distinct(Genotype, Treatment, TotalPrecursors)

p_composition <- ggplot(count_long, aes(x = Treatment, y = Count, fill = LipidClass)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.18) +
  geom_text(
    data = label_df,
    aes(x = Treatment, y = TotalPrecursors + 5, label = paste0("n=", TotalPrecursors)),
    inherit.aes = FALSE,
    size = 3.1,
    color = "grey20"
  ) +
  facet_wrap(~ Genotype, nrow = 1) +
  paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "grey20"),
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey25")
  ) +
  labs(
    title = "Recovery lipid class precursor composition",
    subtitle = "Stack height shows total distinct precursors detected per genotype/treatment.",
    x = NULL,
    y = "Number of distinct precursors"
  )

ggsave(out_main_pdf, p_main, width = 7.2, height = 4.8, useDingbats = FALSE)
ggsave(out_main_png, p_main, width = 7.2, height = 4.8, dpi = 900)
ggsave(out_composition_pdf, p_composition, width = 7.2, height = 4.8, useDingbats = FALSE)
ggsave(out_composition_png, p_composition, width = 7.2, height = 4.8, dpi = 600)

cat("\nSaved files:\n")
cat(out_long_csv, "\n")
cat(out_wide_csv, "\n")
cat(out_main_pdf, "\n")
cat(out_main_png, "\n")
cat(out_composition_pdf, "\n")
cat(out_composition_png, "\n")
