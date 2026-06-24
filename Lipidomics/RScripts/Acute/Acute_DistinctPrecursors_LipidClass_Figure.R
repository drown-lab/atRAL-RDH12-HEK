# =========================================================
# Acute lipidomics: distinct precursor counts per lipid class
# - Reads the per-treatment class summary tables
# - Collapses ether labels to match the manuscript summary table
# - Saves a main-text-ready figure and the plotting tables
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
input_files <- c(
  "Veh" = "Lipidomics/output_txts/AcuteVeh_LipidClassEther_summary.csv",
  "100" = "Lipidomics/output_txts/Acute100_LipidClassEther_summary.csv",
  "200" = "Lipidomics/output_txts/Acute200_LipidClassEther_summary.csv"
)

out_dir <- "Lipidomics/Figures/Acute/Distinct_precursors_by_class/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_long_csv <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_long.csv")
out_wide_csv <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_wide.csv")
out_main_pdf <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_main_text.pdf")
out_main_png <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_main_text.png")
out_composition_pdf <- file.path(out_dir, "Acute_lipid_class_precursor_composition_stacked_bar.pdf")
out_composition_png <- file.path(out_dir, "Acute_lipid_class_precursor_composition_stacked_bar.png")
out_separated_long_csv <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_separated_classes_long.csv")
out_separated_wide_csv <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_separated_classes_wide.csv")
out_separated_main_pdf <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_separated_classes_main_text.pdf")
out_separated_main_png <- file.path(out_dir, "Acute_distinct_precursors_by_lipid_class_separated_classes_main_text.png")
out_separated_composition_pdf <- file.path(out_dir, "Acute_lipid_class_precursor_composition_separated_classes_stacked_bar.pdf")
out_separated_composition_png <- file.path(out_dir, "Acute_lipid_class_precursor_composition_separated_classes_stacked_bar.png")

# ---------------------------------------------------------
# 2. READ AND STANDARDIZE CLASS COUNTS
# ---------------------------------------------------------
read_class_summary <- function(path, treatment) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }
  
  read.csv(path, check.names = TRUE, stringsAsFactors = FALSE) %>%
    transmute(
      Treatment = treatment,
      LipidClassRaw = lipid_class_ether,
      Count = distinct_precursor
    )
}

# The manuscript-facing view groups phosphocholine and phosphoethanolamine
# subclasses into broader PC and PE pools to avoid overinterpreting annotation
# details in the compact stacked bar figure.
collapse_lipid_class_main <- function(x) {
  case_when(
    x %in% c("PC", "PC O-", "PC P-", "PC Ether", "PC Plasmalogen", "LPC", "LPC O-") ~ "PC",
    x %in% c("PE", "PE O-", "PE P-", "PE Ether", "PE Plasmalogen", "LPE", "LPE O-") ~ "PE",
    x == "DG O-" ~ "DG",
    x == "PI O-" ~ "PI",
    TRUE ~ x
  )
}

# The backup view preserves ether/plasmalogen/LPC/LPE detail where the summary
# table has it, while keeping the older DG/PI cleanup behavior.
clean_lipid_class_separated <- function(x) {
  case_when(
    x == "PC O-" ~ "PC Ether",
    x == "PC P-" ~ "PC Plasmalogen",
    x == "PE O-" ~ "PE Ether",
    x == "PE P-" ~ "PE Plasmalogen",
    x == "DG O-" ~ "DG",
    x == "PI O-" ~ "PI",
    TRUE ~ x
  )
}

raw_count_long <- bind_rows(
  lapply(names(input_files), function(treatment) {
    read_class_summary(input_files[[treatment]], treatment)
  })
)

build_count_tables <- function(raw_df, class_mapper) {
  count_long <- raw_df %>%
    mutate(
      LipidClass = class_mapper(LipidClassRaw),
      Treatment = factor(Treatment, levels = c("Veh", "100", "200"))
    ) %>%
    group_by(Treatment, LipidClass) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  class_order <- count_long %>%
    group_by(LipidClass) %>%
    summarise(MaxCount = max(Count), .groups = "drop") %>%
    arrange(desc(MaxCount), LipidClass) %>%
    pull(LipidClass)
  
  count_long <- count_long %>%
    complete(
      Treatment = factor(c("Veh", "100", "200"), levels = c("Veh", "100", "200")),
      LipidClass = class_order,
      fill = list(Count = 0)
    ) %>%
    group_by(Treatment) %>%
    mutate(
      TotalPrecursors = sum(Count),
      Proportion = Count / TotalPrecursors
    ) %>%
    ungroup() %>%
    mutate(
      LipidClass = factor(LipidClass, levels = rev(class_order)),
      TreatmentLabel = factor(
        recode(as.character(Treatment), "100" = "100 uM", "200" = "200 uM"),
        levels = c("Veh", "100 uM", "200 uM")
      )
    )
  
  count_wide_classes <- count_long %>%
    select(LipidClass, Treatment, Count) %>%
    pivot_wider(names_from = Treatment, values_from = Count) %>%
    arrange(desc(Veh), LipidClass) %>%
    mutate(LipidClass = as.character(LipidClass))
  
  count_wide <- bind_rows(
    count_wide_classes,
    count_wide_classes %>%
      summarise(
        LipidClass = "Sum",
        Veh = sum(Veh),
        `100` = sum(`100`),
        `200` = sum(`200`)
      )
  )
  
  list(
    long = count_long,
    wide = count_wide,
    class_order = class_order
  )
}

main_counts <- build_count_tables(raw_count_long, collapse_lipid_class_main)
separated_counts <- build_count_tables(raw_count_long, clean_lipid_class_separated)

count_long <- main_counts$long
count_wide <- main_counts$wide
class_order <- main_counts$class_order

write.csv(count_long, out_long_csv, row.names = FALSE)
write.csv(count_wide, out_wide_csv, row.names = FALSE)
write.csv(separated_counts$long, out_separated_long_csv, row.names = FALSE)
write.csv(separated_counts$wide, out_separated_wide_csv, row.names = FALSE)

# ---------------------------------------------------------
# 3. MAIN TEXT FIGURE: HORIZONTAL DOT PLOT
# ---------------------------------------------------------
treatment_colors <- c(
  "Veh" = "#4D4D4D",
  "100 uM" = "#E69F00",
  "200 uM" = "#0072B2"
)

make_main_plot <- function(plot_df, title, subtitle) {
  ggplot(plot_df, aes(x = Count, y = LipidClass, color = TreatmentLabel)) +
    geom_line(
      aes(group = LipidClass),
      color = "grey82",
      linewidth = 0.45
    ) +
    geom_point(size = 2.8, alpha = 0.95) +
    scale_color_manual(values = treatment_colors) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "grey25"),
      plot.margin = margin(6, 12, 6, 6)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Number of distinct precursors",
      y = NULL
    )
}

p_main <- make_main_plot(
  count_long,
  title = "Detected lipid precursor classes",
  subtitle = "PC-related and PE-related subclasses are grouped for the manuscript view"
)

p_main_separated <- make_main_plot(
  separated_counts$long,
  title = "Detected lipid precursor classes",
  subtitle = "Backup view with ether/plasmalogen subclasses retained where available"
)

ggsave(out_main_pdf, p_main, width = 5.8, height = 4.5, useDingbats = FALSE)
ggsave(out_main_png, p_main, width = 5.8, height = 4.5, dpi = 900)

# ---------------------------------------------------------
# 4. COMPOSITION BAR PLOTS
# ---------------------------------------------------------
make_composition_plot <- function(plot_df, class_order, title, subtitle) {
  composition_df <- plot_df %>%
    mutate(LipidClass = fct_relevel(LipidClass, rev(class_order)))
  
  ggplot(composition_df, aes(x = TreatmentLabel, y = Count, fill = LipidClass)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.2) +
    geom_text(
      data = plot_df %>% distinct(TreatmentLabel, TotalPrecursors),
      aes(x = TreatmentLabel, y = TotalPrecursors + 5, label = paste0("n=", TotalPrecursors)),
      inherit.aes = FALSE,
      size = 4.5,
      color = "grey20"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "grey25")
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Number of distinct precursors"
    )
}

p_composition <- make_composition_plot(
  count_long,
  class_order,
  title = "Lipid precursor class composition",
  subtitle = "Stack height shows total distinct precursors; PC/PE subclasses are grouped"
)

p_composition_separated <- make_composition_plot(
  separated_counts$long,
  separated_counts$class_order,
  title = "Lipid precursor class composition",
  subtitle = "Backup view with ether/plasmalogen subclasses retained where available"
)

ggsave(out_composition_pdf, p_composition, width = 5.8, height = 4.5, useDingbats = FALSE)
ggsave(out_composition_png, p_composition, width = 5.8, height = 4.5, dpi = 600)
ggsave(out_separated_main_pdf, p_main_separated, width = 5.8, height = 4.8, useDingbats = FALSE)
ggsave(out_separated_main_png, p_main_separated, width = 5.8, height = 4.8, dpi = 900)
ggsave(out_separated_composition_pdf, p_composition_separated, width = 6.2, height = 4.8, useDingbats = FALSE)
ggsave(out_separated_composition_png, p_composition_separated, width = 6.2, height = 4.8, dpi = 600)

# ---------------------------------------------------------
# 5. CONSOLE SUMMARY
# ---------------------------------------------------------
print(count_wide)

cat("\nSaved files:\n")
cat(out_long_csv, "\n")
cat(out_wide_csv, "\n")
cat(out_main_pdf, "\n")
cat(out_main_png, "\n")
cat(out_composition_pdf, "\n")
cat(out_composition_png, "\n")
cat(out_separated_long_csv, "\n")
cat(out_separated_wide_csv, "\n")
cat(out_separated_main_pdf, "\n")
cat(out_separated_main_png, "\n")
cat(out_separated_composition_pdf, "\n")
cat(out_separated_composition_png, "\n")
