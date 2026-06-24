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
out_separated_long_csv <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_separated_classes_long.csv")
out_separated_wide_csv <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_separated_classes_wide.csv")
out_separated_main_pdf <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_separated_classes_main_text.pdf")
out_separated_main_png <- file.path(out_dir, "Recovery_distinct_precursors_by_lipid_class_separated_classes_main_text.png")
out_separated_composition_pdf <- file.path(out_dir, "Recovery_lipid_class_precursor_composition_separated_classes_stacked_bar.pdf")
out_separated_composition_png <- file.path(out_dir, "Recovery_lipid_class_precursor_composition_separated_classes_stacked_bar.png")

# ---------------------------------------------------------
# 2. HELPERS
# ---------------------------------------------------------
prepare_detected_precursors <- function(df, group, sample_cols, blank_cutoff = 1.3) {
  missing_cols <- setdiff(c(sample_cols, "blank", "precursor", "lipid_class", "lipid_name"), colnames(df))
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
      LipidName = as.character(lipid_name),
      precursor = precursor
    ) |>
    distinct(Group, Genotype, Treatment, LipidClass, LipidName, precursor)
}

get_first_annotation <- function(x) {
  sub("_.*$", "", x)
}

# The manuscript-facing view intentionally combines PC/LPC/ether-PC and
# PE/LPE/ether-PE into broader PC and PE pools.
collapse_lipid_class_main <- function(lipid_class, lipid_name) {
  first_annotation <- get_first_annotation(lipid_name)
  
  case_when(
    lipid_class %in% c("PC", "LPC") | grepl("^PC O-|^PC P-|^PC\\(O-|^PC\\(P-", first_annotation) ~ "PC",
    lipid_class %in% c("PE", "LPE") | grepl("^PE O-|^PE P-|^PE\\(O-|^PE\\(P-", first_annotation) ~ "PE",
    TRUE ~ lipid_class
  )
}

# The backup view preserves detail only when the first annotation supports it.
# Alternate candidates after underscores are intentionally not used for splitting.
collapse_lipid_class_separated <- function(lipid_class, lipid_name) {
  first_annotation <- get_first_annotation(lipid_name)
  
  case_when(
    grepl("^PC O-|^PC\\(O-", first_annotation) ~ "PC Ether",
    grepl("^PC P-|^PC\\(P-", first_annotation) ~ "PC Plasmalogen",
    grepl("^PE O-|^PE\\(O-", first_annotation) ~ "PE Ether",
    grepl("^PE P-|^PE\\(P-", first_annotation) ~ "PE Plasmalogen",
    TRUE ~ lipid_class
  )
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

build_count_tables <- function(detected_df, class_mapper) {
  count_long <- detected_df |>
    mutate(LipidClass = class_mapper(LipidClass, LipidName)) |>
    distinct(Group, Genotype, Treatment, LipidClass, precursor) |>
    count(Group, Genotype, Treatment, LipidClass, name = "Count") |>
    complete(
      Group = names(group_sample_cols),
      LipidClass = sort(unique(class_mapper(detected_df$LipidClass, detected_df$LipidName))),
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
  
  list(
    long = count_long,
    wide = count_wide,
    class_order = class_order
  )
}

main_counts <- build_count_tables(detected_precursors, collapse_lipid_class_main)
separated_counts <- build_count_tables(detected_precursors, collapse_lipid_class_separated)

count_long <- main_counts$long
count_wide <- main_counts$wide
class_order <- main_counts$class_order

write.csv(count_long, out_long_csv, row.names = FALSE)
write.csv(count_wide, out_wide_csv, row.names = FALSE)
write.csv(separated_counts$long, out_separated_long_csv, row.names = FALSE)
write.csv(separated_counts$wide, out_separated_wide_csv, row.names = FALSE)

# ---------------------------------------------------------
# 4. FIGURES
# ---------------------------------------------------------
treatment_colors <- c("Veh" = "#4D4D4D", "100 uM" = "#E69F00", "200 uM" = "#0072B2")

make_main_plot <- function(plot_df, title, subtitle) {
  ggplot(plot_df, aes(x = Count, y = LipidClass, color = Treatment)) +
    geom_line(aes(group = interaction(Genotype, LipidClass)), color = "grey82", linewidth = 0.4) +
    geom_point(size = 2.4, alpha = 0.95) +
    facet_wrap(~ Genotype, nrow = 1) +
    scale_color_manual(values = treatment_colors) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 14),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "grey25"),
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
  title = "Detected recovery lipid precursor classes",
  subtitle = paste0("Distinct precursors after filtering; PC/PE subclasses are grouped (", blank_cutoff, "x blank)")
)

p_main_separated <- make_main_plot(
  separated_counts$long,
  title = "Detected recovery lipid precursor classes",
  subtitle = paste0("Backup view with first-annotation ether/plasmalogen subclasses retained (", blank_cutoff, "x blank)")
)

make_composition_plot <- function(plot_df, title, subtitle) {
  label_df <- plot_df |>
    distinct(Genotype, Treatment, TotalPrecursors)
  
  ggplot(plot_df, aes(x = Treatment, y = Count, fill = LipidClass)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.18) +
    geom_text(
      data = label_df,
      aes(x = Treatment, y = TotalPrecursors + 5, label = paste0("n=", TotalPrecursors)),
      inherit.aes = FALSE,
      size = 4.2,
      color = "grey20"
    ) +
    facet_wrap(~ Genotype, nrow = 1) +
    paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 18),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 8.5, color = "grey25")
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
  title = "Recovery lipid class precursor composition",
  subtitle = "Stack height shows total distinct precursors; PC/PE subclasses are grouped."
)

p_composition_separated <- make_composition_plot(
  separated_counts$long,
  title = "Recovery lipid class precursor composition",
  subtitle = "Backup view with first-annotation ether/plasmalogen subclasses retained."
)

ggsave(out_main_pdf, p_main, width = 7.2, height = 4.8, useDingbats = FALSE)
ggsave(out_main_png, p_main, width = 7.2, height = 4.8, dpi = 900)
ggsave(out_composition_pdf, p_composition, width = 7.2, height = 4.8, useDingbats = FALSE)
ggsave(out_composition_png, p_composition, width = 7.2, height = 4.8, dpi = 600)
ggsave(out_separated_main_pdf, p_main_separated, width = 7.2, height = 5.2, useDingbats = FALSE)
ggsave(out_separated_main_png, p_main_separated, width = 7.2, height = 5.2, dpi = 900)
ggsave(out_separated_composition_pdf, p_composition_separated, width = 8.0, height = 5.2, useDingbats = FALSE)
ggsave(out_separated_composition_png, p_composition_separated, width = 8.0, height = 5.2, dpi = 600)

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
