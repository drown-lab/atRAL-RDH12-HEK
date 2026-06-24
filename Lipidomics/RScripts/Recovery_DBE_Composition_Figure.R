# =========================================================
# Recovery lipidomics DBE composition figures
# - Counts distinct detected precursors by lipid class and DBE group
# - Main heatmap focuses on treatment-vs-EtOH within genotype
# - Supplemental heatmap shows RDH12-vs-WT matched-dose context
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(tidyr)
})

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------
recovery_file <- "Lipidomics/output_txts/Recovery_filteredlipidtable.csv"
out_dir <- "Lipidomics/Figures/Recovery/DBE_composition/"
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

key_lipid_classes <- c("TG", "DG", "CE", "Cer", "PE", "PC", "PS")
dbe_heatmap_lipid_classes <- c("TG", "DG", "CE", "PE", "PC", "PS")
all_lipid_classes <- c("CAR", "CE", "Cer", "DG", "LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "SM", "TG")

out_long_csv <- file.path(out_dir, "Recovery_DBE_composition_distinct_precursors_long.csv")
out_main_stacked_pdf <- file.path(out_dir, "Recovery_DBE_composition_key_classes_main_text_stacked.pdf")
out_main_stacked_png <- file.path(out_dir, "Recovery_DBE_composition_key_classes_main_text_stacked.png")
out_supp_stacked_pdf <- file.path(out_dir, "Recovery_DBE_composition_all_classes_supplement_stacked.pdf")
out_supp_stacked_png <- file.path(out_dir, "Recovery_DBE_composition_all_classes_supplement_stacked.png")
out_within_delta_csv <- file.path(out_dir, "Recovery_DBE_composition_delta_vs_EtOH_within_genotype.csv")
out_within_delta_pdf <- file.path(out_dir, "Recovery_DBE_composition_delta_vs_EtOH_within_genotype_heatmap.pdf")
out_within_delta_png <- file.path(out_dir, "Recovery_DBE_composition_delta_vs_EtOH_within_genotype_heatmap.png")
out_genotype_delta_csv <- file.path(out_dir, "Recovery_DBE_composition_delta_RDH12_vs_WT_matched_dose.csv")
out_genotype_delta_pdf <- file.path(out_dir, "Recovery_DBE_composition_delta_RDH12_vs_WT_matched_dose_heatmap_supplement.pdf")
out_genotype_delta_png <- file.path(out_dir, "Recovery_DBE_composition_delta_RDH12_vs_WT_matched_dose_heatmap_supplement.png")
out_fisher_csv <- file.path(out_dir, "Recovery_DBE_composition_fisher_tests.csv")

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------
classify_dbe <- function(DB_total) {
  case_when(
    DB_total == 0 ~ "SFA",
    DB_total == 1 ~ "MUFA",
    DB_total >= 2 ~ "PUFA",
    TRUE ~ NA_character_
  )
}

parse_group <- function(group) {
  data.frame(
    Group = group,
    Genotype = ifelse(grepl("^Control", group), "WT", "RDH12"),
    Treatment = case_when(
      grepl("_Veh$", group) ~ "Veh",
      grepl("_100$", group) ~ "100 uM",
      grepl("_200$", group) ~ "200 uM",
      TRUE ~ group
    ),
    stringsAsFactors = FALSE
  )
}

prepare_detected_precursors <- function(df, group, sample_cols, blank_cutoff = 1.3) {
  missing_cols <- setdiff(c(sample_cols, "blank", "precursor", "lipid_class", "DB_total"), colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns for ", group, ":\n", paste(missing_cols, collapse = "\n"))
  }
  
  group_info <- parse_group(group)
  
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
      !is.na(lipid_class),
      !is.na(DB_total)
    ) |>
    transmute(
      Group = group,
      Genotype = group_info$Genotype,
      Treatment = group_info$Treatment,
      lipid_class = as.character(lipid_class),
      unsat_class = classify_dbe(DB_total),
      precursor = precursor
    ) |>
    filter(!is.na(unsat_class)) |>
    distinct(Group, Genotype, Treatment, lipid_class, unsat_class, precursor)
}

complete_dbe_summary <- function(detected_df, lipid_classes) {
  group_info <- bind_rows(lapply(names(group_sample_cols), parse_group))
  
  detected_df |>
    filter(lipid_class %in% lipid_classes) |>
    count(Group, Genotype, Treatment, lipid_class, unsat_class, name = "distinct_precursor") |>
    complete(
      Group = names(group_sample_cols),
      lipid_class = lipid_classes,
      unsat_class = c("SFA", "MUFA", "PUFA"),
      fill = list(distinct_precursor = 0)
    ) |>
    select(-Genotype, -Treatment) |>
    left_join(group_info, by = "Group") |>
    group_by(Group, lipid_class) |>
    mutate(
      class_total = sum(distinct_precursor),
      prop = if_else(class_total > 0, distinct_precursor / class_total, 0)
    ) |>
    ungroup() |>
    mutate(
      Group = factor(Group, levels = names(group_sample_cols)),
      Genotype = factor(Genotype, levels = c("WT", "RDH12")),
      Treatment = factor(Treatment, levels = c("Veh", "100 uM", "200 uM")),
      unsat_class = factor(unsat_class, levels = c("SFA", "MUFA", "PUFA")),
      lipid_class = factor(lipid_class, levels = lipid_classes)
    )
}

make_dbe_stacked_plot <- function(plot_df, title, ncol = 4) {
  total_labels <- plot_df |>
    distinct(Genotype, Treatment, lipid_class, class_total) |>
    mutate(label = if_else(class_total > 0, paste0("n=", class_total), ""))
  
  ggplot(plot_df, aes(x = prop, y = Treatment, fill = unsat_class)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(
      data = total_labels,
      aes(x = 1.03, y = Treatment, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.2,
      color = "grey20"
    ) +
    facet_grid(Genotype ~ lipid_class) +
    scale_x_continuous(limits = c(0, 1.22), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_fill_manual(
      values = c("SFA" = "#4E79A7", "MUFA" = "#59A14F", "PUFA" = "#E15759"),
      labels = c("SFA (0 DBE)", "MUFA (1 DBE)", "PUFA (>=2 DBE)")
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = 8),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      panel.spacing = unit(0.7, "lines"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "grey25"),
      plot.margin = margin(6, 28, 6, 6)
    ) +
    labs(
      title = title,
      subtitle = paste0("Distinct detected precursors; detection threshold = ", blank_cutoff, "x blank"),
      x = "Proportion within lipid class",
      y = NULL
    )
}

get_fisher_stars <- function(padj) {
  case_when(
    is.na(padj) ~ "",
    padj <= 0.01 ~ "***",
    padj <= 0.05 ~ "**",
    padj <= 0.1 ~ "*",
    TRUE ~ ""
  )
}

run_fisher_test <- function(summary_df, group_a, group_b, lipid_class_i) {
  counts_a <- summary_df |>
    filter(Group == group_a, lipid_class == lipid_class_i) |>
    arrange(unsat_class) |>
    pull(distinct_precursor)
  
  counts_b <- summary_df |>
    filter(Group == group_b, lipid_class == lipid_class_i) |>
    arrange(unsat_class) |>
    pull(distinct_precursor)
  
  count_matrix <- rbind(Reference = counts_a, Comparison = counts_b)
  colnames(count_matrix) <- c("SFA", "MUFA", "PUFA")
  
  if (sum(count_matrix[1, ]) > 0 && sum(count_matrix[2, ]) > 0) {
    fisher.test(count_matrix)$p.value
  } else {
    NA_real_
  }
}

build_delta_table <- function(summary_df, comparison_tbl, lipid_classes) {
  delta_rows <- list()
  
  for (i in seq_len(nrow(comparison_tbl))) {
    ref_group <- comparison_tbl$ReferenceGroup[i]
    comp_group <- comparison_tbl$ComparisonGroup[i]
    comparison_label <- comparison_tbl$Comparison[i]
    comparison_type <- comparison_tbl$ComparisonType[i]
    
    ref_df <- summary_df |>
      filter(Group == ref_group, lipid_class %in% lipid_classes) |>
      select(lipid_class, unsat_class, reference_prop = prop, reference_total = class_total)
    
    comp_df <- summary_df |>
      filter(Group == comp_group, lipid_class %in% lipid_classes) |>
      select(lipid_class, unsat_class, comparison_prop = prop, comparison_total = class_total)
    
    delta_rows[[i]] <- ref_df |>
      left_join(comp_df, by = c("lipid_class", "unsat_class")) |>
      mutate(
        Comparison = comparison_label,
        ComparisonType = comparison_type,
        ReferenceGroup = ref_group,
        ComparisonGroup = comp_group,
        delta_prop = comparison_prop - reference_prop,
        delta_label = if_else(abs(delta_prop) < 0.005, "0", sprintf("%+.2f", delta_prop))
      )
  }
  
  bind_rows(delta_rows)
}

run_fisher_tests <- function(summary_df, comparison_tbl, lipid_classes) {
  test_grid <- merge(
    comparison_tbl,
    data.frame(lipid_class = lipid_classes, stringsAsFactors = FALSE),
    by = NULL
  )
  
  test_grid |>
    rowwise() |>
    mutate(p_value = run_fisher_test(summary_df, ReferenceGroup, ComparisonGroup, lipid_class)) |>
    ungroup() |>
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significance = get_fisher_stars(p_adj_BH)
    )
}

make_delta_heatmap <- function(delta_df, title, subtitle) {
  star_df <- delta_df |>
    distinct(lipid_class, Comparison, significance) |>
    filter(significance != "") |>
    mutate(unsat_class = factor("PUFA", levels = c("SFA", "MUFA", "PUFA")))
  
  ggplot(delta_df, aes(x = unsat_class, y = lipid_class, fill = delta_prop)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = delta_label), size = 2.5, color = "grey10") +
    geom_text(
      data = star_df,
      aes(x = unsat_class, y = lipid_class, label = significance),
      inherit.aes = FALSE,
      nudge_x = 0.45,
      size = 3.8,
      fontface = "bold",
      color = "grey10"
    ) +
    scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      limits = c(-0.35, 0.35),
      breaks = c(-0.3, -0.15, 0, 0.15, 0.3),
      labels = c("-0.30", "-0.15", "0", "+0.15", "+0.30"),
      oob = scales::squish
    ) +
    facet_wrap(~ Comparison, nrow = 1) +
    scale_x_discrete(position = "top") +
    coord_fixed(ratio = 0.75, clip = "off", xlim = c(0.5, 3.55)) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      axis.title = element_blank(),
      axis.text.x = element_text(color = "grey15"),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.6, "lines"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "grey25"),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    labs(title = title, subtitle = subtitle, fill = "Delta\nproportion")
}

# ---------------------------------------------------------
# 3. BUILD DBE SUMMARY TABLE
# ---------------------------------------------------------
lipid_df <- read.csv(recovery_file, check.names = TRUE, stringsAsFactors = FALSE)

detected_precursors <- bind_rows(
  lapply(names(group_sample_cols), function(group) {
    prepare_detected_precursors(lipid_df, group, group_sample_cols[[group]], blank_cutoff)
  })
)

dbe_summary <- complete_dbe_summary(detected_precursors, all_lipid_classes)
write.csv(dbe_summary, out_long_csv, row.names = FALSE)

within_comparisons <- tibble::tribble(
  ~Comparison, ~ComparisonType, ~ReferenceGroup, ~ComparisonGroup,
  "WT 100 vs WT EtOH", "Within genotype", "Control_Veh", "Control_100",
  "WT 200 vs WT EtOH", "Within genotype", "Control_Veh", "Control_200",
  "RDH12 100 vs RDH12 EtOH", "Within genotype", "RDH12_Veh", "RDH12_100",
  "RDH12 200 vs RDH12 EtOH", "Within genotype", "RDH12_Veh", "RDH12_200"
)

genotype_comparisons <- tibble::tribble(
  ~Comparison, ~ComparisonType, ~ReferenceGroup, ~ComparisonGroup,
  "RDH12 EtOH vs WT EtOH", "RDH12 vs WT matched dose", "Control_Veh", "RDH12_Veh",
  "RDH12 100 vs WT 100", "RDH12 vs WT matched dose", "Control_100", "RDH12_100",
  "RDH12 200 vs WT 200", "RDH12 vs WT matched dose", "Control_200", "RDH12_200"
)

all_comparisons <- bind_rows(within_comparisons, genotype_comparisons)

fisher_results <- run_fisher_tests(dbe_summary, all_comparisons, dbe_heatmap_lipid_classes)
write.csv(fisher_results, out_fisher_csv, row.names = FALSE)

within_delta <- build_delta_table(dbe_summary, within_comparisons, dbe_heatmap_lipid_classes) |>
  left_join(
    fisher_results |> select(Comparison, lipid_class, p_value, p_adj_BH, significance),
    by = c("Comparison", "lipid_class")
  ) |>
  mutate(
    lipid_class = factor(lipid_class, levels = rev(dbe_heatmap_lipid_classes)),
    unsat_class = factor(unsat_class, levels = c("SFA", "MUFA", "PUFA")),
    Comparison = factor(Comparison, levels = within_comparisons$Comparison)
  )

genotype_delta <- build_delta_table(dbe_summary, genotype_comparisons, dbe_heatmap_lipid_classes) |>
  left_join(
    fisher_results |> select(Comparison, lipid_class, p_value, p_adj_BH, significance),
    by = c("Comparison", "lipid_class")
  ) |>
  mutate(
    lipid_class = factor(lipid_class, levels = rev(dbe_heatmap_lipid_classes)),
    unsat_class = factor(unsat_class, levels = c("SFA", "MUFA", "PUFA")),
    Comparison = factor(Comparison, levels = genotype_comparisons$Comparison)
  )

write.csv(within_delta, out_within_delta_csv, row.names = FALSE)
write.csv(genotype_delta, out_genotype_delta_csv, row.names = FALSE)

# ---------------------------------------------------------
# 4. SAVE FIGURES
# ---------------------------------------------------------
main_df <- dbe_summary |>
  filter(lipid_class %in% key_lipid_classes) |>
  mutate(lipid_class = fct_relevel(lipid_class, key_lipid_classes))

supp_df <- dbe_summary |>
  mutate(lipid_class = fct_relevel(lipid_class, all_lipid_classes))

p_main_stacked <- make_dbe_stacked_plot(
  main_df,
  title = "Recovery DBE composition of key lipid classes"
)

p_supp_stacked <- make_dbe_stacked_plot(
  supp_df,
  title = "Recovery DBE composition of detected lipid classes"
)

p_within_delta <- make_delta_heatmap(
  within_delta,
  title = "Recovery DBE composition shifts within genotype",
  subtitle = "Values: delta proportion vs genotype-matched EtOH; stars: Fisher BH-FDR."
)

p_genotype_delta <- make_delta_heatmap(
  genotype_delta,
  title = "Recovery DBE composition differences by genotype",
  subtitle = "Values: RDH12 minus WT at matched treatment; stars: Fisher BH-FDR."
)

ggsave(out_main_stacked_pdf, p_main_stacked, width = 8.2, height = 4.2, useDingbats = FALSE)
ggsave(out_main_stacked_png, p_main_stacked, width = 8.2, height = 4.2, dpi = 600)
ggsave(out_supp_stacked_pdf, p_supp_stacked, width = 11.5, height = 4.8, useDingbats = FALSE)
ggsave(out_supp_stacked_png, p_supp_stacked, width = 11.5, height = 4.8, dpi = 600)
ggsave(out_within_delta_pdf, p_within_delta, width = 7.2, height = 3.6, useDingbats = FALSE)
ggsave(out_within_delta_png, p_within_delta, width = 7.2, height = 3.6, dpi = 600)
ggsave(out_genotype_delta_pdf, p_genotype_delta, width = 5.8, height = 3.6, useDingbats = FALSE)
ggsave(out_genotype_delta_png, p_genotype_delta, width = 5.8, height = 3.6, dpi = 600)

# ---------------------------------------------------------
# 5. VERIFICATION SUMMARY
# ---------------------------------------------------------
prop_check <- dbe_summary |>
  group_by(Group, lipid_class) |>
  summarise(prop_sum = sum(prop), class_total = first(class_total), .groups = "drop") |>
  filter(class_total > 0, abs(prop_sum - 1) > 1e-8)

if (nrow(prop_check) > 0) {
  stop("Some Group + lipid_class proportions do not sum to 1.")
}

print(
  dbe_summary |>
    distinct(Group, lipid_class, class_total) |>
    tidyr::pivot_wider(names_from = Group, values_from = class_total) |>
    arrange(lipid_class)
)

cat("\nSaved files:\n")
cat(out_long_csv, "\n")
cat(out_main_stacked_pdf, "\n")
cat(out_main_stacked_png, "\n")
cat(out_supp_stacked_pdf, "\n")
cat(out_supp_stacked_png, "\n")
cat(out_within_delta_csv, "\n")
cat(out_within_delta_pdf, "\n")
cat(out_within_delta_png, "\n")
cat(out_genotype_delta_csv, "\n")
cat(out_genotype_delta_pdf, "\n")
cat(out_genotype_delta_png, "\n")
cat(out_fisher_csv, "\n")
