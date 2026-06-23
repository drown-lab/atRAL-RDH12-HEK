# =========================================================
# Acute lipidomics DBE composition figures
# - Counts distinct detected precursors by lipid class and DBE group
# - Generates a key-class main-text panel and an all-class supplement panel
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
acute_file <- "Lipidomics/output_txts/Acute_filteredlipidtable.csv"
out_dir <- "Lipidomics/Figures/Acute/DBE_composition/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

blank_cutoff <- 1.3

treatment_sample_cols <- list(
  "Veh" = c("s1", "s2", "s3"),
  "100 uM" = c("s4", "s5", "s6"),
  "200 uM" = c("s7", "s8", "s9")
)

key_lipid_classes <- c("TG", "DG", "CE", "Cer", "PE", "PC", "PS")
all_lipid_classes <- c("CAR", "CE", "Cer", "DG", "LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "SM", "TG")

out_long_csv <- file.path(out_dir, "Acute_DBE_composition_distinct_precursors_long.csv")
out_main_pdf <- file.path(out_dir, "Acute_DBE_composition_key_classes_main_text.pdf")
out_main_png <- file.path(out_dir, "Acute_DBE_composition_key_classes_main_text.png")
out_supp_pdf <- file.path(out_dir, "Acute_DBE_composition_all_classes_supplement.pdf")
out_supp_png <- file.path(out_dir, "Acute_DBE_composition_all_classes_supplement.png")
out_delta_csv <- file.path(out_dir, "Acute_DBE_composition_delta_vs_vehicle.csv")
out_delta_pdf <- file.path(out_dir, "Acute_DBE_composition_delta_vs_vehicle_heatmap.pdf")
out_delta_png <- file.path(out_dir, "Acute_DBE_composition_delta_vs_vehicle_heatmap.png")
out_fisher_csv <- file.path(out_dir, "Acute_DBE_composition_fisher_tests_vs_vehicle.csv")

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------
# DBE groups are composition bins for detected precursors, not abundance calls.
classify_dbe <- function(DB_total) {
  case_when(
    DB_total == 0 ~ "SFA",
    DB_total == 1 ~ "MUFA",
    DB_total >= 2 ~ "PUFA",
    TRUE ~ NA_character_
  )
}

prepare_detected_precursors <- function(df, treatment, sample_cols, blank_cutoff = 1.3) {
  missing_cols <- setdiff(c(sample_cols, "blank", "precursor", "lipid_class", "DB_total"), colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "These required columns are missing for ", treatment, ":\n",
      paste(missing_cols, collapse = "\n")
    )
  }
  
  df %>%
    rowwise() %>%
    mutate(
      max_treatment_signal = max(c_across(all_of(sample_cols)), na.rm = TRUE),
      signal_to_blank = max_treatment_signal / blank
    ) %>%
    ungroup() %>%
    filter(
      is.finite(signal_to_blank),
      signal_to_blank >= blank_cutoff,
      !is.na(precursor),
      !is.na(lipid_class),
      !is.na(DB_total)
    ) %>%
    transmute(
      Treatment = treatment,
      lipid_class = as.character(lipid_class),
      unsat_class = classify_dbe(DB_total),
      precursor = precursor
    ) %>%
    filter(!is.na(unsat_class)) %>%
    distinct(Treatment, lipid_class, unsat_class, precursor)
}

complete_dbe_summary <- function(detected_df, lipid_classes) {
  detected_df %>%
    filter(lipid_class %in% lipid_classes) %>%
    count(Treatment, lipid_class, unsat_class, name = "distinct_precursor") %>%
    complete(
      Treatment = names(treatment_sample_cols),
      lipid_class = lipid_classes,
      unsat_class = c("SFA", "MUFA", "PUFA"),
      fill = list(distinct_precursor = 0)
    ) %>%
    group_by(Treatment, lipid_class) %>%
    mutate(
      class_total = sum(distinct_precursor),
      prop = if_else(class_total > 0, distinct_precursor / class_total, 0)
    ) %>%
    ungroup() %>%
    mutate(
      Treatment = factor(Treatment, levels = names(treatment_sample_cols)),
      unsat_class = factor(unsat_class, levels = c("SFA", "MUFA", "PUFA")),
      lipid_class = factor(lipid_class, levels = lipid_classes)
    )
}

make_dbe_plot <- function(plot_df, title, ncol = 4) {
  total_labels <- plot_df %>%
    distinct(Treatment, lipid_class, class_total) %>%
    mutate(label = if_else(class_total > 0, paste0("n=", class_total), ""))
  
  ggplot(plot_df, aes(x = prop, y = Treatment, fill = unsat_class)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(
      data = total_labels,
      aes(x = 1.03, y = Treatment, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.5,
      color = "grey20"
    ) +
    facet_wrap(~ lipid_class, ncol = ncol) +
    scale_x_continuous(
      limits = c(0, 1.18),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1")
    ) +
    scale_fill_manual(
      values = c(
        "SFA" = "#4E79A7",
        "MUFA" = "#59A14F",
        "PUFA" = "#E15759"
      ),
      labels = c("SFA (0 DBE)", "MUFA (1 DBE)", "PUFA (>=2 DBE)")
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = 9),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      panel.spacing = unit(0.9, "lines"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "grey25"),
      plot.margin = margin(6, 26, 6, 6)
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

run_dbe_fisher_tests <- function(summary_df) {
  test_grid <- expand.grid(
    lipid_class = levels(summary_df$lipid_class),
    Treatment = setdiff(names(treatment_sample_cols), "Veh"),
    stringsAsFactors = FALSE
  )
  
  test_results <- lapply(seq_len(nrow(test_grid)), function(i) {
    lipid_class_i <- test_grid$lipid_class[i]
    treatment_i <- test_grid$Treatment[i]
    
    vehicle_counts <- summary_df %>%
      filter(lipid_class == lipid_class_i, Treatment == "Veh") %>%
      arrange(unsat_class) %>%
      pull(distinct_precursor)
    
    treatment_counts <- summary_df %>%
      filter(lipid_class == lipid_class_i, Treatment == treatment_i) %>%
      arrange(unsat_class) %>%
      pull(distinct_precursor)
    
    count_matrix <- rbind(
      Veh = vehicle_counts,
      Treatment = treatment_counts
    )
    colnames(count_matrix) <- c("MUFA", "PUFA", "SFA")[order(factor(c("MUFA", "PUFA", "SFA"), levels = c("SFA", "MUFA", "PUFA")))]
    colnames(count_matrix) <- c("SFA", "MUFA", "PUFA")
    
    raw_p <- if (sum(count_matrix[1, ]) > 0 && sum(count_matrix[2, ]) > 0) {
      fisher.test(count_matrix)$p.value
    } else {
      NA_real_
    }
    
    data.frame(
      lipid_class = lipid_class_i,
      Treatment = treatment_i,
      veh_SFA = count_matrix["Veh", "SFA"],
      veh_MUFA = count_matrix["Veh", "MUFA"],
      veh_PUFA = count_matrix["Veh", "PUFA"],
      treatment_SFA = count_matrix["Treatment", "SFA"],
      treatment_MUFA = count_matrix["Treatment", "MUFA"],
      treatment_PUFA = count_matrix["Treatment", "PUFA"],
      p_value = raw_p,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(test_results) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significance = get_fisher_stars(p_adj_BH)
    )
}

make_delta_heatmap <- function(delta_df) {
  star_df <- delta_df %>%
    distinct(lipid_class, Treatment, significance) %>%
    filter(significance != "") %>%
    mutate(unsat_class = factor("PUFA", levels = c("SFA", "MUFA", "PUFA")))
  
  ggplot(delta_df, aes(x = unsat_class, y = lipid_class, fill = delta_prop)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(
      aes(label = delta_label),
      size = 2.7,
      color = "grey10"
    ) +
    geom_text(
      data = star_df,
      aes(x = unsat_class, y = lipid_class, label = significance),
      inherit.aes = FALSE,
      nudge_x = 0.45,
      size = 4,
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
    facet_wrap(~ Treatment, nrow = 1) +
    scale_x_discrete(position = "top") +
    coord_fixed(ratio = 0.75, clip = "off", xlim = c(0.5, 3.55)) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      axis.title = element_blank(),
      axis.text.x = element_text(color = "grey15"),
      axis.text.y = element_text(color = "grey15"),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.8, "lines"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "grey25"),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    labs(
      title = "DBE composition shifts in key lipid classes",
      subtitle = "Values: delta proportion vs vehicle; stars: Fisher BH-FDR",
      fill = "Delta\nproportion"
    )
}

# ---------------------------------------------------------
# 3. BUILD DBE SUMMARY TABLE
# ---------------------------------------------------------
lipid_df <- read.csv(
  acute_file,
  check.names = TRUE,
  stringsAsFactors = FALSE
)

detected_precursors <- bind_rows(
  lapply(names(treatment_sample_cols), function(treatment) {
    prepare_detected_precursors(
      df = lipid_df,
      treatment = treatment,
      sample_cols = treatment_sample_cols[[treatment]],
      blank_cutoff = blank_cutoff
    )
  })
)

dbe_summary <- complete_dbe_summary(detected_precursors, all_lipid_classes)

write.csv(dbe_summary, out_long_csv, row.names = FALSE)

fisher_results <- run_dbe_fisher_tests(dbe_summary)
write.csv(fisher_results, out_fisher_csv, row.names = FALSE)

delta_summary <- dbe_summary %>%
  filter(lipid_class %in% key_lipid_classes) %>%
  select(Treatment, lipid_class, unsat_class, prop, class_total) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(prop, class_total),
    names_sep = "__"
  ) %>%
  pivot_longer(
    cols = c("prop__100 uM", "prop__200 uM"),
    names_to = "Treatment",
    values_to = "treatment_prop"
  ) %>%
  mutate(
    Treatment = sub("^prop__", "", Treatment),
    vehicle_prop = `prop__Veh`,
    delta_prop = treatment_prop - vehicle_prop,
    treatment_total = if_else(Treatment == "100 uM", `class_total__100 uM`, `class_total__200 uM`),
    vehicle_total = `class_total__Veh`,
    dbe_treatment = paste(unsat_class, Treatment, sep = "\n"),
    delta_label = if_else(
      abs(delta_prop) < 0.005,
      "0",
      sprintf("%+.2f", delta_prop)
    ),
    lipid_class = factor(lipid_class, levels = rev(key_lipid_classes)),
    unsat_class = factor(unsat_class, levels = c("SFA", "MUFA", "PUFA")),
    Treatment = factor(Treatment, levels = c("100 uM", "200 uM")),
    dbe_treatment = factor(dbe_treatment)
  ) %>%
  left_join(
    fisher_results %>%
      select(lipid_class, Treatment, p_value, p_adj_BH, significance),
    by = c("lipid_class", "Treatment")
  ) %>%
  arrange(lipid_class, unsat_class, Treatment)

write.csv(delta_summary, out_delta_csv, row.names = FALSE)

# ---------------------------------------------------------
# 4. SAVE MAIN-TEXT AND SUPPLEMENTAL FIGURES
# ---------------------------------------------------------
main_df <- dbe_summary %>%
  filter(lipid_class %in% key_lipid_classes) %>%
  mutate(lipid_class = fct_relevel(lipid_class, key_lipid_classes))

supp_df <- dbe_summary %>%
  mutate(lipid_class = fct_relevel(lipid_class, all_lipid_classes))

p_main <- make_dbe_plot(
  plot_df = main_df,
  title = "DBE composition of key acute lipid classes",
  ncol = 4
)

p_supp <- make_dbe_plot(
  plot_df = supp_df,
  title = "DBE composition of detected acute lipid classes",
  ncol = 5
)

p_delta <- make_delta_heatmap(delta_summary)

ggsave(out_main_pdf, p_main, width = 6.6, height = 4.2, useDingbats = FALSE)
ggsave(out_main_png, p_main, width = 6.6, height = 4.2, dpi = 600)

ggsave(out_supp_pdf, p_supp, width = 8.4, height = 5.2, useDingbats = FALSE)
ggsave(out_supp_png, p_supp, width = 8.4, height = 5.2, dpi = 600)

ggsave(out_delta_pdf, p_delta, width = 5.2, height = 3.6, useDingbats = FALSE)
ggsave(out_delta_png, p_delta, width = 5.2, height = 3.6, dpi = 600)

# ---------------------------------------------------------
# 5. VERIFICATION SUMMARY
# ---------------------------------------------------------
prop_check <- dbe_summary %>%
  group_by(Treatment, lipid_class) %>%
  summarise(prop_sum = sum(prop), class_total = first(class_total), .groups = "drop") %>%
  filter(class_total > 0, abs(prop_sum - 1) > 1e-8)

if (nrow(prop_check) > 0) {
  stop("Some Treatment + lipid_class proportions do not sum to 1.")
}

print(
  dbe_summary %>%
    distinct(Treatment, lipid_class, class_total) %>%
    tidyr::pivot_wider(names_from = Treatment, values_from = class_total) %>%
    arrange(lipid_class)
)

cat("\nSaved files:\n")
cat(out_long_csv, "\n")
cat(out_main_pdf, "\n")
cat(out_main_png, "\n")
cat(out_supp_pdf, "\n")
cat(out_supp_png, "\n")
cat(out_delta_csv, "\n")
cat(out_delta_pdf, "\n")
cat(out_delta_png, "\n")
cat(out_fisher_csv, "\n")
