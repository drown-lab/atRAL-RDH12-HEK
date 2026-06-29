# =========================================================
# Protein-level treatment-response ratio scatterplots
# - Genotype comparison plots: RDH12 treatment vs vehicle on x,
#   GFP treatment vs vehicle on y
# - Within-genotype plots: 100 uM treatment vs vehicle on x,
#   200 uM treatment vs vehicle on y
# =========================================================

# -------------------------
# 1. Load packages
# -------------------------
packages_needed <- c("dplyr", "ggplot2", "ggrepel", "readr", "stringr", "svglite", "tibble")

missing_packages <- packages_needed[!vapply(packages_needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(stringr)

# =========================================================
# 2. User-editable parameters
# =========================================================

input_file <- "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv"
out_dir <- "Proteomic_Figs/Treatment_response_ratio_scatters"

padj_cutoff <- 0.05
lfc_cutoff <- log2(1.5)
n_top_labels <- 15

out_table <- file.path(out_dir, "treatment_response_ratio_scatter_table.csv")

# =========================================================
# 3. Define contrasts and plots
# =========================================================

contrast_cols <- tibble::tribble(
  ~contrast_id, ~contrast_name,
  "rdh12_100_recovery", "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
  "rdh12_200_recovery", "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
  "gfp_100_recovery", "GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr",
  "gfp_200_recovery", "GFP_200_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr"
) %>%
  mutate(
    ratio_col = paste0(contrast_name, "_ratio"),
    padj_col = paste0(contrast_name, "_p.adj")
  )

plot_specs <- tibble::tribble(
  ~plot_id, ~plot_title, ~x_contrast_id, ~y_contrast_id, ~x_label, ~y_label,
  "genotype_comparison_100uM",
  "100 uM treatment response: RDH12 vs GFP",
  "rdh12_100_recovery",
  "gfp_100_recovery",
  "RDH12 100 uM vs vehicle log2 ratio",
  "GFP 100 uM vs vehicle log2 ratio",

  "genotype_comparison_200uM",
  "200 uM treatment response: RDH12 vs GFP",
  "rdh12_200_recovery",
  "gfp_200_recovery",
  "RDH12 200 uM vs vehicle log2 ratio",
  "GFP 200 uM vs vehicle log2 ratio",

  "rdh12_within_genotype_100uM_vs_200uM",
  "RDH12 treatment response: 100 uM vs 200 uM",
  "rdh12_100_recovery",
  "rdh12_200_recovery",
  "RDH12 100 uM vs vehicle log2 ratio",
  "RDH12 200 uM vs vehicle log2 ratio",

  "gfp_within_genotype_100uM_vs_200uM",
  "GFP treatment response: 100 uM vs 200 uM",
  "gfp_100_recovery",
  "gfp_200_recovery",
  "GFP 100 uM vs vehicle log2 ratio",
  "GFP 200 uM vs vehicle log2 ratio"
)

# =========================================================
# 4. Import data and check required columns
# =========================================================

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dep_results <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c("name", "ID", contrast_cols$ratio_col, contrast_cols$padj_col)
missing_cols <- setdiff(required_cols, colnames(dep_results))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# =========================================================
# 5. Prepare plotting data
# =========================================================

call_axis_significance <- function(ratio, padj, padj_cutoff, lfc_cutoff) {
  !is.na(padj) &
    !is.na(ratio) &
    padj <= padj_cutoff &
    abs(ratio) >= lfc_cutoff
}

make_plot_df <- function(plot_spec, dep_results, contrast_cols, padj_cutoff, lfc_cutoff) {
  x_info <- contrast_cols %>%
    filter(.data$contrast_id == plot_spec$x_contrast_id)
  y_info <- contrast_cols %>%
    filter(.data$contrast_id == plot_spec$y_contrast_id)

  if (nrow(x_info) != 1 || nrow(y_info) != 1) {
    stop("Plot spec has an unknown contrast id: ", plot_spec$plot_id)
  }

  plot_df <- dep_results %>%
    transmute(
      plot_id = plot_spec$plot_id,
      plot_title = plot_spec$plot_title,
      x_label = plot_spec$x_label,
      y_label = plot_spec$y_label,
      name = as.character(.data$name),
      ID = as.character(.data$ID),
      x_ratio = suppressWarnings(as.numeric(.data[[x_info$ratio_col]])),
      y_ratio = suppressWarnings(as.numeric(.data[[y_info$ratio_col]])),
      x_padj = suppressWarnings(as.numeric(.data[[x_info$padj_col]])),
      y_padj = suppressWarnings(as.numeric(.data[[y_info$padj_col]]))
    ) %>%
    filter(!is.na(.data$x_ratio), !is.na(.data$y_ratio)) %>%
    mutate(
      significant_x = call_axis_significance(
        ratio = .data$x_ratio,
        padj = .data$x_padj,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff
      ),
      significant_y = call_axis_significance(
        ratio = .data$y_ratio,
        padj = .data$y_padj,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff
      ),
      significance_class = case_when(
        .data$significant_x & .data$significant_y ~ "Significant on both axes",
        .data$significant_x & !.data$significant_y ~ "Significant on x-axis only",
        !.data$significant_x & .data$significant_y ~ "Significant on y-axis only",
        TRUE ~ "Not significant"
      ),
      significant_either_axis = .data$significant_x | .data$significant_y,
      combined_abs_ratio = abs(.data$x_ratio) + abs(.data$y_ratio)
    )

  if (nrow(plot_df) == 0) {
    stop("No rows remain after removing missing x/y ratios for plot: ", plot_spec$plot_id)
  }

  plot_df
}

plot_tables <- lapply(
  seq_len(nrow(plot_specs)),
  function(i) {
    make_plot_df(
      plot_spec = plot_specs[i, ],
      dep_results = dep_results,
      contrast_cols = contrast_cols,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    )
  }
)

plot_df_all <- bind_rows(plot_tables) %>%
  mutate(
    significance_class = factor(
      .data$significance_class,
      levels = c(
        "Significant on both axes",
        "Significant on x-axis only",
        "Significant on y-axis only",
        "Not significant"
      )
    )
  )

readr::write_csv(plot_df_all, out_table)

# =========================================================
# 6. Plot helper and export figures
# =========================================================

significance_colors <- c(
  "Significant on both axes" = "#7B3294",
  "Significant on x-axis only" = "#D73027",
  "Significant on y-axis only" = "#4575B4",
  "Not significant" = "grey78"
)

make_scatter <- function(plot_df, padj_cutoff, lfc_cutoff, n_top_labels) {
  label_df <- plot_df %>%
    filter(.data$significant_either_axis) %>%
    arrange(desc(.data$combined_abs_ratio)) %>%
    slice_head(n = n_top_labels)

  axis_limit <- max(abs(c(plot_df$x_ratio, plot_df$y_ratio, lfc_cutoff)), na.rm = TRUE)
  axis_limit <- ceiling(axis_limit * 10) / 10

  ggplot(plot_df, aes(x = .data$x_ratio, y = .data$y_ratio)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
    geom_hline(
      yintercept = c(-lfc_cutoff, lfc_cutoff),
      color = "grey35",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      color = "grey35",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dotted", linewidth = 0.5) +
    geom_point(
      data = plot_df %>% filter(.data$significance_class == "Not significant"),
      aes(color = .data$significance_class),
      alpha = 0.55,
      size = 2.1
    ) +
    geom_point(
      data = plot_df %>% filter(.data$significance_class != "Not significant"),
      aes(color = .data$significance_class),
      alpha = 0.9,
      size = 2.6
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = .data$name),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.35,
      point.padding = 0.25,
      seed = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = significance_colors, drop = FALSE) +
    coord_equal(xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit)) +
    labs(
      title = unique(plot_df$plot_title),
      subtitle = paste0(
        "Significant if padj <= ",
        padj_cutoff,
        " and abs(log2 ratio) >= ",
        round(lfc_cutoff, 3)
      ),
      x = unique(plot_df$x_label),
      y = unique(plot_df$y_label),
      color = "Class"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

exported_files <- c()

for (plot_id in plot_specs$plot_id) {
  plot_df <- plot_df_all %>%
    filter(.data$plot_id == .env$plot_id)

  scatter_plot <- make_scatter(
    plot_df = plot_df,
    padj_cutoff = padj_cutoff,
    lfc_cutoff = lfc_cutoff,
    n_top_labels = n_top_labels
  )

  out_svg <- file.path(out_dir, paste0(plot_id, ".svg"))
  out_pdf <- file.path(out_dir, paste0(plot_id, ".pdf"))

  ggsave(out_svg, plot = scatter_plot, width = 8.5, height = 7.2, units = "in")
  ggsave(out_pdf, plot = scatter_plot, width = 8.5, height = 7.2, units = "in")

  exported_files <- c(exported_files, out_svg, out_pdf)
}

# =========================================================
# 7. Completion summary
# =========================================================

summary_counts <- plot_df_all %>%
  count(.data$plot_id, .data$significance_class, name = "n") %>%
  arrange(.data$plot_id, .data$significance_class)

print(summary_counts)

message("\nTreatment-response ratio scatterplots completed.")
message("Rows in combined plotting table: ", nrow(plot_df_all))
message("Table written to: ", out_table)
message("Figure files written:")
for (file in exported_files) {
  message("  ", file)
}
