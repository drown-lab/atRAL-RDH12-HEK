# Supplemental acute RDH12 GO Biological Process and Reactome enrichment figure
#
# This script rebuilds six enrichment dotplot panels from the curated
# manual-filtered batch enrichment tables and exports a grouped supplemental
# figure as PDF and SVG.

# ---------------------------------------------------------
# 1. LOAD PACKAGES
# ---------------------------------------------------------

required_packages <- c("dplyr", "forcats", "ggplot2", "patchwork", "readr", "stringr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(dplyr)
library(forcats)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)

# ---------------------------------------------------------
# 2. USER-EDITABLE SETTINGS
# ---------------------------------------------------------

output_dir <- "Proteomic_Figs/Supplemental_Acute_RDH12_GO_Reactome"
output_base <- "acute_RDH12_5hr_GO_Reactome_supplemental"

figure_width <- 11
figure_height <- 14

max_terms_per_panel <- 16

panel_specs <- tibble::tribble(
  ~panel, ~panel_title, ~input_file,
  "A", "GO BP: Veh enriched, 100 uM atRAL vs Veh",
  "GOanalysis/output/clusterProfiler_batch/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_down/enrichGO_BP_manual_filtered.csv",
  "B", "GO BP: Veh enriched, 200 uM atRAL vs Veh",
  "GOanalysis/output/clusterProfiler_batch/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_down/enrichGO_BP_manual_filtered.csv",
  "C", "Reactome: 100 uM atRAL enriched",
  "PathwayAnalysis/output/Reactome_batch/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_up/Reactome_manual_filtered.csv",
  "D", "Reactome: 200 uM atRAL enriched",
  "PathwayAnalysis/output/Reactome_batch/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_up/Reactome_manual_filtered.csv",
  "E", "Reactome: Veh enriched, 100 uM atRAL vs Veh",
  "PathwayAnalysis/output/Reactome_batch/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_down/Reactome_manual_filtered.csv",
  "F", "Reactome: Veh enriched, 200 uM atRAL vs Veh",
  "PathwayAnalysis/output/Reactome_batch/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_down/Reactome_manual_filtered.csv"
)

# ---------------------------------------------------------
# 3. HELPER FUNCTIONS
# ---------------------------------------------------------

check_input_files <- function(specs) {
  missing_files <- specs$input_file[!file.exists(specs$input_file)]

  if (length(missing_files) > 0) {
    stop(
      "Missing input files:\n",
      paste(missing_files, collapse = "\n")
    )
  }
}

pick_column <- function(df, candidates, label, file) {
  found <- candidates[candidates %in% colnames(df)][1]

  if (is.na(found) || length(found) == 0) {
    stop(
      "Could not find ", label, " column in: ", file,
      "\nExpected one of: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(df), collapse = ", ")
    )
  }

  found
}

read_panel_table <- function(panel, panel_title, input_file, max_terms) {
  enrich_df <- readr::read_csv(input_file, show_col_types = FALSE)

  if (nrow(enrich_df) == 0) {
    stop("Input table has no rows: ", input_file)
  }

  term_col <- pick_column(
    enrich_df,
    candidates = c("Description", "Gene Set Name", "Term", "Pathway"),
    label = "term/description",
    file = input_file
  )
  fold_col <- pick_column(
    enrich_df,
    candidates = c("FoldEnrichment_num", "FoldEnrichment", "Fold Enrichment"),
    label = "fold enrichment",
    file = input_file
  )
  fdr_col <- pick_column(
    enrich_df,
    candidates = c("p.adjust", "FDR", "Benjamini & Hochberg", "qvalue"),
    label = "FDR/BH adjusted p-value",
    file = input_file
  )
  count_col <- pick_column(
    enrich_df,
    candidates = c("Count", "Count Overlap Gene", "overlap_count"),
    label = "overlap count",
    file = input_file
  )

  plot_df <- enrich_df %>%
    transmute(
      panel = panel,
      panel_title = panel_title,
      term = .data[[term_col]],
      fold_enrichment = suppressWarnings(as.numeric(.data[[fold_col]])),
      fdr = suppressWarnings(as.numeric(.data[[fdr_col]])),
      count = suppressWarnings(as.numeric(.data[[count_col]]))
    ) %>%
    filter(
      !is.na(term),
      !is.na(fold_enrichment),
      !is.na(fdr),
      !is.na(count)
    ) %>%
    arrange(fdr, desc(fold_enrichment)) %>%
    slice_head(n = max_terms) %>%
    mutate(
      term_wrapped = stringr::str_wrap(term, width = 36),
      term_wrapped = forcats::fct_reorder(term_wrapped, fold_enrichment)
    )

  if (nrow(plot_df) == 0) {
    stop("No plottable rows remained after numeric column checks: ", input_file)
  }

  plot_df
}

make_panel_plot <- function(plot_df) {
  fdr_limits <- range(plot_df$fdr, na.rm = TRUE)
  count_limits <- range(plot_df$count, na.rm = TRUE)

  ggplot(
    plot_df,
    aes(
      x = fold_enrichment,
      y = term_wrapped,
      size = count,
      color = fdr
    )
  ) +
    geom_point(alpha = 0.95) +
    scale_color_viridis_c(
      option = "plasma",
      direction = -1,
      name = "p.adjust",
      limits = fdr_limits
    ) +
    scale_size_continuous(
      name = "Count",
      range = c(1.4, 5.2),
      limits = count_limits
    ) +
    labs(
      title = NULL,
      x = "Fold enrichment",
      y = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.y = element_text(size = 6.8, color = "black", lineheight = 0.9),
      axis.ticks.y = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length.y = grid::unit(1.5, "mm"),
      axis.text.x = element_text(size = 7, color = "black"),
      axis.title.x = element_text(size = 8),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.key.width = grid::unit(0.36, "cm"),
      plot.margin = margin(8, 12, 8, 26)
    )
}

# ---------------------------------------------------------
# 4. IMPORT AND STANDARDIZE INPUTS
# ---------------------------------------------------------

check_input_files(panel_specs)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_data <- panel_specs %>%
  rowwise() %>%
  mutate(
    data = list(read_panel_table(panel, panel_title, input_file, max_terms_per_panel))
  ) %>%
  ungroup()

panel_plots <- panel_data %>%
  mutate(
    plot = lapply(data, make_panel_plot)
  )

named_plots <- stats::setNames(panel_plots$plot, panel_plots$panel)

# ---------------------------------------------------------
# 5. BUILD AND SAVE GROUPED FIGURE
# ---------------------------------------------------------

supplemental_figure <- (
  named_plots$A + named_plots$B
) / (
  named_plots$C + named_plots$D
) / (
  named_plots$E + named_plots$F
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 13),
    plot.tag.position = c(0, 1)
  )

pdf_file <- file.path(output_dir, paste0(output_base, ".pdf"))
svg_file <- file.path(output_dir, paste0(output_base, ".svg"))

ggsave(
  filename = pdf_file,
  plot = supplemental_figure,
  width = figure_width,
  height = figure_height,
  units = "in",
  device = grDevices::cairo_pdf
)

ggsave(
  filename = svg_file,
  plot = supplemental_figure,
  width = figure_width,
  height = figure_height,
  units = "in",
  device = grDevices::svg
)

# ---------------------------------------------------------
# 6. COMPLETION SUMMARY
# ---------------------------------------------------------

summary_tbl <- panel_data %>%
  transmute(
    panel,
    panel_title,
    input_file,
    n_terms_plotted = vapply(data, nrow, integer(1))
  )

print(summary_tbl)
message("Saved supplemental figure PDF: ", pdf_file)
message("Saved supplemental figure SVG: ", svg_file)
