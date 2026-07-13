# =========================================================
# Acute RDH12 prefoldin heatmap from DEPresults_v2 centered values
# - Shows PFDN1, PFDN2, VBP1/PFDN3, PFDN4, PFDN5, and PFDN6
# - Uses the same centered columns and blue-white-red scale as HeatMap_Fig1.R
# - Keeps all requested prefoldin proteins, regardless of significance
# =========================================================

# ---------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------

packages_needed <- c("dplyr", "tidyr", "pheatmap", "stringr", "tibble")

missing_packages <- packages_needed[!vapply(packages_needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(dplyr)
library(tidyr)
library(pheatmap)
library(stringr)
library(tibble)

# ---------------------------------------------------------
# 2. User inputs
# ---------------------------------------------------------

input_file <- "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv"
out_dir <- "Proteomic_Figs/enriched"

out_file_pdf <- file.path(out_dir, "prefoldin_heatmap_centered_acute.pdf")
out_file_svg <- file.path(out_dir, "prefoldin_heatmap_centered_acute.svg")
out_file_csv <- file.path(out_dir, "prefoldin_heatmap_input_table.csv")

# Acute centered columns to display
acute_centered_cols <- c(
  "RDH12_control_atRAL5hr_centered",
  "RDH12_100_atRAL5hr_centered",
  "RDH12_200_atRAL5hr_centered"
)

acute_col_labels <- c(
  "RDH12_control_atRAL5hr_centered" = "Veh",
  "RDH12_100_atRAL5hr_centered" = "100uM",
  "RDH12_200_atRAL5hr_centered" = "200uM"
)

# Acute contrasts used for row annotations
acute_contrasts <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr"
)

acute_contrast_labels <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr" = "100uM vs Veh",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr" = "200uM vs 100uM",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr" = "200uM vs Veh"
)

# Match HeatMap_Fig1.R thresholds for the row annotation flags.
padj_cutoff <- 0.05
lfc_cutoff <- log2(1.5)

prefoldin_tbl <- tibble::tibble(
  Gene = c("PFDN1", "PFDN2", "VBP1", "PFDN4", "PFDN5", "PFDN6"),
  DisplayGene = c("PFDN1", "PFDN2", "VBP1/PFDN3", "PFDN4", "PFDN5", "PFDN6"),
  order_index = seq_len(6)
)

# ---------------------------------------------------------
# 3. Import data and check required columns
# ---------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

DEPresults_v2 <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

sig_required_cols <- unlist(lapply(
  acute_contrasts,
  function(x) c(paste0(x, "_p.adj"), paste0(x, "_ratio"))
))

required_cols <- c("Gene", "name", "ID", acute_centered_cols, sig_required_cols)
missing_cols <- setdiff(required_cols, colnames(DEPresults_v2))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from DEPresults_v2:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 4. Flag acute contrast significance
# ---------------------------------------------------------

DEPresults_sig <- DEPresults_v2 %>%
  mutate(
    GeneClean = dplyr::coalesce(
      na_if(as.character(.data$Gene), ""),
      na_if(as.character(.data$name), "")
    )
  )

for (ct in acute_contrasts) {
  padj_col <- paste0(ct, "_p.adj")
  ratio_col <- paste0(ct, "_ratio")
  sig_col <- paste0(ct, "_sigflag")

  DEPresults_sig[[sig_col]] <-
    !is.na(DEPresults_sig[[padj_col]]) &
    !is.na(DEPresults_sig[[ratio_col]]) &
    suppressWarnings(as.numeric(DEPresults_sig[[padj_col]])) <= padj_cutoff &
    abs(suppressWarnings(as.numeric(DEPresults_sig[[ratio_col]]))) >= lfc_cutoff
}

sigflag_cols <- paste0(acute_contrasts, "_sigflag")

DEPresults_sig <- DEPresults_sig %>%
  mutate(
    SigIn = apply(
      select(., all_of(sigflag_cols)),
      1,
      function(x) {
        hit_names <- acute_contrast_labels[acute_contrasts[as.logical(x)]]
        if (length(hit_names) == 0) {
          "None"
        } else {
          paste(hit_names, collapse = "; ")
        }
      }
    )
  )

# ---------------------------------------------------------
# 5. Subset to prefoldin proteins
# ---------------------------------------------------------

plot_df <- prefoldin_tbl %>%
  left_join(DEPresults_sig, by = c("Gene" = "GeneClean")) %>%
  arrange(order_index)

missing_prefoldins <- plot_df %>%
  filter(is.na(.data$name) & is.na(.data$ID)) %>%
  pull(.data$DisplayGene)

if (length(missing_prefoldins) > 0) {
  stop(
    "These requested prefoldin proteins were not found in the input table:\n",
    paste(missing_prefoldins, collapse = "\n")
  )
}

plot_df <- plot_df %>%
  mutate(
    across(all_of(c(acute_centered_cols, sig_required_cols)), ~ suppressWarnings(as.numeric(.x)))
  )

# ---------------------------------------------------------
# 6. Build heatmap matrix
# ---------------------------------------------------------

heat_mat <- plot_df %>%
  select(DisplayGene, all_of(acute_centered_cols)) %>%
  as.data.frame()

rownames(heat_mat) <- heat_mat$DisplayGene
heat_mat$DisplayGene <- NULL
heat_mat <- as.matrix(heat_mat)
colnames(heat_mat) <- acute_col_labels[colnames(heat_mat)]

if (any(rowSums(!is.na(heat_mat)) == 0)) {
  stop("At least one prefoldin row has all NA centered values.")
}

# ---------------------------------------------------------
# 7. Row annotation
# ---------------------------------------------------------

row_annot <- plot_df %>%
  transmute(
    DisplayGene = .data$DisplayGene,
    `100uM vs Veh` = if_else(.data[[paste0(acute_contrasts[1], "_sigflag")]], "Yes", "No"),
    `200uM vs 100uM` = if_else(.data[[paste0(acute_contrasts[2], "_sigflag")]], "Yes", "No"),
    `200uM vs Veh` = if_else(.data[[paste0(acute_contrasts[3], "_sigflag")]], "Yes", "No")
  ) %>%
  as.data.frame()

rownames(row_annot) <- row_annot$DisplayGene
row_annot$DisplayGene <- NULL

annotation_colors <- list(
  `100uM vs Veh` = c("Yes" = "#4DAF4A", "No" = "grey90"),
  `200uM vs 100uM` = c("Yes" = "#984EA3", "No" = "grey90"),
  `200uM vs Veh` = c("Yes" = "#E41A1C", "No" = "grey90")
)

# ---------------------------------------------------------
# 8. Match HeatMap_Fig1.R color scale
# ---------------------------------------------------------

max_abs <- max(abs(heat_mat), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) {
  max_abs <- 1
}

my_breaks <- seq(-max_abs, max_abs, length.out = 101)
my_colors <- colorRampPalette(c("#2400D8", "#FFFFEA", "#A50021"))(100)

# ---------------------------------------------------------
# 9. Draw and save heatmap
# ---------------------------------------------------------

draw_heatmap <- function() {
  pheatmap(
    mat = heat_mat,
    annotation_row = row_annot,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    labels_row = rownames(heat_mat),
    fontsize_row = 10,
    fontsize_col = 12,
    border_color = "grey80",
    angle_col = 0,
    color = my_colors,
    breaks = my_breaks,
    main = paste0(
      "Acute RDH12 prefoldin proteins\n",
      "Centered intensities"
    )
  )
}

draw_heatmap()

pdf(out_file_pdf, width = 7, height = 4.5)
draw_heatmap()
dev.off()

svg(out_file_svg, width = 7, height = 4.5)
draw_heatmap()
dev.off()

# ---------------------------------------------------------
# 10. Export input table in display order
# ---------------------------------------------------------

plot_df_export <- plot_df %>%
  mutate(
    `100uM vs Veh significant` = .data[[paste0(acute_contrasts[1], "_sigflag")]],
    `200uM vs 100uM significant` = .data[[paste0(acute_contrasts[2], "_sigflag")]],
    `200uM vs Veh significant` = .data[[paste0(acute_contrasts[3], "_sigflag")]]
  ) %>%
  select(
    DisplayGene,
    Gene,
    name,
    ID,
    SigIn,
    `100uM vs Veh significant`,
    `200uM vs 100uM significant`,
    `200uM vs Veh significant`,
    all_of(acute_centered_cols),
    all_of(sig_required_cols)
  )

write.csv(plot_df_export, file = out_file_csv, row.names = FALSE)

# ---------------------------------------------------------
# 11. Completion summary
# ---------------------------------------------------------

cat("Prefoldin proteins requested:", nrow(prefoldin_tbl), "\n")
cat("Prefoldin proteins plotted:", nrow(heat_mat), "\n")
cat("PDF written to:", out_file_pdf, "\n")
cat("SVG written to:", out_file_svg, "\n")
cat("CSV written to:", out_file_csv, "\n")

plot_df_export %>%
  select(DisplayGene, Gene, ID, SigIn) %>%
  print(n = Inf)
