# =========================================================
# Main-text acute RDH12 pathway heatmap
# - Treatments are rows: Veh, 100, 200
# - Proteins are columns, labeled by gene ID
# - Protein pathway/group is shown as a column annotation strip
# - Keeps proteins significant in at least one acute contrast
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
})

# ---------------------------------------------------------
# 1. User inputs
# ---------------------------------------------------------

input_file <- "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv"
out_dir <- "Proteomic_Figs/enriched"

out_file_pdf <- file.path(out_dir, "Fig1_horizontal_maintext_heatmap_centered_sigOnly.pdf")
out_file_svg <- file.path(out_dir, "Fig1_horizontal_maintext_heatmap_centered_sigOnly.svg")
out_file_csv <- file.path(out_dir, "Fig1_horizontal_maintext_heatmap_input_table.csv")

padj_cutoff <- 0.1
lfc_cutoff <- log2(1.4)

acute_centered_cols <- c(
  "RDH12_control_atRAL5hr_centered",
  "RDH12_100_atRAL5hr_centered",
  "RDH12_200_atRAL5hr_centered"
)

condition_labels <- c(
  "RDH12_control_atRAL5hr_centered" = "Veh",
  "RDH12_100_atRAL5hr_centered" = "100",
  "RDH12_200_atRAL5hr_centered" = "200"
)

acute_contrasts <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr"
)

acute_contrast_labels <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr" = "100 vs Veh",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr" = "200 vs 100",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr" = "200 vs Veh"
)

cluster_within_group <- TRUE
row_font_size <- 10
col_font_size <- 7
figure_width <- 6.9
figure_height <- 3.8

# ---------------------------------------------------------
# 2. Pathway / protein lists
# ---------------------------------------------------------

Nglycan_prcrsr_biosynthsis <- c(
  "ALG3", "DPAGT1", "ALG1", "ALG2", "ALG11", "ALG5", "ALG8",
  "ALG9", "RFT1", "DPM1", "RPN1", "RPN2", "STT3A", "STT3B",
  "UGGT1"
)

CCT_TriC_genes <- c(
  "TUBB4B", "TUBA1C", "TCP1", "CCT2", "TUBB2A", "TUBAL3",
  "CCT7", "CCT3", "CCT8", "TUBB4A", "TUBB6", "CCT6A",
  "CCT4", "TUBA4A", "CCT5"
)

MetalIon_SLCtransporters <- c(
  "SLC39A7", "SLC39A1", "SLC31A1", "SLC9A6", "SLC39A10",
  "SLC9A1", "SLC39A14", "SLC30A5", "SLC11A2", "SLC7A11"
)

Calcium_handling <- c("ORAI1", "ATP2A2", "ATP2B1")

Prefoldin <- c("PFDN1", "PFDN2", "VBP1", "PFDN4", "PFDN5", "PFDN6")

Synthesis_PA <- c(
  "LPCAT1", "PLA2G4A", "LPCAT4", "AGPAT3", "DDHD1", "GNPAT",
  "MIGA1", "AGPAT4", "GPD2", "GPAT4", "AGPAT5", "DGAT1"
)

pathway_tbl_all <- bind_rows(
  tibble(Gene = Nglycan_prcrsr_biosynthsis, Pathway = "N-glycan"),
  tibble(Gene = CCT_TriC_genes, Pathway = "CCT/TriC"),
  tibble(Gene = MetalIon_SLCtransporters, Pathway = "Transporters"),
  tibble(Gene = Calcium_handling, Pathway = "Calcium handling"),
  tibble(Gene = Prefoldin, Pathway = "Prefoldin"),
  tibble(Gene = Synthesis_PA, Pathway = "Lipid Metabolism")
)

pathway_order <- c(
  "N-glycan",
  "CCT/TriC",
  "Transporters",
  "Calcium handling",
  "Prefoldin",
  "Lipid Metabolism"
)

pathway_tbl <- pathway_tbl_all %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(
    Pathway = factor(Pathway, levels = pathway_order),
    DisplayGene = if_else(Gene == "VBP1", "VBP1/PFDN3", Gene)
  )

dup_genes <- pathway_tbl_all %>%
  count(Gene) %>%
  filter(n > 1)

if (nrow(dup_genes) > 0) {
  message("Genes listed in more than one pathway were assigned to the first pathway:")
  print(dup_genes)
}

# ---------------------------------------------------------
# 3. Import data and check columns
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
# 4. Apply acute significance filter
# ---------------------------------------------------------

DEPresults_sig <- DEPresults_v2 %>%
  mutate(
    GeneClean = coalesce(
      na_if(as.character(.data$Gene), ""),
      na_if(as.character(.data$name), "")
    ),
    across(all_of(c(acute_centered_cols, sig_required_cols)), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  filter(!is.na(.data$GeneClean), .data$GeneClean != "") %>%
  distinct(.data$GeneClean, .keep_all = TRUE)

for (ct in acute_contrasts) {
  padj_col <- paste0(ct, "_p.adj")
  ratio_col <- paste0(ct, "_ratio")
  sig_col <- paste0(ct, "_sigflag")

  DEPresults_sig[[sig_col]] <-
    !is.na(DEPresults_sig[[padj_col]]) &
    !is.na(DEPresults_sig[[ratio_col]]) &
    DEPresults_sig[[padj_col]] < padj_cutoff &
    abs(DEPresults_sig[[ratio_col]]) >= lfc_cutoff
}

sigflag_cols <- paste0(acute_contrasts, "_sigflag")

DEPresults_sig <- DEPresults_sig %>%
  mutate(
    keep_sig_any = if_any(all_of(sigflag_cols), identity),
    SigIn = apply(
      dplyr::select(., all_of(sigflag_cols)),
      1,
      function(x) {
        hit_names <- acute_contrast_labels[acute_contrasts[as.logical(x)]]
        if (length(hit_names) == 0) "None" else paste(hit_names, collapse = "; ")
      }
    )
  )

plot_df <- pathway_tbl %>%
  left_join(DEPresults_sig, by = c("Gene" = "GeneClean")) %>%
  filter(.data$keep_sig_any) %>%
  filter(rowSums(!is.na(across(all_of(acute_centered_cols)))) > 0)

if (nrow(plot_df) == 0) {
  stop("No proteins from the pathway lists passed the acute significance filter.")
}

# ---------------------------------------------------------
# 5. Order proteins by pathway, then by similar treatment pattern
# ---------------------------------------------------------

cluster_block <- function(block_df) {
  if (!cluster_within_group || nrow(block_df) <= 1) {
    return(block_df$Gene)
  }

  mat_block <- block_df %>%
    dplyr::select(all_of(acute_centered_cols)) %>%
    as.matrix()

  rownames(mat_block) <- block_df$Gene

  for (i in seq_len(nrow(mat_block))) {
    row_vals <- mat_block[i, ]
    if (anyNA(row_vals)) {
      row_mean <- mean(row_vals, na.rm = TRUE)
      if (is.nan(row_mean)) row_mean <- 0
      row_vals[is.na(row_vals)] <- row_mean
      mat_block[i, ] <- row_vals
    }
  }

  rownames(mat_block)[hclust(dist(mat_block), method = "complete")$order]
}

ordered_genes <- character(0)

for (pw in pathway_order) {
  block_df <- plot_df %>%
    filter(.data$Pathway == pw)

  if (nrow(block_df) == 0) next

  ordered_genes <- c(ordered_genes, cluster_block(block_df))
}

plot_df <- plot_df %>%
  mutate(Gene = factor(.data$Gene, levels = ordered_genes)) %>%
  arrange(.data$Gene)

# ---------------------------------------------------------
# 6. Build treatment-by-protein heatmap matrix
# ---------------------------------------------------------

protein_mat <- plot_df %>%
  dplyr::select(DisplayGene, all_of(acute_centered_cols)) %>%
  as.data.frame()

rownames(protein_mat) <- protein_mat$DisplayGene
protein_mat$DisplayGene <- NULL
protein_mat <- as.matrix(protein_mat)

heat_mat <- t(protein_mat)
rownames(heat_mat) <- unname(condition_labels[rownames(heat_mat)])

annotation_col <- plot_df %>%
  transmute(DisplayGene, Pathway = factor(Pathway, levels = pathway_order)) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$DisplayGene
annotation_col$DisplayGene <- NULL

pathway_counts <- table(factor(as.character(annotation_col$Pathway), levels = pathway_order))
pathway_counts <- pathway_counts[pathway_counts > 0]
gaps_col <- cumsum(pathway_counts)
gaps_col <- gaps_col[-length(gaps_col)]

# ---------------------------------------------------------
# 7. Colors
# ---------------------------------------------------------

annotation_colors <- list(
  Pathway = c(
    "N-glycan" = "#E69F00",
    "CCT/TriC" = "#56B4E9",
    "Transporters" = "#009E73",
    "Calcium handling" = "#0072B2",
    "Prefoldin" = "#CC79A7",
    "Lipid Metabolism" = "#F0E442"
  )
)

max_abs <- max(abs(heat_mat), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) {
  max_abs <- 1
}

my_breaks <- seq(-max_abs, max_abs, length.out = 101)
my_colors <- colorRampPalette(c("#2400D8", "#FFFFEA", "#A50021"))(100)

# ---------------------------------------------------------
# 8. Draw and save heatmap
# ---------------------------------------------------------

draw_heatmap <- function() {
  pheatmap(
    mat = heat_mat,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_col = gaps_col,
    labels_col = colnames(heat_mat),
    fontsize_row = row_font_size,
    fontsize_col = col_font_size,
    border_color = "grey80",
    angle_col = 90,
    color = my_colors,
    breaks = my_breaks,
    legend = TRUE,
    main = "Acute RDH12 pathway proteins"
  )
}

draw_heatmap()

pdf(out_file_pdf, width = figure_width, height = figure_height)
draw_heatmap()
dev.off()

svg(out_file_svg, width = figure_width, height = figure_height)
draw_heatmap()
dev.off()

# ---------------------------------------------------------
# 9. Export plotted input table in display order
# ---------------------------------------------------------

plot_df_export <- plot_df %>%
  mutate(
    `100 vs Veh significant` = .data[[paste0(acute_contrasts[1], "_sigflag")]],
    `200 vs 100 significant` = .data[[paste0(acute_contrasts[2], "_sigflag")]],
    `200 vs Veh significant` = .data[[paste0(acute_contrasts[3], "_sigflag")]]
  ) %>%
  dplyr::select(
    Gene,
    DisplayGene,
    Pathway,
    name,
    ID,
    SigIn,
    `100 vs Veh significant`,
    `200 vs 100 significant`,
    `200 vs Veh significant`,
    all_of(acute_centered_cols),
    all_of(sig_required_cols)
  )

write.csv(plot_df_export, file = out_file_csv, row.names = FALSE)

# ---------------------------------------------------------
# 10. Completion checks
# ---------------------------------------------------------

cat("Proteins requested in pathway table:", nrow(pathway_tbl), "\n")
cat("Proteins plotted:", ncol(heat_mat), "\n")
cat("Heatmap rows:", paste(rownames(heat_mat), collapse = ", "), "\n")
cat("Heatmap columns include VBP1/PFDN3:", "VBP1/PFDN3" %in% colnames(heat_mat), "\n")
cat("PDF written to:", out_file_pdf, "\n")
cat("SVG written to:", out_file_svg, "\n")
cat("CSV written to:", out_file_csv, "\n\n")

plot_df_export %>%
  count(Pathway, sort = FALSE) %>%
  print()
