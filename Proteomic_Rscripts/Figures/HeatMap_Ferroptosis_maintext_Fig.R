# =========================================================
# Ferroptosis-associated proteomic heatmap
#
# Displays RDH12 acute + RDH12 recovery centered intensities
# Breaks ferroptosis-associated proteins into functional modules:
#   1. KEAP1-NRF2 oxidative stress axis
#   2. Cystine import / glutathione synthesis
#   3. GPX / thioredoxin / peroxide detoxification
#   4. GST / MGST lipid electrophile detoxification
#   5. Iron import / storage / export / chaperoning
#   6. Heme / metal / copper redox handling
#   7. Lipid activation / PUFA remodeling / lipoxygenase axis
#   8. Mitochondrial ROS / VDAC / respiration-associated stress
#   9. Autophagy / ferritinophagy-related stress
#   10. Stress/apoptosis-adjacent regulators
#
# Uses color scale:
#   c("#2400D8", "#FFFFEA", "#A50021")
#
# By default, this script keeps ALL genes from the ferroptosis-related lists
# that are found in the DEP table, regardless of significance.
#
# Optional: set keep_only_significant <- TRUE to keep only proteins
# significant in at least one acute contrast.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(pheatmap)
})

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------

DEPresults_v2 <- read.csv(
  "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

out_dir <- "Proteomic_Figs/ferroptosis_heatmap/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file_pdf <- file.path(out_dir, "1Ferroptosis_process_heatmap_RDH12_acute_recovery.pdf")
out_file_svg <- file.path(out_dir, "1Ferroptosis_process_heatmap_RDH12_acute_recovery.svg")
out_file_png <- file.path(out_dir, "1Ferroptosis_process_heatmap_RDH12_acute_recovery.png")
out_file_csv <- file.path(out_dir, "Ferroptosis_process_heatmap_input_table.csv")

# Significance thresholds
padj_cutoff <- 0.1
lfc_cutoff  <- log2(1.4)

# FALSE = keep all ferroptosis-related proteins found in DEP table
# TRUE  = keep only proteins significant in at least one acute contrast
keep_only_significant <- FALSE

cluster_cols <- FALSE
cluster_rows_within_process <- TRUE
remove_all_na_rows <- TRUE

row_distance_method <- "euclidean"
row_clustering_method <- "complete"

# FALSE = use centered intensities directly
# TRUE  = z-score each protein across displayed columns
scale_rows <- FALSE

# ---------------------------------------------------------
# 2. DISPLAY COLUMNS
# ---------------------------------------------------------

display_centered_cols <- c(
  "RDH12_control_atRAL5hr_centered",
  "RDH12_100_atRAL5hr_centered",
  "RDH12_200_atRAL5hr_centered",
  "RDH12_control_atRAL5hr.24h_recvr_centered",
  "RDH12_100_atRAL5hr.24h_recvr_centered",
  "RDH12_200_atRAL5hr.24h_recvr_centered"
)

display_col_labels <- c(
  "RDH12_control_atRAL5hr_centered"            = "Acute Veh",
  "RDH12_100_atRAL5hr_centered"                = "Acute 100 µM",
  "RDH12_200_atRAL5hr_centered"                = "Acute 200 µM",
  "RDH12_control_atRAL5hr.24h_recvr_centered"  = "Recovery Veh",
  "RDH12_100_atRAL5hr.24h_recvr_centered"      = "Recovery 100 µM",
  "RDH12_200_atRAL5hr.24h_recvr_centered"      = "Recovery 200 µM"
)

# ---------------------------------------------------------
# 3. ACUTE CONTRASTS FOR OPTIONAL SIGNIFICANCE FILTERING
# ---------------------------------------------------------

acute_contrasts <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr"
)

acute_contrast_labels <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr" = "100 µM vs Veh",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr"     = "200 µM vs 100 µM",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr" = "200 µM vs Veh"
)

sig_required_cols <- unlist(lapply(
  acute_contrasts,
  function(x) c(paste0(x, "_p.adj"), paste0(x, "_ratio"))
))

required_cols <- c(
  "Gene",
  "name",
  "ID",
  display_centered_cols,
  sig_required_cols
)

missing_cols <- setdiff(required_cols, colnames(DEPresults_v2))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from DEPresults_v2:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 4. FERROPTOSIS-ASSOCIATED FUNCTIONAL GROUPS
# ---------------------------------------------------------
# These groups intentionally separate proteins into mechanistic roles.
# Some proteins could plausibly belong to multiple categories.
# To avoid duplicated rows, genes are assigned to the first category
# where they appear.

keap1_nrf2_axis <- c(
  # NRF2 / electrophilic-stress response
  "NFE2L2", "KEAP1", "HMOX1",
  
  # Cystine uptake and glutathione synthesis / recycling
  "SLC7A11", "GCLC", "GCLM", "GSS", "GSR",
  
  # Peroxide detoxification capacity
  "GPX4", "GPX1", "TXNRD1", "TXNRD2", "PRDX6",
  
  # Lipid-electrophile detoxification
  "MGST2", "MGST3", "GSTA4", "GSTZ1",
  
  # Iron uptake and storage
  "TFRC", "SLC11A2", "FTH1", "FTL",
  
  # Lipid remodeling / ROS context
  "ACSL4", "ACSL3", "CYBA"
)

cystine_gsh_synthesis <- c(
  "SLC7A"
)

gpx_txn_prdx_peroxide_detox <- c(
  "GP1" 
 
)

gst_mgst_lipid_electrophile_detox <- c(
  "GSTO2")

iron_import_storage_export <- c(
  "TF"
)

heme_metal_copper_redox <- c(
  "HMO1")
 

lipid_activation_peroxidation_axis <- c(
  "ACSL1"
)

mitochondrial_ros_vdac_respiration <- c(
  "VDAC2"
)

autophagy_ferritinophagy_related <- c(
  "ATG7"
  
)

stress_apoptosis_adjacent <- c(
  "TP53"
  
)

# ---------------------------------------------------------
# 5. BUILD FUNCTIONAL PROCESS TABLE
# ---------------------------------------------------------

process_tbl_all <- bind_rows(
  tibble(Gene = keap1_nrf2_axis, Process = "KEAP1-NRF2 oxidative stress axis"),
  tibble(Gene = cystine_gsh_synthesis, Process = "Cystine import / glutathione synthesis"),
  tibble(Gene = gpx_txn_prdx_peroxide_detox, Process = "GPX / TXN / PRDX peroxide detoxification"),
  tibble(Gene = gst_mgst_lipid_electrophile_detox, Process = "GST / MGST lipid-electrophile detoxification"),
  tibble(Gene = iron_import_storage_export, Process = "Iron import / storage / export"),
  tibble(Gene = heme_metal_copper_redox, Process = "Heme / metal / copper redox handling"),
  tibble(Gene = lipid_activation_peroxidation_axis, Process = "Lipid activation / peroxidation axis"),
  tibble(Gene = mitochondrial_ros_vdac_respiration, Process = "Mitochondrial ROS / VDAC stress"),
  tibble(Gene = autophagy_ferritinophagy_related, Process = "Autophagy / ferritinophagy-related"),
  tibble(Gene = stress_apoptosis_adjacent, Process = "Stress / apoptosis-adjacent regulators")
) %>%
  mutate(
    Gene = str_trim(Gene),
    Process = str_trim(Process)
  ) %>%
  filter(Gene != "")

# Optional diagnostic: genes present in more than one functional group
dup_genes <- process_tbl_all %>%
  count(Gene, sort = TRUE) %>%
  filter(n > 1)

if (nrow(dup_genes) > 0) {
  message("\nGenes appearing in more than one functional group; assigned to first group only:\n")
  print(dup_genes)
}

process_tbl <- process_tbl_all %>%
  distinct(Gene, .keep_all = TRUE)

process_order <- c(
  "KEAP1-NRF2 oxidative stress axis",
  "Cystine import / glutathione synthesis",
  "GPX / TXN / PRDX peroxide detoxification",
  "GST / MGST lipid-electrophile detoxification",
  "Iron import / storage / export",
  "Heme / metal / copper redox handling",
  "Lipid activation / peroxidation axis",
  "Mitochondrial ROS / VDAC stress",
  "Autophagy / ferritinophagy-related",
  "Stress / apoptosis-adjacent regulators"
)

process_tbl <- process_tbl %>%
  mutate(Process = factor(Process, levels = process_order))

# ---------------------------------------------------------
# 6. FLAG ACUTE SIGNIFICANCE
# ---------------------------------------------------------

DEPresults_sig <- DEPresults_v2 %>%
  mutate(
    Gene = dplyr::coalesce(na_if(Gene, ""), na_if(name, ""))
  ) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Gene, .keep_all = TRUE)

for (ct in acute_contrasts) {
  padj_col <- paste0(ct, "_p.adj")
  ratio_col <- paste0(ct, "_ratio")
  sig_col <- paste0(ct, "_sigflag")
  
  DEPresults_sig[[sig_col]] <-
    !is.na(DEPresults_sig[[padj_col]]) &
    !is.na(DEPresults_sig[[ratio_col]]) &
    DEPresults_sig[[padj_col]] <= padj_cutoff &
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
        if (length(hit_names) == 0) {
          "Not significant in acute contrasts"
        } else {
          paste(hit_names, collapse = "; ")
        }
      }
    )
  )

# ---------------------------------------------------------
# 7. SUBSET TO FERROPTOSIS-ASSOCIATED PROTEINS
# ---------------------------------------------------------

plot_df <- process_tbl %>%
  left_join(DEPresults_sig, by = "Gene")

# Keep all requested proteins found in the table
plot_df <- plot_df %>%
  filter(!is.na(ID) | !is.na(name))

# Optional significance filter
if (keep_only_significant) {
  plot_df <- plot_df %>%
    filter(keep_sig_any)
}

if (nrow(plot_df) == 0) {
  stop("No ferroptosis-associated proteins from the process lists were found in DEPresults_v2.")
}

# ---------------------------------------------------------
# 8. BUILD HEATMAP MATRIX
# ---------------------------------------------------------

heat_df <- plot_df %>%
  dplyr::select(
    Gene,
    Process,
    SigIn,
    all_of(display_centered_cols),
    all_of(sigflag_cols),
    all_of(sig_required_cols)
  )

if (remove_all_na_rows) {
  heat_df <- heat_df %>%
    filter(rowSums(!is.na(across(all_of(display_centered_cols)))) > 0)
}

if (nrow(heat_df) == 0) {
  stop("All ferroptosis-associated rows were removed after NA filtering.")
}

heat_mat <- heat_df %>%
  dplyr::select(all_of(display_centered_cols)) %>%
  as.data.frame()

rownames(heat_mat) <- heat_df$Gene
heat_mat <- as.matrix(heat_mat)

unmapped_cols <- setdiff(colnames(heat_mat), names(display_col_labels))
if (length(unmapped_cols) > 0) {
  stop(
    "These heatmap columns do not have labels in display_col_labels:\n",
    paste(unmapped_cols, collapse = "\n")
  )
}

colnames(heat_mat) <- unname(display_col_labels[colnames(heat_mat)])

# ---------------------------------------------------------
# 9. ROW ANNOTATIONS
# ---------------------------------------------------------

row_annot <- heat_df %>%
  transmute(
    Gene = Gene,
    Process = Process,
    `100 µM vs Veh` = if_else(.data[[paste0(acute_contrasts[1], "_sigflag")]], "Yes", "No"),
    `200 µM vs 100 µM` = if_else(.data[[paste0(acute_contrasts[2], "_sigflag")]], "Yes", "No"),
    `200 µM vs Veh` = if_else(.data[[paste0(acute_contrasts[3], "_sigflag")]], "Yes", "No")
  ) %>%
  distinct() %>%
  as.data.frame()

rownames(row_annot) <- row_annot$Gene
row_annot$Gene <- NULL

row_annot$Process <- factor(row_annot$Process, levels = process_order)

# ---------------------------------------------------------
# 10. OPTIONAL ROW SCALING
# ---------------------------------------------------------

if (scale_rows) {
  heat_mat <- t(scale(t(heat_mat)))
}

# ---------------------------------------------------------
# 11. CLUSTER WITHIN EACH PROCESS
# ---------------------------------------------------------

cluster_within_group <- function(mat_block,
                                 distance_method = "euclidean",
                                 clustering_method = "complete") {
  if (nrow(mat_block) <= 1) {
    return(rownames(mat_block))
  }
  
  mat_for_clust <- mat_block
  
  for (i in seq_len(nrow(mat_for_clust))) {
    row_vals <- mat_for_clust[i, ]
    
    if (anyNA(row_vals)) {
      row_mean <- mean(row_vals, na.rm = TRUE)
      if (is.nan(row_mean)) row_mean <- 0
      row_vals[is.na(row_vals)] <- row_mean
      mat_for_clust[i, ] <- row_vals
    }
  }
  
  hc <- hclust(
    dist(mat_for_clust, method = distance_method),
    method = clustering_method
  )
  
  rownames(mat_block)[hc$order]
}

ordered_genes <- character(0)

for (proc in process_order) {
  genes_in_proc <- rownames(row_annot)[row_annot$Process == proc]
  
  if (length(genes_in_proc) == 0) next
  
  mat_block <- heat_mat[genes_in_proc, , drop = FALSE]
  
  if (cluster_rows_within_process) {
    block_order <- cluster_within_group(
      mat_block = mat_block,
      distance_method = row_distance_method,
      clustering_method = row_clustering_method
    )
  } else {
    block_order <- rownames(mat_block)
  }
  
  ordered_genes <- c(ordered_genes, block_order)
}

heat_mat <- heat_mat[ordered_genes, , drop = FALSE]
row_annot <- row_annot[ordered_genes, , drop = FALSE]
display_labels <- ordered_genes

# ---------------------------------------------------------
# 12. GAPS BETWEEN PROCESS BLOCKS
# ---------------------------------------------------------

process_counts <- table(factor(row_annot$Process, levels = process_order))
process_counts <- process_counts[process_counts > 0]

gaps_row <- cumsum(process_counts)
gaps_row <- gaps_row[-length(gaps_row)]

# ---------------------------------------------------------
# 13. COLORS
# ---------------------------------------------------------

process_colors <- c(
  "KEAP1-NRF2 oxidative stress axis" = "#004949",
  "Cystine import / glutathione synthesis" = "#1B9E77",
  "GPX / TXN / PRDX peroxide detoxification" = "#7570B3",
  "GST / MGST lipid-electrophile detoxification" = "#E7298A",
  "Iron import / storage / export" = "#D95F02",
  "Heme / metal / copper redox handling" = "#A6761D",
  "Lipid activation / peroxidation axis" = "#66A61E",
  "Mitochondrial ROS / VDAC stress" = "#E6AB02",
  "Autophagy / ferritinophagy-related" = "#666666",
  "Stress / apoptosis-adjacent regulators" = "#A6CEE3"
)

sig_colors <- c(
  "Yes" = "#D73027",
  "No" = "grey90"
)

annotation_colors <- list(
  Process = process_colors,
  `100 µM vs Veh` = sig_colors,
  `200 µM vs 100 µM` = sig_colors,
  `200 µM vs Veh` = sig_colors
)

# User-requested heatmap color scale
max_abs <- max(abs(heat_mat), na.rm = TRUE)

my_breaks <- seq(-max_abs, max_abs, length.out = 101)
my_colors <- colorRampPalette(c("#2400D8", "#FFFFEA", "#A50021"))(100)

legend_breaks <- c(-max_abs, -max_abs / 2, 0, max_abs / 2, max_abs)
legend_labels <- round(legend_breaks, 2)

# ---------------------------------------------------------
# 14. DRAW HEATMAP
# ---------------------------------------------------------

draw_heatmap <- function() {
  pheatmap(
    mat = heat_mat,
    annotation_row = row_annot,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    gaps_row = gaps_row,
    labels_row = display_labels,
    fontsize_row = 9,
    fontsize_col = 11,
    border_color = "grey85",
    angle_col = "45",
    color = my_colors,
    breaks = my_breaks,
    legend_breaks = legend_breaks,
    legend_labels = legend_labels,
    main = paste0(
      "Ferroptosis-associated proteomic response in RDH12 cells\n",
      "Functional modules shown across acute and recovery atRAL exposure"
    )
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 15. SAVE HEATMAP
# ---------------------------------------------------------

pdf(out_file_pdf, width = 11, height = 12)
draw_heatmap()
dev.off()

svg(out_file_svg, width = 11, height = 12)
draw_heatmap()
dev.off()

png(out_file_png, width = 3300, height = 3600, res = 300)
draw_heatmap()
dev.off()

# ---------------------------------------------------------
# 16. EXPORT INPUT TABLE IN DISPLAY ORDER
# ---------------------------------------------------------

plot_df_export <- heat_df %>%
  mutate(
    Gene = factor(Gene, levels = ordered_genes)
  ) %>%
  arrange(Gene) %>%
  dplyr::select(
    Gene,
    Process,
    SigIn,
    all_of(display_centered_cols),
    all_of(sig_required_cols)
  )

write.csv(
  plot_df_export,
  file = out_file_csv,
  row.names = FALSE
)

# ---------------------------------------------------------
# 17. QUICK CHECKS
# ---------------------------------------------------------

cat("\nRequested ferroptosis-associated genes:", nrow(process_tbl), "\n")
cat("Genes found and displayed:", nrow(heat_mat), "\n")
cat("Columns displayed:", ncol(heat_mat), "\n")

cat("\nHeatmap column labels:\n")
print(colnames(heat_mat))

cat("\nProteins per functional process shown:\n")
print(table(row_annot$Process))

cat("\nProteins shown and functional assignments:\n")
heat_df %>%
  dplyr::select(Gene, Process, SigIn) %>%
  arrange(Process, Gene) %>%
  print(n = 300)

cat("\nSaved files:\n")
cat(out_file_pdf, "\n")
cat(out_file_svg, "\n")
cat(out_file_png, "\n")
cat(out_file_csv, "\n")