# =========================================================
# KEAP1/NRF2 pathway heatmap from DEPresults_v2 centered intensities
# Displays acute + recovery centered intensities
# Keeps only proteins significant in at least one acute contrast
# Clusters similar-behaving proteins WITHIN each pathway
# Adds row annotations for which acute contrast(s) each protein is significant in
# =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(stringr)
library(forcats)
library(tibble)

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------
DEPresults_v2 <- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")
colnames(DEPresults_v2)

# Columns to display in the heatmap
# FIXED: removed duplicated acute columns
display_centered_cols <- c(
  "RDH12_control_atRAL5hr_centered",
  "RDH12_100_atRAL5hr_centered",
  "RDH12_200_atRAL5hr_centered",
  "RDH12_control_atRAL5hr.24h_recvr_centered",
  "RDH12_100_atRAL5hr.24h_recvr_centered",
  "RDH12_200_atRAL5hr.24h_recvr_centered"
)

# Acute contrasts used to decide whether a protein is kept
acute_contrasts <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr"
)

# Pretty names for those contrasts
acute_contrast_labels <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr" = "100uM vs Veh",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr"     = "200uM vs 100uM",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr" = "200uM vs Veh"
)

# Thresholds for filtering proteins
padj_cutoff <- 0.1
lfc_cutoff  <- log2(1.4)

# Pretty labels for heatmap columns
# FIXED: recovery names now exactly match the actual *_centered column names
display_col_labels <- c(
  "RDH12_control_atRAL5hr_centered"            = "Acute Veh",
  "RDH12_100_atRAL5hr_centered"                = "Acute 100uM",
  "RDH12_200_atRAL5hr_centered"                = "Acute 200uM",
  "RDH12_control_atRAL5hr.24h_recvr_centered" = "Recvr Veh",
  "RDH12_100_atRAL5hr.24h_recvr_centered"     = "Recvr 100uM",
  "RDH12_200_atRAL5hr.24h_recvr_centered"     = "Recvr 200uM"
)

# Output files
out_file_pdf <- "Proteomic_Figs/enriched/KEAP1_NRF2_pathway_heatmap_centered_sigOnly_clusteredWithinPathway_wContrastAnno.pdf"
out_file_svg <- "Proteomic_Figs/enriched/KEAP1_NRF2_pathway_heatmap_centered_sigOnly_clusteredWithinPathway_wContrastAnno.svg"
out_file_csv <- "Proteomic_Figs/enriched/KEAP1_NRF2_pathway_heatmap_input_table_sigOnly_clusteredWithinPathway_wContrastAnno.csv"

# Cluster columns?
cluster_cols <- FALSE

# Remove genes with all NA across displayed columns
remove_all_na_rows <- TRUE

# Distance and linkage for within-pathway clustering
row_distance_method <- "euclidean"
row_clustering_method <- "complete"

# Scale rows before clustering/plotting?
# FALSE = use centered intensities directly
# TRUE  = z-score each row across displayed columns
scale_rows <- FALSE

# ---------------------------------------------------------
# 2. PATHWAY / GENE LISTS
# ---------------------------------------------------------

KEAP1Nrf2 <- c(
  "KEAP1", "NRF2", "NFE2L2", "NQO1", "HMOX1", "GSTO2", "GSTA5", "GSTA4",
  "GSTCD", "GSTK1", "GSTM2", "GSTM3", "GSTO1", "GSTZ1", "MGST1", "MGST2", "MGST3",
  "TXNRD1", "TXNRD2", "AKR1A1", "AKR1B1", "AKR1E2", "AKR7A2","CYBA",
  "SOD2", "MET", "PTPN1", "MET", "GPX1", "GPX3", "GPX4", "GPX8", "CUL1", "ATOX1", "CCS", "GSR", 
  "AQP8", "CYBB", "CYBC", "P4HB", "TXN", "NCF1", "NCF2", "PRDX6", "PRDX5", "PRDX3", "PRDX2", "NUDT2", "CYCS", 
  "ATP7A", "PRDX1", "NCF4", "ERO1A", "NOX5", "TXN2","CYP2S1","CYP51A1","CYP20A1","ESD","RFESD","ACY1","GSS","GCLC",
  "GCLM", "BACH1", "BACH2", "CHAC2", "CNDP2", "CHAC1","HIGD1A", "HDAC3", "MAFK", "TBL1X", "SIN3B", "NCOR1", "MT-CO2", 
  "ALB", "FABP1", "RXRA", "BLVRB", "HMOX2", "ABCC1", "STAT3", "BLVRA", "HBB",
  "HBA1", "CYCS", "PPARA", "PTK6", "NCOA6", "NCOA2", "MED1", "CHD9", "SMARCD3", "CARM1", "HM13", "CREBBP",
  "NLRP3","TGS1", "SIN3A", "HELZ2","TBL1XR1","TXNIP", "STAP2", "NCOR2", "ACOX1","ACOX2","COX11","COX14",
  "COX15", "COX16","COX17","COX19","COX20","COX4I1","COX5A","COX5B", "COX6A1","COX6B1","COX6C",
  "COX7A2","COX7A2L","COX7C", "NCOA4","UBE2Q1", "NEDD4L", "TCP11L1","YAP1", "ATG4","TAX1BP1","TAX1BP3", "ATG4", "ATG7",
  "ATG3","ATG9A","HPCAL1","HPCAL4","CDH2","SFXN1"
  
)

# ---------------------------------------------------------
# 3. BUILD PATHWAY TABLE
# ---------------------------------------------------------

pathway_tbl <- bind_rows(
  tibble(
    Gene = KEAP1Nrf2,
    Pathway = "KEAP1Nrf2"
  )
) %>%
  distinct(Gene, .keep_all = TRUE)

pathway_order <- c("KEAP1Nrf2")

pathway_tbl <- pathway_tbl %>%
  mutate(Pathway = factor(Pathway, levels = pathway_order))

# ---------------------------------------------------------
# 4. CHECK REQUIRED COLUMNS
# ---------------------------------------------------------

sig_required_cols <- unlist(lapply(
  acute_contrasts,
  function(x) c(paste0(x, "_p.adj"), paste0(x, "_ratio"))
))

required_cols <- c("Gene", "name", "ID", display_centered_cols, sig_required_cols)

missing_cols <- setdiff(required_cols, colnames(DEPresults_v2))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from DEPresults_v2:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 5. FLAG PROTEINS SIGNIFICANT IN >=1 ACUTE CONTRAST
# ---------------------------------------------------------

DEPresults_sig <- DEPresults_v2 %>%
  mutate(
    Gene = dplyr::coalesce(Gene, name)
  )

# Create one TRUE/FALSE flag per contrast
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
# 6. SUBSET TO GENES OF INTEREST AND SIG-ONLY
# ---------------------------------------------------------

plot_df <- DEPresults_sig %>%
  filter(keep_sig_any) %>%
  inner_join(pathway_tbl, by = "Gene") %>%
  select(
    Gene,
    name,
    ID,
    Pathway,
    keep_sig_any,
    SigIn,
    all_of(sigflag_cols),
    all_of(display_centered_cols),
    all_of(sig_required_cols)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

if (nrow(plot_df) == 0) {
  stop("No significant genes from pathway_tbl were found in DEPresults_v2.")
}

# ---------------------------------------------------------
# 7. BUILD HEATMAP MATRIX
# ---------------------------------------------------------

heat_mat <- plot_df %>%
  select(Gene, all_of(display_centered_cols)) %>%
  as.data.frame()

rownames(heat_mat) <- heat_mat$Gene
heat_mat$Gene <- NULL
heat_mat <- as.matrix(heat_mat)

# Safety check for unmapped column labels
unmapped_cols <- setdiff(colnames(heat_mat), names(display_col_labels))
if (length(unmapped_cols) > 0) {
  stop(
    "These heatmap columns do not have labels in display_col_labels:\n",
    paste(unmapped_cols, collapse = "\n")
  )
}

colnames(heat_mat) <- unname(display_col_labels[colnames(heat_mat)])

# ---------------------------------------------------------
# 8. ROW ANNOTATION
# ---------------------------------------------------------

row_annot <- plot_df %>%
  transmute(
    Gene = Gene,
    Pathway = Pathway,
    `100uM vs Veh` = if_else(.data[[paste0(acute_contrasts[1], "_sigflag")]], "Yes", "No"),
    `200uM vs 100uM` = if_else(.data[[paste0(acute_contrasts[2], "_sigflag")]], "Yes", "No"),
    `200uM vs Veh` = if_else(.data[[paste0(acute_contrasts[3], "_sigflag")]], "Yes", "No")
  ) %>%
  distinct() %>%
  as.data.frame()

rownames(row_annot) <- row_annot$Gene
row_annot$Gene <- NULL

# ---------------------------------------------------------
# 9. OPTIONAL: REMOVE ROWS WITH ALL NA
# ---------------------------------------------------------

if (remove_all_na_rows) {
  keep_rows <- rowSums(!is.na(heat_mat)) > 0
  heat_mat <- heat_mat[keep_rows, , drop = FALSE]
  row_annot <- row_annot[keep_rows, , drop = FALSE]
}

# ---------------------------------------------------------
# 10. OPTIONAL: SCALE EACH ROW
# ---------------------------------------------------------

if (scale_rows) {
  heat_mat <- t(scale(t(heat_mat)))
}

# ---------------------------------------------------------
# 11. CLUSTER WITHIN EACH PATHWAY
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

ordered_genes <- c()

for (pw in pathway_order) {
  genes_in_pw <- rownames(row_annot)[row_annot$Pathway == pw]
  
  if (length(genes_in_pw) == 0) next
  
  mat_block <- heat_mat[genes_in_pw, , drop = FALSE]
  
  block_order <- cluster_within_group(
    mat_block = mat_block,
    distance_method = row_distance_method,
    clustering_method = row_clustering_method
  )
  
  ordered_genes <- c(ordered_genes, block_order)
}

heat_mat <- heat_mat[ordered_genes, , drop = FALSE]
row_annot <- row_annot[ordered_genes, , drop = FALSE]
display_labels <- ordered_genes

# ---------------------------------------------------------
# 12. GAPS BETWEEN PATHWAY BLOCKS
# ---------------------------------------------------------

pathway_counts <- table(factor(row_annot$Pathway, levels = pathway_order))
pathway_counts <- pathway_counts[pathway_counts > 0]
gaps_row <- cumsum(pathway_counts)
gaps_row <- gaps_row[-length(gaps_row)]

# ---------------------------------------------------------
# 13. COLORS
# ---------------------------------------------------------

pathway_colors <- list(
  Pathway = c(
    "KEAP1Nrf2" = "#004949"
  ),
  `100uM vs Veh` = c("Yes" = "#4DAF4A", "No" = "grey90"),
  `200uM vs 100uM` = c("Yes" = "#984EA3", "No" = "grey90"),
  `200uM vs Veh` = c("Yes" = "#E41A1C", "No" = "grey90")
)

# ---------------------------------------------------------
# 14. DRAW HEATMAP
# ---------------------------------------------------------

draw_heatmap <- function() {
  pheatmap(
    mat = heat_mat,
    annotation_row = row_annot,
    annotation_colors = pathway_colors,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    gaps_row = gaps_row,
    labels_row = display_labels,
    fontsize_row = 12,
    fontsize_col = 12,
    border_color = "grey80",
    angle_col = 0,
    main = paste0(
      "KEAP1/NRF2 pathway proteins in RDH12\n",
      "Shown: acute + recovery centered intensities; kept if significant in >=1 acute contrast"
    )
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 15. SAVE HEATMAP
# ---------------------------------------------------------

pdf(out_file_pdf, width = 9.5, height = 15)
draw_heatmap()
dev.off()

svg(out_file_svg, width = 9.5, height = 15)
draw_heatmap()
dev.off()

# ---------------------------------------------------------
# 16. EXPORT INPUT TABLE IN DISPLAY ORDER
# ---------------------------------------------------------

plot_df_export <- plot_df %>%
  mutate(
    Gene = factor(Gene, levels = ordered_genes)
  ) %>%
  arrange(Gene)

write.csv(
  plot_df_export,
  file = out_file_csv,
  row.names = FALSE
)

# ---------------------------------------------------------
# 17. QUICK CHECKS
# ---------------------------------------------------------

cat("Number of unique genes requested:", nrow(pathway_tbl), "\n")
cat("Number of genes significant in >=1 acute contrast and found in pathway lists:", nrow(plot_df), "\n")
cat("Rows displayed in heatmap:", nrow(heat_mat), "\n")
cat("Columns displayed in heatmap:", ncol(heat_mat), "\n")
cat("Heatmap column labels:", paste(colnames(heat_mat), collapse = ", "), "\n\n")

plot_df %>%
  count(Pathway, sort = FALSE) %>%
  print()

plot_df %>%
  select(Gene, Pathway, SigIn) %>%
  arrange(Pathway, Gene) %>%
  print(n = 200)