# =========================================================
# Multi-pathway heatmap from DEPresults_v2 centered intensities
#
# Shows ALL 9 groups:
#   Acute RDH12: Veh, 100uM, 200uM
#   Recovery RDH12: Veh, 100uM, 200uM
#   Recovery WT: Veh, 100uM, 200uM
#
# Keeps only proteins significant in at least one acute contrast
# BUT removes significance labels/markings from the plot
#
# Labels pathways directly on the heatmap as row blocks
# Adds column annotations for Phase / Cell Type / Treatment
# WT recovery columns are placed at the far right
# Column labels are rotated 90 degrees at the bottom
# Exports PDF, SVG, and CSV of plotted data
# Drops proteins that failed XIC curation (xic_curated_exclusions.R) and
# marks single-peptide protein groups with "*" in the row label
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("Proteomic_Rscripts/Figures/xic_curated_exclusions.R")
mark_single_peptide <- TRUE

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------

DEPresults_v2 <- read.csv(
  "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dir.create("Proteomic_Figs/enriched", recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.1
lfc_cutoff  <- log2(1.4)

cluster_cols <- FALSE
remove_all_na_rows <- TRUE
row_distance_method <- "euclidean"
row_clustering_method <- "complete"
scale_rows <- FALSE

out_file_pdf <- "Proteomic_Figs/enriched/Fig4a_pathway_heatmap_centered_noSigAnno_WTfarRight.pdf"
out_file_svg <- "Proteomic_Figs/enriched/Fig4a_pathway_heatmap_centered_noSigAnno_WTfarRight.svg"
out_file_csv <- "Proteomic_Figs/enriched/Fig4a_pathway_heatmap_centered_noSigAnno_WTfarRight_input.csv"

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------

pick_first_existing <- function(candidates, all_cols, label_for_error) {
  hit <- candidates[candidates %in% all_cols]
  
  if (length(hit) == 0) {
    stop(
      paste0(
        "Could not find a column for: ", label_for_error, "\n",
        "Tried:\n", paste(candidates, collapse = "\n")
      )
    )
  }
  
  hit[1]
}

cluster_within_group <- function(mat_block,
                                 distance_method = "euclidean",
                                 clustering_method = "complete") {
  if (nrow(mat_block) <= 1) {
    return(rownames(mat_block))
  }
  
  mat_for_clust <- mat_block
  
  for (i in seq_len(nrow(mat_for_clust))) {
    vals <- mat_for_clust[i, ]
    
    if (anyNA(vals)) {
      row_mean <- mean(vals, na.rm = TRUE)
      if (is.nan(row_mean)) row_mean <- 0
      vals[is.na(vals)] <- row_mean
      mat_for_clust[i, ] <- vals
    }
  }
  
  hc <- hclust(
    dist(mat_for_clust, method = distance_method),
    method = clustering_method
  )
  
  rownames(mat_block)[hc$order]
}

# ---------------------------------------------------------
# 3. FIND DISPLAY COLUMNS
# ---------------------------------------------------------

all_cols <- colnames(DEPresults_v2)

acute_veh_col <- pick_first_existing(
  c("RDH12_control_atRAL5hr_centered"),
  all_cols,
  "Acute RDH12 Veh centered"
)

acute_100_col <- pick_first_existing(
  c("RDH12_100_atRAL5hr_centered"),
  all_cols,
  "Acute RDH12 100uM centered"
)

acute_200_col <- pick_first_existing(
  c("RDH12_200_atRAL5hr_centered"),
  all_cols,
  "Acute RDH12 200uM centered"
)

rec_rdh12_veh_col <- pick_first_existing(
  c("RDH12_control_atRAL5hr.24h_recvr_centered"),
  all_cols,
  "Recovery RDH12 Veh centered"
)

rec_rdh12_100_col <- pick_first_existing(
  c("RDH12_100_atRAL5hr.24h_recvr_centered"),
  all_cols,
  "Recovery RDH12 100uM centered"
)

rec_rdh12_200_col <- pick_first_existing(
  c("RDH12_200_atRAL5hr.24h_recvr_centered"),
  all_cols,
  "Recovery RDH12 200uM centered"
)

rec_wt_veh_col <- pick_first_existing(
  c(
    "GFP_control_atRAL5hr.24h_recvr_centered",
    "WT_control_atRAL5hr.24h_recvr_centered",
    "Control_control_atRAL5hr.24h_recvr_centered"
  ),
  all_cols,
  "Recovery WT Veh centered"
)

rec_wt_100_col <- pick_first_existing(
  c(
    "GFP_100_atRAL5hr.24h_recvr_centered",
    "WT_100_atRAL5hr.24h_recvr_centered",
    "Control_100_atRAL5hr.24h_recvr_centered"
  ),
  all_cols,
  "Recovery WT 100uM centered"
)

rec_wt_200_col <- pick_first_existing(
  c(
    "GFP_200_atRAL5hr.24h_recvr_centered",
    "WT_200_atRAL5hr.24h_recvr_centered",
    "Control_200_atRAL5hr.24h_recvr_centered"
  ),
  all_cols,
  "Recovery WT 200uM centered"
)

# Column order:
# Acute RDH12 first, Recovery RDH12 second, Recovery WT far right
display_centered_cols <- c(
  acute_veh_col,
  acute_100_col,
  acute_200_col,
  rec_rdh12_veh_col,
  rec_rdh12_100_col,
  rec_rdh12_200_col,
  rec_wt_veh_col,
  rec_wt_100_col,
  rec_wt_200_col
)

display_col_labels <- setNames(
  c(
    "Acute RDH12 Veh",
    "Acute RDH12 100uM",
    "Acute RDH12 200uM",
    "Recovery RDH12 Veh",
    "Recovery RDH12 100uM",
    "Recovery RDH12 200uM",
    "Recovery WT Veh",
    "Recovery WT 100uM",
    "Recovery WT 200uM"
  ),
  display_centered_cols
)

# ---------------------------------------------------------
# 4. CONTRASTS USED TO KEEP PROTEINS
# ---------------------------------------------------------

acute_contrasts <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr"
)

acute_contrast_labels <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr" = "100uM vs Veh",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr"     = "200uM vs 100uM",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr" = "200uM vs Veh"
)

sig_required_cols <- unlist(lapply(
  acute_contrasts,
  function(x) c(paste0(x, "_p.adj"), paste0(x, "_ratio"))
))

required_cols <- c("Gene", "name", display_centered_cols, sig_required_cols)

missing_cols <- setdiff(required_cols, colnames(DEPresults_v2))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from DEPresults_v2:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 5. PATHWAY / GENE LISTS
# ---------------------------------------------------------

Synthesis_PA <- c(
    
  "AGPAT4", "GPAT4", "AGPAT5",  "DGKE",
  "CHKA"
)

Heme_related <- c(
  "ABCB10", "SLC11A2", "UROS", "IBA57"
)

CCT_TriC_genes <- c(
  "CLIP1", "FLNA"
)

ERK_MAPK <- c(
  "FN1", "WNT5A", "ERBB2", "FGFR2", "SPRY2", "ATF3"
)

KEAP1Nrf2 <- c(
  "HMOX1",
  "GCLC", "GCLM", "GSS",
  "TXNRD1", "GSR",
  "ABCC1",
  "MGST2", "MGST3",
  "GSTO1", "GSTZ1",
  "AKR1A1",
  "BLVRB",
  "NFE2L2", "HMOX1",
  "GCLC", "GCLM",
  "TXNRD1", "ABCC1",
  "BACH1",
  "CREBBP",
  "BACH2",
  "RXRA",
  "SLC39A7",
  "FTH1", "FTL", "TF", "CP",
  "SLC11A2", "SLC40A1", "STEAP3",
  "PCBP1", "PCBP2", "FTMT",
  "GPX4", "ACSL4", "ACSL5", "ACSL6",
  "LPCAT3", "ALOX15",
  "SAT1", "SAT2",
  "VDAC2", "VDAC3",
  "CYBB", "SLC3A2"
  
)

# ---------------------------------------------------------
# 6. BUILD PATHWAY TABLE
# ---------------------------------------------------------

pathway_tbl_all <- bind_rows(
  tibble(Gene = Synthesis_PA,    Pathway = "Lipid synthesis"),
  tibble(Gene = Heme_related,    Pathway = "Heme-related proteins"),
  tibble(Gene = CCT_TriC_genes,  Pathway = "CCT/TRiC"),
  tibble(Gene = ERK_MAPK,        Pathway = "ERK/MAPK"),
  tibble(Gene = KEAP1Nrf2,       Pathway = "Stress response")
)

dup_genes <- pathway_tbl_all %>%
  count(Gene) %>%
  filter(n > 1)

if (nrow(dup_genes) > 0) {
  message("These genes appeared in more than one pathway list and were assigned to the first pathway only:")
  print(dup_genes)
}

pathway_tbl <- pathway_tbl_all %>%
  distinct(Gene, .keep_all = TRUE) %>%
  apply_xic_curation(gene_col = "Gene", label = "Fig 5E")

pathway_order <- c(
  "Lipid synthesis",
  "Heme-related proteins",
  "CCT/TRiC",
  "ERK/MAPK",
  "Stress response"
)

pathway_tbl <- pathway_tbl %>%
  mutate(Pathway = factor(Pathway, levels = pathway_order))

# ---------------------------------------------------------
# 7. FLAG PROTEINS SIGNIFICANT IN >=1 ACUTE CONTRAST
# ---------------------------------------------------------

plot_df <- DEPresults_v2 %>%
  mutate(
    Gene = dplyr::coalesce(na_if(Gene, ""), na_if(name, ""))
  ) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Gene, .keep_all = TRUE)

for (ct in acute_contrasts) {
  padj_col  <- paste0(ct, "_p.adj")
  ratio_col <- paste0(ct, "_ratio")
  out_col   <- acute_contrast_labels[[ct]]
  
  plot_df[[out_col]] <- ifelse(
    !is.na(plot_df[[padj_col]]) &
      !is.na(plot_df[[ratio_col]]) &
      plot_df[[padj_col]] <= padj_cutoff &
      abs(plot_df[[ratio_col]]) >= lfc_cutoff,
    "Yes",
    "No"
  )
}

plot_df <- plot_df %>%
  mutate(
    keep_sig_any = if_else(
      (`100uM vs Veh` == "Yes") |
        (`200uM vs 100uM` == "Yes") |
        (`200uM vs Veh` == "Yes"),
      TRUE,
      FALSE
    )
  )

# ---------------------------------------------------------
# 8. SUBSET TO REQUESTED PATHWAYS + SIG PROTEINS
# ---------------------------------------------------------

plot_df <- pathway_tbl %>%
  left_join(plot_df, by = "Gene") %>%
  filter(keep_sig_any)

if (nrow(plot_df) == 0) {
  stop("No proteins from the pathway lists passed the acute significance filter.")
}

# ---------------------------------------------------------
# 9. BUILD MATRIX + PATHWAY ANNOTATION ONLY
# ---------------------------------------------------------

heat_df <- plot_df %>%
  dplyr::select(
    Gene,
    ID,
    Pathway,
    all_of(display_centered_cols)
  )

if (mark_single_peptide) {
  heat_df <- heat_df %>%
    mutate(Gene = label_single_peptide(Gene, ID))
}

if (remove_all_na_rows) {
  heat_df <- heat_df %>%
    filter(rowSums(!is.na(across(all_of(display_centered_cols)))) > 0)
}

if (nrow(heat_df) == 0) {
  stop("All rows were removed after NA filtering.")
}

heat_mat <- heat_df %>%
  dplyr::select(all_of(display_centered_cols)) %>%
  as.data.frame()

rownames(heat_mat) <- heat_df$Gene
heat_mat <- as.matrix(heat_mat)

colnames(heat_mat) <- unname(display_col_labels[colnames(heat_mat)])

row_annot_df <- heat_df %>%
  dplyr::select(
    Gene,
    Pathway
  ) %>%
  as.data.frame()

rownames(row_annot_df) <- row_annot_df$Gene
row_annot_df$Gene <- NULL

row_annot_df$Pathway <- factor(row_annot_df$Pathway, levels = pathway_order)

annotation_col <- data.frame(
  Phase = c(
    "Acute", "Acute", "Acute",
    "Recovery", "Recovery", "Recovery",
    "Recovery", "Recovery", "Recovery"
  ),
  CellType = c(
    "RDH12", "RDH12", "RDH12",
    "RDH12", "RDH12", "RDH12",
    "WT", "WT", "WT"
  ),
  Treatment = c(
    "Veh", "100uM", "200uM",
    "Veh", "100uM", "200uM",
    "Veh", "100uM", "200uM"
  ),
  row.names = colnames(heat_mat),
  stringsAsFactors = FALSE
)

annotation_col$Phase <- factor(annotation_col$Phase, levels = c("Acute", "Recovery"))
annotation_col$CellType <- factor(annotation_col$CellType, levels = c("RDH12", "WT"))
annotation_col$Treatment <- factor(annotation_col$Treatment, levels = c("Veh", "100uM", "200uM"))

if (scale_rows) {
  heat_mat <- t(scale(t(heat_mat)))
}

# ---------------------------------------------------------
# 10. ORDER ROWS BY PATHWAY, THEN CLUSTER WITHIN PATHWAY
# ---------------------------------------------------------

ordered_genes <- character(0)

for (pw in pathway_order) {
  genes_in_pw <- rownames(row_annot_df)[row_annot_df$Pathway == pw]
  
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
row_annot_df <- row_annot_df[ordered_genes, , drop = FALSE]

pathway_split <- factor(as.character(row_annot_df$Pathway), levels = pathway_order)
pathway_split <- droplevels(pathway_split)

# ---------------------------------------------------------
# 11. COLORS
# ---------------------------------------------------------

pathway_colors <- c(
  "Lipid synthesis" = "#E69F00",
  "Heme-related proteins" = "#56B4E9",
  "CCT/TRiC" = "#F0E442",
  "ERK/MAPK" = "#999999",
  "Stress response" = "#004949"
)

phase_colors <- c(
  "Acute" = "#377EB8",
  "Recovery" = "#FF7F00"
)

celltype_colors <- c(
  "RDH12" = "#1B9E77",
  "WT" = "#999999"
)

treatment_colors <- c(
  "Veh" = "#F0F0F0",
  "100uM" = "#B2DF8A",
  "200uM" = "#FB9A99"
)

# Force 0 to be the center of the color scale
max_abs <- max(abs(heat_mat), na.rm = TRUE)

heat_col_fun <- colorRamp2(
  c(-max_abs, 0, max_abs),
  c("#2400D8", "#FFFFEA", "#A50021")
)

# ---------------------------------------------------------
# 12. ANNOTATIONS
# ---------------------------------------------------------

# Only pathway annotation on the left.
# Significance annotation bars have been removed.
left_anno <- rowAnnotation(
  Pathway = row_annot_df$Pathway,
  col = list(
    Pathway = pathway_colors
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(4, "mm")
)

top_anno <- HeatmapAnnotation(
  Phase = annotation_col$Phase,
  CellType = annotation_col$CellType,
  Treatment = annotation_col$Treatment,
  col = list(
    Phase = phase_colors,
    CellType = celltype_colors,
    Treatment = treatment_colors
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(4, "mm")
)

# ---------------------------------------------------------
# 13. DRAW HEATMAP
# ---------------------------------------------------------

ht <- Heatmap(
  heat_mat,
  name = "Centered\nintensity",
  col = heat_col_fun,
  na_col = "grey95",
  left_annotation = left_anno,
  top_annotation = top_anno,
  cluster_rows = FALSE,
  cluster_columns = cluster_cols,
  row_split = pathway_split,
  cluster_row_slices = FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_gap = unit(3, "mm"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 11),
  row_names_side = "right",
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 90,
  border = TRUE,
  rect_gp = gpar(col = "grey85"),
  column_title = "Pathway-centered proteomic heatmap",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
    at = c(-max_abs, -max_abs / 2, 0, max_abs / 2, max_abs),
    labels = round(c(-max_abs, -max_abs / 2, 0, max_abs / 2, max_abs), 2),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)

# ---------------------------------------------------------
# 14. SAVE PDF
# ---------------------------------------------------------

pdf(out_file_pdf, width = 13, height = 18)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)
dev.off()

# ---------------------------------------------------------
# 15. SAVE SVG
# ---------------------------------------------------------

svg(out_file_svg, width = 13, height = 11)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)
dev.off()

# ---------------------------------------------------------
# 16. EXPORT PLOTTED INPUT TABLE
# ---------------------------------------------------------

out_tbl <- data.frame(
  Gene = rownames(row_annot_df),
  Pathway = as.character(row_annot_df$Pathway),
  heat_mat,
  check.names = FALSE
)

write.csv(out_tbl, out_file_csv, row.names = FALSE)

# ---------------------------------------------------------
# 17. QUICK CHECKS
# ---------------------------------------------------------

cat("Rows displayed:", nrow(heat_mat), "\n")
cat("Columns displayed:", ncol(heat_mat), "\n")

cat("\nDisplayed columns:\n")
print(colnames(heat_mat))

cat("\nProteins per pathway shown:\n")
print(table(pathway_split))

cat("\nSingle-peptide protein groups (marked *):",
    paste(grep("\\*$", rownames(heat_mat), value = TRUE), collapse = ", "), "\n")

cat("\nSaved files:\n")
cat(out_file_pdf, "\n")
cat(out_file_svg, "\n")
cat(out_file_csv, "\n")
