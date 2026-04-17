# =========================================================
# Multi-pathway heatmap from DEPresults_v2 centered intensities
# - Shows ALL 9 groups:
#     Acute RDH12: Veh, 100uM, 200uM
#     Recovery WT/GFP: Veh, 100uM, 200uM
#     Recovery RDH12: Veh, 100uM, 200uM
# - Keeps only proteins significant in at least one acute contrast
# - Labels PATHWAYS directly on the heatmap as row blocks
# - Adds row annotations for acute significance
# - Adds column annotations for Phase / Cell Type / Treatment
# - Exports PDF, SVG, and CSV of plotted data
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

out_file_pdf <- "Proteomic_Figs/enriched/Fig1_pathway_heatmap_centered_allGroups_rowSplit.pdf"
out_file_svg <- "Proteomic_Figs/enriched/Fig1_pathway_heatmap_centered_allGroups_rowSplit.svg"
out_file_csv <- "Proteomic_Figs/enriched/Fig1_pathway_heatmap_centered_allGroups_rowSplit_input.csv"

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
# 3. FIND THE 9 DISPLAY COLUMNS
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

display_centered_cols <- c(
  acute_veh_col,
  acute_100_col,
  acute_200_col,
  rec_wt_veh_col,
  rec_wt_100_col,
  rec_wt_200_col,
  rec_rdh12_veh_col,
  rec_rdh12_100_col,
  rec_rdh12_200_col
)

display_col_labels <- c(
  acute_veh_col     = "Acute\nRDH12 Veh",
  acute_100_col     = "Acute\nRDH12 100uM",
  acute_200_col     = "Acute\nRDH12 200uM",
  rec_wt_veh_col    = "Recovery\nWT Veh",
  rec_wt_100_col    = "Recovery\nWT 100uM",
  rec_wt_200_col    = "Recovery\nWT 200uM",
  rec_rdh12_veh_col = "Recovery\nRDH12 Veh",
  rec_rdh12_100_col = "Recovery\nRDH12 100uM",
  rec_rdh12_200_col = "Recovery\nRDH12 200uM"
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

required_cols <- c("Gene", display_centered_cols, sig_required_cols)
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
stress <- c("MOCS3", "PAN2", "ZDHHC18", "ZDHHC13")

Nglycan_prcrsr_biosynthsis <- c(
  "ALG3","DPAGT1","ALG1","ALG2","ALG11","ALG5","ALG8","ALG9",
  "RFT1","DPM1","STT3A","STT3B","RPN2","RPN1"
)

HSP <- c("HSPA6","HSPA12A")

Synthesis_PA <- c(
  "LPCAT1","PLA2G4A","LPCAT4","AGPAT3","DDHD1","GNPAT","MIGA1",
  "AGPAT4","GPD2","GPAT4","AGPAT5","DGAT1","PTDSS1","DGKE",
  "CHKA","ETNK1"
)

MetalIon_SLCtransporters <- c(
  "SLC39A7","SLC39A1","SLC31A1","SLC9A6","SLC39A10","SLC9A1",
  "SLC39A14","SLC30A5","SLC11A2","SLC7A11"
)

Ca <- c("ORAI1", "ATP2A2", "ATP2B1")

CCT_TriC_genes <- c(
  "TUBB4B","TUBA1C","TCP1","CCT2","TUBB2A","TUBAL3","CCT7",
  "CCT3","CCT8","TUBB4A","TUBB6","CCT6A","CCT4","TUBA4A","CCT5"
)

Prefoldin <- c("PFDN1", "PFDN2", "PFDN4", "PFDN5", "PFDN6")

KEAP1Nrf2 <- c(
  "TXNIP","CYBA","BACH2","ATP7A","TBL1X","NFE2L2","ATOX1","COX19",
  "BACH1","CREBBP","PRDX6","COX7A2","ACY1","COX16","SMARCD3","ABCC1",
  "CYP2S1","CHD9","GSTCD","TXN2","COX20","HMOX2","COX6C","ACOX1",
  "SOD2","PRDX2","COX14","BLVRB","CHAC2","GSTO1","PRDX5","AKR1A1",
  "CCS","CYP20A1","MT-CO2","STAT3","TGS1","CYP51A1","HM13","COX15",
  "HIGD1A","HMOX1","MGST2","KEAP1","MGST3","GCLC","BLVRA","CARM1",
  "GPX4","AKR1B1","GSTZ1","TXNRD1","GSR","GSS","COX5B","COX6B1",
  "TXNRD2","GPX1","PRDX3","MAFK","RXRA","NCOR1","NCOR2","COX17",
  "NCOA6","FTH1","FTL","TFRC"
)

# ---------------------------------------------------------
# 6. BUILD PATHWAY TABLE
# ---------------------------------------------------------
pathway_tbl <- bind_rows(
  tibble(Gene = stress,                   Pathway = "Stress"),
  tibble(Gene = Nglycan_prcrsr_biosynthsis, Pathway = "N-glycan precursor biosynthesis"),
  tibble(Gene = HSP,                      Pathway = "HSP"),
  tibble(Gene = Synthesis_PA,             Pathway = "PA / phospholipid synthesis"),
  tibble(Gene = MetalIon_SLCtransporters, Pathway = "Metal-ion / SLC transporters"),
  tibble(Gene = Ca,                       Pathway = "Calcium handling"),
  tibble(Gene = CCT_TriC_genes,           Pathway = "CCT/TRiC"),
  tibble(Gene = Prefoldin,                Pathway = "Prefoldin"),
  tibble(Gene = KEAP1Nrf2,                Pathway = "KEAP1/NRF2")
) %>%
  distinct(Gene, .keep_all = TRUE)

pathway_order <- c(
  "Stress",
  "N-glycan precursor biosynthesis",
  "HSP",
  "PA / phospholipid synthesis",
  "Metal-ion / SLC transporters",
  "Calcium handling",
  "CCT/TRiC",
  "Prefoldin",
  "KEAP1/NRF2"
)

pathway_tbl <- pathway_tbl %>%
  mutate(Pathway = factor(Pathway, levels = pathway_order))

dup_genes <- bind_rows(
  tibble(Gene = stress,                   Pathway = "Stress"),
  tibble(Gene = Nglycan_prcrsr_biosynthsis, Pathway = "N-glycan precursor biosynthesis"),
  tibble(Gene = HSP,                      Pathway = "HSP"),
  tibble(Gene = Synthesis_PA,             Pathway = "PA / phospholipid synthesis"),
  tibble(Gene = MetalIon_SLCtransporters, Pathway = "Metal-ion / SLC transporters"),
  tibble(Gene = Ca,                       Pathway = "Calcium handling"),
  tibble(Gene = CCT_TriC_genes,           Pathway = "CCT/TRiC"),
  tibble(Gene = Prefoldin,                Pathway = "Prefoldin"),
  tibble(Gene = KEAP1Nrf2,                Pathway = "KEAP1/NRF2")
) %>%
  count(Gene) %>%
  filter(n > 1)

if (nrow(dup_genes) > 0) {
  message("These genes appeared in more than one pathway list and were assigned to the first pathway only:")
  print(dup_genes)
}

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
    "Yes", "No"
  )
}

plot_df <- plot_df %>%
  mutate(
    keep_sig_any = if_else(
      (`100uM vs Veh` == "Yes") |
        (`200uM vs 100uM` == "Yes") |
        (`200uM vs Veh` == "Yes"),
      TRUE, FALSE
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
# 9. BUILD MATRIX + ANNOTATIONS
# ---------------------------------------------------------
heat_df <- plot_df %>%
  select(
    Gene, Pathway,
    all_of(names(acute_contrast_labels)),
    all_of(display_centered_cols)
  )

if (remove_all_na_rows) {
  heat_df <- heat_df %>%
    filter(rowSums(!is.na(across(all_of(display_centered_cols)))) > 0)
}

if (nrow(heat_df) == 0) {
  stop("All rows were removed after NA filtering.")
}

heat_mat <- heat_df %>%
  select(all_of(display_centered_cols)) %>%
  as.matrix()

rownames(heat_mat) <- heat_df$Gene
colnames(heat_mat) <- unname(display_col_labels[colnames(heat_mat)])

row_annot_df <- heat_df %>%
  select(Gene, Pathway, all_of(names(acute_contrast_labels))) %>%
  as.data.frame()

rownames(row_annot_df) <- row_annot_df$Gene
row_annot_df$Gene <- NULL

row_annot_df$Pathway <- factor(row_annot_df$Pathway, levels = pathway_order)
row_annot_df$`100uM vs Veh` <- factor(row_annot_df$`100uM vs Veh`, levels = c("Yes", "No"))
row_annot_df$`200uM vs 100uM` <- factor(row_annot_df$`200uM vs 100uM`, levels = c("Yes", "No"))
row_annot_df$`200uM vs Veh` <- factor(row_annot_df$`200uM vs Veh`, levels = c("Yes", "No"))

annotation_col <- data.frame(
  Phase = c(
    "Acute", "Acute", "Acute",
    "Recovery", "Recovery", "Recovery",
    "Recovery", "Recovery", "Recovery"
  ),
  CellType = c(
    "RDH12", "RDH12", "RDH12",
    "WT", "WT", "WT",
    "RDH12", "RDH12", "RDH12"
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
annotation_col$CellType <- factor(annotation_col$CellType, levels = c("WT", "RDH12"))
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
  "Stress" = "#E69F00",
  "N-glycan precursor biosynthesis" = "#56B4E9",
  "HSP" = "#CC79A7",
  "PA / phospholipid synthesis" = "#009E73",
  "Metal-ion / SLC transporters" = "#0072B2",
  "Calcium handling" = "#D55E00",
  "CCT/TRiC" = "#F0E442",
  "Prefoldin" = "#999999",
  "KEAP1/NRF2" = "#004949"
)

sig_colors <- c("Yes" = "#D73027", "No" = "grey90")

phase_colors <- c("Acute" = "#377EB8", "Recovery" = "#FF7F00")
celltype_colors <- c("WT" = "#999999", "RDH12" = "#1B9E77")
treatment_colors <- c("Veh" = "#F0F0F0", "100uM" = "#B2DF8A", "200uM" = "#FB9A99")

heat_col_fun <- colorRamp2(
  c(-2.5, 0, 2.5),
  c("#3B4CC0", "#F7F7F7", "#B40426")
)

# ---------------------------------------------------------
# 12. ANNOTATIONS
# ---------------------------------------------------------
left_anno <- rowAnnotation(
  Pathway = row_annot_df$Pathway,
  `100uM vs Veh` = row_annot_df$`100uM vs Veh`,
  `200uM vs 100uM` = row_annot_df$`200uM vs 100uM`,
  `200uM vs Veh` = row_annot_df$`200uM vs Veh`,
  col = list(
    Pathway = pathway_colors,
    `100uM vs Veh` = sig_colors,
    `200uM vs 100uM` = c("Yes" = "#984EA3", "No" = "grey90"),
    `200uM vs Veh` = c("Yes" = "#4DAF4A", "No" = "grey90")
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
  row_split = pathway_split,              # THIS prints pathway blocks/titles
  cluster_row_slices = FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_gap = unit(3, "mm"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 11),
  row_names_side = "right",
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 0,
  border = TRUE,
  rect_gp = gpar(col = "grey85"),
  column_title = "Pathway-centered proteomic heatmap",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
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
pdf(out_file_pdf, width = 14, height = 10)
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
svg(out_file_svg, width = 14, height = 10)
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
  Pathway = row_annot_df$Pathway,
  `100uM vs Veh` = row_annot_df$`100uM vs Veh`,
  `200uM vs 100uM` = row_annot_df$`200uM vs 100uM`,
  `200uM vs Veh` = row_annot_df$`200uM vs Veh`,
  heat_mat,
  check.names = FALSE
)

write.csv(out_tbl, out_file_csv, row.names = FALSE)

# ---------------------------------------------------------
# 17. QUICK CHECKS
# ---------------------------------------------------------
cat("Rows displayed:", nrow(heat_mat), "\n")
cat("Columns displayed:", ncol(heat_mat), "\n")
cat("Displayed columns:\n")
print(colnames(heat_mat))

cat("\nProteins per pathway shown:\n")
print(table(pathway_split))

cat("\nSaved files:\n")
cat(out_file_pdf, "\n")
cat(out_file_svg, "\n")
cat(out_file_csv, "\n")