# =========================================================
# Acute pathway heatmap from DEPresults_v2 centered intensities
# Includes all user-defined gene lists, grouped by pathway/list
# and clusters similar-behaving proteins WITHIN each pathway
# =========================================================

# DEPresults_v2 <- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")

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

# Acute centered columns to display
acute_centered_cols <- c(
  "RDH12_control_atRAL5hr_centered",
  "RDH12_100_atRAL5hr_centered",
  "RDH12_200_atRAL5hr_centered"
)

# Pretty labels for heatmap columns
acute_col_labels <- c(
  "RDH12_control_atRAL5hr_centered" = "Veh",
  "RDH12_100_atRAL5hr_centered"     = "100uM",
  "RDH12_200_atRAL5hr_centered"     = "200uM"
)

# Output files
out_file_pdf <- "Proteomic_Figs/acute_pathway_heatmap_centered_allLists_clusteredWithinPathway.pdf"
out_file_svg <- "Proteomic_Figs/acute_pathway_heatmap_centered_allLists_clusteredWithinPathway.svg"
out_file_csv <- "Proteomic_Figs/acute_pathway_heatmap_input_table_allLists_clusteredWithinPathway.csv"

# Cluster columns?
cluster_cols <- FALSE

# Remove genes with all NA across displayed columns
remove_all_na_rows <- TRUE

# Distance and linkage for within-pathway clustering
row_distance_method <- "euclidean"
row_clustering_method <- "complete"

# Scale rows before clustering/plotting?
# FALSE = use centered intensities directly
# TRUE  = z-score each row across Veh/100uM/200uM
scale_rows <- FALSE

# ---------------------------------------------------------
# 2. PATHWAY / GENE LISTS
# ---------------------------------------------------------

# https://www.kegg.jp/pathway/map04141

er_genes <- c(
   "EIF2AK3", "ERN1", "DNAJB9",
  "SEC61A1", "SEC61A2", "SEC63", "SEC62", "SEC61B", "SEC61G",
  "OSTC", "DDOST", "RPN1", "RPN2", "DAD1", "STT3A", "STT3B", "MAGT1", "TMEM258", "KRTCAP2", "CKAP4",
  "RRBP1", "HYOU1", "DNAJA2", "DNAJA1", "DNAJB5", "DNAJC10", "DNAJC11", "DNAJC9", "DNAJA3",
  "DNAJB1", "DNAJB11", "DNAJB12", "DNAJB14", "DNAJB2", "DNAJB4", "DNAJB6", "DNAJB9", "DNAJC1",
  "SIL1",  "PRKCSH",
  "GANAB", "MOGS", "CALR", "CANX", "PDIA3", "PDIA4", "PDIA5", "PDIA6",
  "MAN1A1", "MAN1B1", "LMAN1", "LMAN2", "LMAN2L",
  "PREB", "SAR1A", "SAR1B", "SEC13", "SEC31A", "SEC23A", "SEC23B"
)

hsp_genes <- c(
  "HSPA1L", "HSPA12A", "HSPA14", "HSPA13", "HSPA5",
  "HSP90B1", "HSP90AA1", "HSP90AA4P", "HSP90AB1", "HSP90AB4P",
  "HSPA2", "HSPA4", "HSPA4L", "HSPA6", "HSPA8", "HSPA9",
  "HSPB1", "HSPBAP1", "HSPBP1", "HSPD1", "HSPE1", "HSPH1"
)

Ubi_genes <- c(
  "VCPIP1", "ATXN3", "MARCHF5", "MARCHF6", "UBE2G2", "UBE2J1", "UBE2J2"
)

ERAD_genes <- c(
  "EDEM2", "EDEM3", "ERO1A", "ERO1B", "PDIA3", "PDIA4", "PDIA5", "PDIA6",
  "ERP29", "ERP44", "OS9", "ERLEC1", "BCAP29", "BCAP31",
  "UBXN1", "UBXN2A", "UBXN2B", "UBXN4", "UBXN6", "UBXN7", "UBXN8",
  "DERL1", "DERL2", "TRAM1", "BAG4", "BAG2", "BAG1", "SVIP",
  "VCP", "NPLOC4", "UFD1", "PLAA", "UBQLN1", "UBQLN2", "UBQLN4"
)

er_stress_genes <- c(
  "TMEM259", "UGGT1", "UGGT2", "ATF6B", "ATF6", "WFS1", "TRAF2", "XBP1"
)

upr_genes <- c(
  "ATF4", "DDIT3", "PPP1R15A", "XBP1", "HERPUD1", "PDIA4",
  "EIF2AK3", "EIF2S1", "PPP1R15B", "NRF1", "DDI2"
)

jnk1_genes <- c(
  "MAP3K1", "MAP3K20", "MAPK8", "MAP2K4", "MAP3K4", "MAP2K7",
  "MAP3K21", "MAP3K7", "JUN", "ATF2", "FOS", "MAPK8IP2", "SH3RF1", "ITCH",
  "PIAS4",  "TRAF2", "LRR1", "WWOX", "GRIPAP1", "MINK1", "PRKD1", "LYN","GNA12",
  "STK24", "BCL2", "TAB1", "YAP1", "SIRT1", "HSF1", "RAC1","DAG1","DAXX","EGFR","BECN1"
  
)

# ---------------------------------------------------------
# 3. BUILD PATHWAY TABLE
# ---------------------------------------------------------

# A gene can appear in more than one list.
# To keep one row per gene in the heatmap, assign each gene to the first list it appears in.
# Order below controls priority and display grouping.

pathway_tbl <- bind_rows(
  tibble(Gene = er_genes,         Pathway = "ER-associated"),
  tibble(Gene = hsp_genes,        Pathway = "HSP/chaperones"),
  tibble(Gene = Ubi_genes,        Pathway = "DUB"),
  tibble(Gene = ERAD_genes,       Pathway = "ERAD"),
  tibble(Gene = er_stress_genes,  Pathway = "ER stress"),
  tibble(Gene = upr_genes,        Pathway = "UPR"),
  tibble(Gene = jnk1_genes,       Pathway = "JNK1 pathway")
) %>%
  distinct(Gene, .keep_all = TRUE)

# Display order of pathway groups
pathway_order <- c(
  "ER-associated",
  "HSP/chaperones",
  "DUB",
  "ERAD",
  "ER stress",
  "UPR",
  "JNK1 pathway"
)

pathway_tbl <- pathway_tbl %>%
  mutate(Pathway = factor(Pathway, levels = pathway_order))

# ---------------------------------------------------------
# 4. CHECK REQUIRED COLUMNS
# ---------------------------------------------------------

required_cols <- c("Gene", "name", "ID", acute_centered_cols)

missing_cols <- setdiff(required_cols, colnames(DEPresults_v2))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from DEPresults_v2:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 5. SUBSET TO GENES OF INTEREST
# ---------------------------------------------------------

plot_df <- DEPresults_v2 %>%
  mutate(
    Gene = dplyr::coalesce(Gene, name)
  ) %>%
  inner_join(pathway_tbl, by = "Gene") %>%
  select(
    Gene,
    name,
    ID,
    Pathway,
    all_of(acute_centered_cols)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

if (nrow(plot_df) == 0) {
  stop("No genes from pathway_tbl were found in DEPresults_v2.")
}

# ---------------------------------------------------------
# 6. BUILD HEATMAP MATRIX
# ---------------------------------------------------------

heat_mat <- plot_df %>%
  select(Gene, all_of(acute_centered_cols)) %>%
  as.data.frame()

rownames(heat_mat) <- heat_mat$Gene
heat_mat$Gene <- NULL
heat_mat <- as.matrix(heat_mat)

colnames(heat_mat) <- acute_col_labels[colnames(heat_mat)]

# ---------------------------------------------------------
# 7. ROW ANNOTATION
# ---------------------------------------------------------

row_annot <- plot_df %>%
  select(Gene, Pathway) %>%
  distinct() %>%
  as.data.frame()

rownames(row_annot) <- row_annot$Gene
row_annot$Gene <- NULL

# ---------------------------------------------------------
# 8. OPTIONAL: REMOVE ROWS WITH ALL NA
# ---------------------------------------------------------

if (remove_all_na_rows) {
  keep_rows <- rowSums(!is.na(heat_mat)) > 0
  heat_mat <- heat_mat[keep_rows, , drop = FALSE]
  row_annot <- row_annot[keep_rows, , drop = FALSE]
}

# ---------------------------------------------------------
# 9. OPTIONAL: SCALE EACH ROW
# ---------------------------------------------------------

if (scale_rows) {
  heat_mat <- t(scale(t(heat_mat)))
}

# ---------------------------------------------------------
# 10. CLUSTER WITHIN EACH PATHWAY
# ---------------------------------------------------------

cluster_within_group <- function(mat_block,
                                 distance_method = "euclidean",
                                 clustering_method = "complete") {
  if (nrow(mat_block) <= 1) {
    return(rownames(mat_block))
  }
  
  # For clustering only: fill NAs with row mean if needed
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
# 11. GAPS BETWEEN PATHWAY BLOCKS
# ---------------------------------------------------------

pathway_counts <- table(factor(row_annot$Pathway, levels = pathway_order))
pathway_counts <- pathway_counts[pathway_counts > 0]
gaps_row <- cumsum(pathway_counts)
gaps_row <- gaps_row[-length(gaps_row)]

# ---------------------------------------------------------
# 12. COLORS
# ---------------------------------------------------------

pathway_colors <- list(
  Pathway = c(
    "ER-associated"  = "#E69F00",
    "HSP/chaperones" = "#56B4E9",
    "DUB"            = "#009E73",
    "ERAD"           = "#F0E442",
    "ER stress"      = "#D55E00",
    "UPR"            = "#CC79A7",
    "JNK1 pathway"   = "#0072B2"
  )
)

# ---------------------------------------------------------
# 13. DRAW HEATMAP
# ---------------------------------------------------------

draw_heatmap <- function() {
  pheatmap(
    mat = heat_mat,
    annotation_row = row_annot,
    annotation_colors = pathway_colors,
    cluster_rows = FALSE,   # already clustered within pathway manually
    cluster_cols = cluster_cols,
    gaps_row = gaps_row,
    labels_row = display_labels,
    fontsize_row = 8,
    fontsize_col = 12,
    border_color = "grey80",
    angle_col = 0,
    main = paste0(
      "Acute RDH12 pathway proteins (centered intensities)\n",
      "rows clustered within pathway groups"
    )
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 14. SAVE HEATMAP
# ---------------------------------------------------------

pdf(out_file_pdf, width = 7, height = 14)
draw_heatmap()
dev.off()

svg(out_file_svg, width = 7, height = 14)
draw_heatmap()
dev.off()

# ---------------------------------------------------------
# 15. EXPORT INPUT TABLE IN DISPLAY ORDER
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
# 16. QUICK CHECKS
# ---------------------------------------------------------

cat("Number of unique genes requested:", nrow(pathway_tbl), "\n")
cat("Number of genes found in DEPresults_v2:", nrow(plot_df), "\n")
cat("Rows displayed in heatmap:", nrow(heat_mat), "\n")

plot_df %>%
  count(Pathway, sort = FALSE) %>%
  print()