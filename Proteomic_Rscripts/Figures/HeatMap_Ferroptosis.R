# =========================================================
# Acute pathway heatmap from DEPresults_v2 centered intensities
# Includes all user-defined gene lists, grouped by pathway/list
# and clusters similar-behaving proteins WITHIN each pathway
# ONLY keeps proteins significant in at least one acute contrast
# ADDS ROW ANNOTATIONS FOR WHICH CONTRAST(S) EACH PROTEIN IS SIGNIFICANT IN
# =========================================================

DEPresults_v2 <- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")

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
padj_cutoff <- 0.05
lfc_cutoff  <- log2(1.5)

# Pretty labels for heatmap columns
acute_col_labels <- c(
  "RDH12_control_atRAL5hr_centered" = "Veh",
  "RDH12_100_atRAL5hr_centered"     = "100uM",
  "RDH12_200_atRAL5hr_centered"     = "200uM"
)

# Output files
out_file_pdf <- "Proteomic_Figs/enriched/acute_pathway_heatmap_centered_sigOnly_Lists_clusteredWithinPathway_wContrastAnnoferrp.pdf"
out_file_svg <- "Proteomic_Figs/enriched/acute_pathway_heatmap_centered_sigOnly_Lists_clusteredWithinPathway_wContrastAnnoFerrop.svg"
out_file_csv <- "Proteomic_Figs/enriched/acute_pathway_heatmap_input_table_sigOnly_clusteredWithinPathway_wContrastAnnoferro.csv"

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

#pos_regCytokine_prdtn_inflamresponse<-c("MAPK9","IL17RA","STAT3","APPL1","KPNA6","TBK1")

ferroptosis <- c("SLC39A14","SLC7A11","ACSL1","TFRC","TF","SLC11A2", "GCLC","HMOX1","ACSL3","ATG7",
                 "SCL3A2", "TP53", "GCLM","GSS", "GPX4","CP", "SLC40A1", "SLC39A8", "PRNP",
                 "PCBP2", "PCBP1", "FTH1" , "FTL", "STEAP3", "VDAC2", "VDAC3", "CYBB",
                 "FTMT", "ACSL4", "ACSL5", "ACSL6","LPCAT3", "SAT1", "SAT2", "ALOX15"
                 
                 )
# ---------------------------------------------------------
# 3. BUILD PATHWAY TABLE
# ---------------------------------------------------------

pathway_tbl <- bind_rows(

  tibble(Gene = ferroptosis,       Pathway = "ferroptosis"
         
  )
) %>%
  distinct(Gene, .keep_all = TRUE)

pathway_order <- c(
  
  "ferroptosis"
)

pathway_tbl <- pathway_tbl %>%
  mutate(Pathway = factor(Pathway, levels = pathway_order))

# ---------------------------------------------------------
# 4. CHECK REQUIRED COLUMNS
# ---------------------------------------------------------

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
    all_of(acute_centered_cols),
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
  select(Gene, all_of(acute_centered_cols)) %>%
  as.data.frame()

rownames(heat_mat) <- heat_mat$Gene
heat_mat$Gene <- NULL
heat_mat <- as.matrix(heat_mat)

colnames(heat_mat) <- acute_col_labels[colnames(heat_mat)]

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
    
    "ferroptosis"   = "#004949"
    
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
      "Acute RDH12 pathway proteins (sig in >=1 acute contrast)\n",
      "Centered intensities; rows clustered within pathway groups"
    )
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 15. SAVE HEATMAP
# ---------------------------------------------------------

pdf(out_file_pdf, width = 8.5, height = 8)
draw_heatmap()
dev.off()

svg(out_file_svg, width = 8.5, height = 8)
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

plot_df %>%
  count(Pathway, sort = FALSE) %>%
  print()

# Optional: inspect which genes were significant in which contrasts
plot_df %>%
  select(Gene, Pathway, SigIn) %>%
  arrange(Pathway, Gene) %>%
  print(n = 200)