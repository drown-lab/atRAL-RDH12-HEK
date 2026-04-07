# =========================================================
# Acute transporter heatmap from DEPresults_v2 centered intensities
# - Uses transporter genes from Transporterproteins_signficant.csv
# - Keeps genes significant in >=1 acute contrast
# - Shows centered values for Veh / 100uM / 200uM
# - Adds row annotations for which contrast(s) each transporter
#   was significant in
# =========================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(forcats)
library(pheatmap)

# ---------------------------------------------------------
# 1. INPUT FILES
# ---------------------------------------------------------

DEPresults_v2 <- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")
transporter_df <- read.csv("Proteomic_output_txts/Transporterproteins_signficant.csv")

# ---------------------------------------------------------
# 2. USER INPUTS
# ---------------------------------------------------------

acute_centered_cols <- c(
  "RDH12_control_atRAL5hr_centered",
  "RDH12_100_atRAL5hr_centered",
  "RDH12_200_atRAL5hr_centered"
)

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

acute_col_labels <- c(
  "RDH12_control_atRAL5hr_centered" = "Veh",
  "RDH12_100_atRAL5hr_centered"     = "100uM",
  "RDH12_200_atRAL5hr_centered"     = "200uM"
)

padj_cutoff <- 0.05
lfc_cutoff  <- log2(1.5)

cluster_cols <- FALSE
remove_all_na_rows <- TRUE
row_distance_method <- "euclidean"
row_clustering_method <- "complete"
scale_rows <- FALSE

out_file_pdf <- "Proteomic_Figs/acute_transporters_heatmap_centered_sigOnly.pdf"
out_file_svg <- "Proteomic_Figs/acute_transporters_heatmap_centered_sigOnly.svg"
out_file_csv <- "Proteomic_Figs/acute_transporters_heatmap_input_table_sigOnly.csv"

# ---------------------------------------------------------
# 3. BUILD TRANSPORTER GENE TABLE
# ---------------------------------------------------------

# Uses Gene column if present, otherwise falls back to name
if (!"Gene" %in% colnames(transporter_df)) {
  stop("Transporter file must contain a 'Gene' column.")
}

transporter_tbl <- transporter_df %>%
  transmute(
    Gene = Gene,
    Group = case_when(
      str_detect(Gene, "^SLC") ~ "SLC transporters",
      str_detect(Gene, "^ABC") ~ "ABC transporters",
      str_detect(Gene, "^ATP") ~ "ATPase-related",
      TRUE ~ "Other transporters"
    )
  ) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Gene, .keep_all = TRUE)

group_order <- c(
  "SLC transporters",
  "ABC transporters",
  "ATPase-related",
  "Other transporters"
)

transporter_tbl <- transporter_tbl %>%
  mutate(Group = factor(Group, levels = group_order))

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
        if (length(hit_names) == 0) "None" else paste(hit_names, collapse = "; ")
      }
    )
  )

# ---------------------------------------------------------
# 6. SUBSET TO TRANSPORTERS AND SIG-ONLY
# ---------------------------------------------------------

plot_df <- DEPresults_sig %>%
  filter(keep_sig_any) %>%
  inner_join(transporter_tbl, by = "Gene") %>%
  select(
    Gene,
    name,
    ID,
    Group,
    keep_sig_any,
    SigIn,
    all_of(sigflag_cols),
    all_of(acute_centered_cols),
    all_of(sig_required_cols)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

if (nrow(plot_df) == 0) {
  stop("No significant transporter genes were found in DEPresults_v2.")
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
    Group = Group,
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
# 11. CLUSTER WITHIN EACH GROUP
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

for (grp in group_order) {
  genes_in_grp <- rownames(row_annot)[row_annot$Group == grp]
  if (length(genes_in_grp) == 0) next
  
  mat_block <- heat_mat[genes_in_grp, , drop = FALSE]
  
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
# 12. GAPS BETWEEN GROUP BLOCKS
# ---------------------------------------------------------

group_counts <- table(factor(row_annot$Group, levels = group_order))
group_counts <- group_counts[group_counts > 0]
gaps_row <- cumsum(group_counts)
gaps_row <- gaps_row[-length(gaps_row)]

# ---------------------------------------------------------
# 13. COLORS
# ---------------------------------------------------------

annotation_colors <- list(
  Group = c(
    "SLC transporters" = "#009E73",
    "ABC transporters" = "#56B4E9",
    "ATPase-related"   = "#CC79A7",
    "Other transporters" = "#E69F00"
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
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    gaps_row = gaps_row,
    labels_row = display_labels,
    fontsize_row = 8,
    fontsize_col = 12,
    border_color = "grey80",
    angle_col = 0,
    main = paste0(
      "Acute treatment transporter proteins\n",
      "Centered intensities; rows clustered within transporter groups"
    )
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 15. SAVE HEATMAP
# ---------------------------------------------------------

pdf(out_file_pdf, width = 8.5, height = 12)
draw_heatmap()
dev.off()

svg(out_file_svg, width = 8.5, height = 12)
draw_heatmap()
dev.off()

# ---------------------------------------------------------
# 16. EXPORT INPUT TABLE IN DISPLAY ORDER
# ---------------------------------------------------------

plot_df_export <- plot_df %>%
  mutate(Gene = factor(Gene, levels = ordered_genes)) %>%
  arrange(Gene)

write.csv(
  plot_df_export,
  file = out_file_csv,
  row.names = FALSE
)

# ---------------------------------------------------------
# 17. QUICK CHECKS
# ---------------------------------------------------------

cat("Number of transporter genes in transporter file:", nrow(transporter_tbl), "\n")
cat("Number of significant transporter genes in >=1 acute contrast:", nrow(plot_df), "\n")
cat("Rows displayed in heatmap:", nrow(heat_mat), "\n")

plot_df %>%
  count(Group, sort = FALSE) %>%
  print()

plot_df %>%
  select(Gene, Group, SigIn) %>%
  arrange(Group, Gene) %>%
  print(n = 200)