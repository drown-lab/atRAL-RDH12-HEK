# =========================================================
# Cluster heatmap of proteins using centered values
# - Uses all 9 *_centered columns
# - Keeps proteins significant in at least one contrast
#   with adjusted p <= 0.05 and abs(FC) >= 2
# - Clusters proteins by similar behavior across conditions
# - Adds column annotations for phase / cell type / treatment
# - Saves PDF, SVG, and CSV of plotted matrix
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(stringr)
})

# ---------------------------------------------------------
# 1. INPUT
# ---------------------------------------------------------
DEPresults_v2 <- read.csv(
  "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(DEPresults_v2)

# =========================================================
# 2. USER SETTINGS
# =========================================================
padj_cutoff <- 0.05
fc_cutoff   <- 2
lfc_cutoff  <- log2(fc_cutoff)

# Heatmap behavior
remove_all_na_rows <- TRUE
cluster_cols <- FALSE
scale_rows <- FALSE   # leave FALSE because values are already centered
show_rownames <- FALSE
fontsize_row <- 6
fontsize_col <- 10

# If you want only the top N most variable proteins in the heatmap, set a number.
# Example: top_n_variable <- 300
# To keep all significant proteins, leave as NULL
top_n_variable <- NULL

# Output
out_dir <- "Proteomic_Figs/ClusterHeatmaps/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "ClusterHeatmap_centered_sigAnyContrast.pdf")
out_svg <- file.path(out_dir, "ClusterHeatmap_centered_sigAnyContrast.svg")
out_csv <- file.path(out_dir, "ClusterHeatmap_centered_sigAnyContrast_matrix.csv")

# ---------------------------------------------------------
# 3. DEFINE CENTERED COLUMNS
# ---------------------------------------------------------
display_centered_cols <- c(
  "GFP_control_atRAL5hr.24h_recvr_centered",
  "GFP_100_atRAL5hr.24h_recvr_centered",
  "GFP_200_atRAL5hr.24h_recvr_centered",
  "RDH12_control_atRAL5hr_centered",
  "RDH12_100_atRAL5hr_centered",
  "RDH12_200_atRAL5hr_centered",
  "RDH12_control_atRAL5hr.24h_recvr_centered",
  "RDH12_100_atRAL5hr.24h_recvr_centered",
  "RDH12_200_atRAL5hr.24h_recvr_centered"
)

missing_centered <- setdiff(display_centered_cols, colnames(DEPresults_v2))
if (length(missing_centered) > 0) {
  stop(
    "These centered columns are missing:\n",
    paste(missing_centered, collapse = "\n")
  )
}

pretty_col_labels <- c(
  "GFP_control_atRAL5hr.24h_recvr_centered"   = "WT Rec Veh",
  "GFP_100_atRAL5hr.24h_recvr_centered"       = "WT Rec 100uM",
  "GFP_200_atRAL5hr.24h_recvr_centered"       = "WT Rec 200uM",
  "RDH12_control_atRAL5hr_centered"           = "RDH12 Acute Veh",
  "RDH12_100_atRAL5hr_centered"               = "RDH12 Acute 100uM",
  "RDH12_200_atRAL5hr_centered"               = "RDH12 Acute 200uM",
  "RDH12_control_atRAL5hr.24h_recvr_centered" = "RDH12 Rec Veh",
  "RDH12_100_atRAL5hr.24h_recvr_centered"     = "RDH12 Rec 100uM",
  "RDH12_200_atRAL5hr.24h_recvr_centered"     = "RDH12 Rec 200uM"
)

# ---------------------------------------------------------
# 4. DETECT ALL CONTRASTS
# ---------------------------------------------------------
all_contrasts <- colnames(DEPresults_v2) %>%
  .[str_detect(., "_ratio$")] %>%
  str_remove("_ratio$")

if (length(all_contrasts) == 0) {
  stop("No contrast columns ending in '_ratio' were found.")
}

# confirm corresponding p.adj columns exist
valid_contrasts <- all_contrasts[
  paste0(all_contrasts, "_p.adj") %in% colnames(DEPresults_v2)
]

if (length(valid_contrasts) == 0) {
  stop("No valid contrasts with both _ratio and _p.adj were found.")
}

# ---------------------------------------------------------
# 5. BUILD SIGNIFICANCE FLAGS
# ---------------------------------------------------------
plot_df <- DEPresults_v2 %>%
  mutate(
    Gene_use = dplyr::coalesce(na_if(Gene, ""), na_if(name, ""))
  ) %>%
  filter(!is.na(Gene_use), Gene_use != "") %>%
  distinct(Gene_use, .keep_all = TRUE)

for (ct in valid_contrasts) {
  ratio_col <- paste0(ct, "_ratio")
  padj_col  <- paste0(ct, "_p.adj")
  flag_col  <- paste0(ct, "_sigflag")
  
  plot_df[[flag_col]] <- !is.na(plot_df[[ratio_col]]) &
    !is.na(plot_df[[padj_col]]) &
    plot_df[[padj_col]] <= padj_cutoff &
    abs(plot_df[[ratio_col]]) >= lfc_cutoff
}

sigflag_cols <- paste0(valid_contrasts, "_sigflag")

plot_df <- plot_df %>%
  mutate(
    keep_sig_any = if_any(all_of(sigflag_cols), identity)
  )

# ---------------------------------------------------------
# 6. SUBSET TO SIGNIFICANT PROTEINS
# ---------------------------------------------------------
heat_df <- plot_df %>%
  filter(keep_sig_any) %>%
  select(Gene = Gene_use, all_of(display_centered_cols))

if (nrow(heat_df) == 0) {
  stop("No proteins passed the significance filter.")
}

# ---------------------------------------------------------
# 7. BUILD MATRIX
# ---------------------------------------------------------
heat_mat <- heat_df %>%
  as.data.frame()

rownames(heat_mat) <- heat_mat$Gene
heat_mat$Gene <- NULL
heat_mat <- as.matrix(heat_mat)

if (remove_all_na_rows) {
  keep_rows <- rowSums(!is.na(heat_mat)) > 0
  heat_mat <- heat_mat[keep_rows, , drop = FALSE]
}

if (nrow(heat_mat) == 0) {
  stop("All rows were removed after NA filtering.")
}

# Optional: keep top N most variable rows
if (!is.null(top_n_variable)) {
  row_var <- apply(heat_mat, 1, function(x) var(x, na.rm = TRUE))
  row_var[is.na(row_var)] <- -Inf
  keep_top <- names(sort(row_var, decreasing = TRUE))[seq_len(min(top_n_variable, nrow(heat_mat)))]
  heat_mat <- heat_mat[keep_top, , drop = FALSE]
}

# Since these are centered values already, do not re-center unless requested
if (scale_rows) {
  heat_mat <- t(scale(t(heat_mat)))
}

# Replace column names with pretty labels
colnames(heat_mat) <- unname(pretty_col_labels[colnames(heat_mat)])

# ---------------------------------------------------------
# 8. COLUMN ANNOTATION
# ---------------------------------------------------------
annotation_col <- data.frame(
  Phase = c(
    "Recovery", "Recovery", "Recovery",
    "Acute", "Acute", "Acute",
    "Recovery", "Recovery", "Recovery"
  ),
  CellType = c(
    "WT", "WT", "WT",
    "RDH12", "RDH12", "RDH12",
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

annotation_colors <- list(
  Phase = c("Acute" = "#377EB8", "Recovery" = "#FF7F00"),
  CellType = c("WT" = "#999999", "RDH12" = "#1B9E77"),
  Treatment = c("Veh" = "#F0F0F0", "100uM" = "#B2DF8A", "200uM" = "#FB9A99")
)

# ---------------------------------------------------------
# 9. OPTIONAL: HANDLE MISSING VALUES FOR CLUSTERING
# ---------------------------------------------------------
# pheatmap cannot cluster properly with many NAs
# replace NA by row mean for clustering/plotting
heat_mat_plot <- heat_mat

for (i in seq_len(nrow(heat_mat_plot))) {
  vals <- heat_mat_plot[i, ]
  if (anyNA(vals)) {
    row_mean <- mean(vals, na.rm = TRUE)
    if (is.nan(row_mean)) row_mean <- 0
    vals[is.na(vals)] <- row_mean
    heat_mat_plot[i, ] <- vals
  }
}

# ---------------------------------------------------------
# 10. DRAW HEATMAP
# ---------------------------------------------------------
draw_heatmap <- function() {
  pheatmap(
    mat = heat_mat_plot,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = TRUE,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    border_color = NA,
    main = paste0(
      "Clustered heatmap of significant proteins\n",
      "Using centered values across acute and recovery conditions"
    )
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 11. SAVE OUTPUTS
# ---------------------------------------------------------
pdf(out_pdf, width = 10, height = 12)
draw_heatmap()
dev.off()

svg(out_svg, width = 10, height = 12)
draw_heatmap()
dev.off()

write.csv(
  data.frame(Gene = rownames(heat_mat_plot), heat_mat_plot, check.names = FALSE),
  out_csv,
  row.names = FALSE
)

# ---------------------------------------------------------
# 12. QUICK CHECKS
# ---------------------------------------------------------
cat("Number of valid contrasts:", length(valid_contrasts), "\n")
cat("Proteins passing significance in at least one contrast:", nrow(heat_mat_plot), "\n")
cat("Centered columns used:\n")
print(display_centered_cols)
cat("\nSaved files:\n")
cat(out_pdf, "\n")
cat(out_svg, "\n")
cat(out_csv, "\n")