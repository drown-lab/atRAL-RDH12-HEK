# =========================================================
# Recovery targeted lipid heatmap for main-text figure
# Uses: Recovery_DEA_results_v1_corrected.csv
# Goal: compact heatmap of biologically relevant recovery lipids
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(pheatmap)
  library(tidyr)
})

# ---------------------------------------------------------
# 1. INPUT
# ---------------------------------------------------------
recovery_file <- "Lipidomics/output_txts/Recovery_DEA_results_v1_corrected.csv"

out_dir <- "Lipidomics/Figures/Recovery/Recovery_targeted_heatmap/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Recovery_targeted_lipid_heatmap.pdf")
out_png <- file.path(out_dir, "Recovery_targeted_lipid_heatmap.png")
out_csv <- file.path(out_dir, "Recovery_targeted_lipid_heatmap_table.csv")

# ---------------------------------------------------------
# 2. READ DATA
# ---------------------------------------------------------
lipid_df <- read.csv(
  recovery_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nColumn names in recovery file:\n")
print(colnames(lipid_df))

# ---------------------------------------------------------
# 3. HELPER TO FIND COLUMN NAMES ROBUSTLY
# ---------------------------------------------------------
find_first_col <- function(candidates, df_names, required = TRUE, label = NULL) {
  hit <- candidates[candidates %in% df_names][1]
  if (length(hit) == 0 || is.na(hit)) {
    if (required) {
      stop(
        paste0(
          "Could not find column for: ",
          ifelse(is.null(label), "unknown", label),
          "\nTried:\n",
          paste(candidates, collapse = "\n")
        )
      )
    } else {
      return(NA_character_)
    }
  }
  return(hit)
}

# ---------------------------------------------------------
# 4. DETECT CENTERED INTENSITY COLUMNS
# ---------------------------------------------------------
df_names <- colnames(lipid_df)

name_col <- find_first_col(
  candidates = c("name", "Name", "Lipid", "lipid"),
  df_names = df_names,
  label = "lipid name"
)

rdh12_veh_col <- find_first_col(
  c("RDH12_Veh_centered"),
  df_names, label = "RDH12 Veh centered"
)

rdh12_100_col <- find_first_col(
  c("RDH12_100_centered", "RDH12_100uM_centered"),
  df_names, label = "RDH12 100 centered"
)

rdh12_200_col <- find_first_col(
  c("RDH12_200_centered", "RDH12_200uM_centered"),
  df_names, label = "RDH12 200 centered"
)

wt_veh_col <- find_first_col(
  c("Ctrl_Veh_centered", "WT_Veh_centered", "GFP_Veh_centered", "Control_Veh_centered"),
  df_names, label = "WT/Ctrl Veh centered"
)

wt_100_col <- find_first_col(
  c("Ctrl_100_centered", "WT_100_centered", "GFP_100_centered", "Control_100_centered"),
  df_names, label = "WT/Ctrl 100 centered"
)

wt_200_col <- find_first_col(
  c("Ctrl_200_centered", "WT_200_centered", "GFP_200_centered", "Control_200_centered"),
  df_names, label = "WT/Ctrl 200 centered"
)

# ---------------------------------------------------------
# 5. DETECT RELEVANT STATISTICS COLUMNS FOR RANKING
#    (we use these only to help choose representative lipids)
# ---------------------------------------------------------
candidate_padj_cols <- c(
  # RDH12 treatment effects
  "RDH12_100_vs_RDH12_Veh_p.adj",
  "RDH12_200_vs_RDH12_Veh_p.adj",
  
  # WT/Ctrl treatment effects
  "Ctrl_100_vs_Ctrl_Veh_p.adj",
  "Ctrl_200_vs_Ctrl_Veh_p.adj",
  "WT_100_vs_WT_Veh_p.adj",
  "WT_200_vs_WT_Veh_p.adj",
  
  # RDH12 vs WT/Ctrl within each recovery condition
  "RDH12_Veh_vs_Ctrl_Veh_p.adj",
  "RDH12_100_vs_Ctrl_100_p.adj",
  "RDH12_200_vs_Ctrl_200_p.adj",
  "RDH12_Veh_vs_WT_Veh_p.adj",
  "RDH12_100_vs_WT_100_p.adj",
  "RDH12_200_vs_WT_200_p.adj"
)

candidate_ratio_cols <- c(
  "RDH12_100_vs_RDH12_Veh_ratio",
  "RDH12_200_vs_RDH12_Veh_ratio",
  
  "Ctrl_100_vs_Ctrl_Veh_ratio",
  "Ctrl_200_vs_Ctrl_Veh_ratio",
  "WT_100_vs_WT_Veh_ratio",
  "WT_200_vs_WT_Veh_ratio",
  
  "RDH12_Veh_vs_Ctrl_Veh_ratio",
  "RDH12_100_vs_Ctrl_100_ratio",
  "RDH12_200_vs_Ctrl_200_ratio",
  "RDH12_Veh_vs_WT_Veh_ratio",
  "RDH12_100_vs_WT_100_ratio",
  "RDH12_200_vs_WT_200_ratio"
)

padj_cols_found <- candidate_padj_cols[candidate_padj_cols %in% df_names]
ratio_cols_found <- candidate_ratio_cols[candidate_ratio_cols %in% df_names]

if (length(padj_cols_found) == 0) {
  stop("No expected recovery adjusted p-value columns were found.")
}
if (length(ratio_cols_found) == 0) {
  stop("No expected recovery ratio columns were found.")
}

# ---------------------------------------------------------
# 6. USER-DEFINED TARGETS
# ---------------------------------------------------------
# Exact lipids explicitly discussed in your text.
# Edit these after first run if names differ slightly in your file.

exact_targets <- c(
  "PC 34:1",
  "PC 36:7",
  "DG 30:3",
  "DG 32:5",
  "DG 34:2",
  "PS 34:4",
  "PC(18:0/Aze)",
  "PC (18:0/Aze)",
  "PC(16:0/8:(COOH))",
  "PC (16:0/8:(COOH))",
  "PC(Azelaoyl-PAF)"
)

# Pattern-based additions for representative species if exact names differ
pattern_targets <- c(
  "^PC O-",      # ether-PC / plasmalogen-like
  "^PC P-",      # plasmalogen naming
  "^PC\\s",      # additional PCs
  "^PE\\s",      # PE species
  "^PI\\s",      # PI species
  "^DG\\s",      # DG species
  "^PS\\s",      # PS species
  "^TG\\s",      # optional TG representatives
  "^CE\\s"       # optional CE representatives
)

# Number to keep from each group, ranked by strongest significance
n_keep_per_pattern <- list(
  "^PC O-" = 4,
  "^PC P-" = 2,
  "^PC\\s" = 4,
  "^PE\\s" = 3,
  "^PI\\s" = 3,
  "^DG\\s" = 3,
  "^PS\\s" = 2,
  "^TG\\s" = 1,
  "^CE\\s" = 1
)

# ---------------------------------------------------------
# 7. PREP AND RANK LIPIDS
# ---------------------------------------------------------
lipid_df2 <- lipid_df %>%
  mutate(
    name = str_trim(.data[[name_col]])
  )

# compute minimum adjusted p-value across available relevant contrasts
padj_matrix <- as.matrix(lipid_df2[, padj_cols_found, drop = FALSE])
ratio_matrix <- as.matrix(lipid_df2[, ratio_cols_found, drop = FALSE])

lipid_df2$rank_padj <- apply(padj_matrix, 1, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
})

lipid_df2$rank_abs_ratio <- apply(ratio_matrix, 1, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else max(abs(x), na.rm = TRUE)
})

# ---------------------------------------------------------
# 8. GET EXACT TARGETS
# ---------------------------------------------------------
exact_df <- lipid_df2 %>%
  filter(name %in% exact_targets)

missing_exact <- setdiff(exact_targets, exact_df$name)
if (length(missing_exact) > 0) {
  message("\nExact targets not found:")
  print(missing_exact)
}

# ---------------------------------------------------------
# 9. GET PATTERN-BASED REPRESENTATIVE TARGETS
# ---------------------------------------------------------
pattern_matches_list <- list()

for (pat in pattern_targets) {
  n_keep <- n_keep_per_pattern[[pat]]
  if (is.null(n_keep)) next
  
  tmp <- lipid_df2 %>%
    filter(str_detect(name, regex(pat, ignore_case = FALSE))) %>%
    arrange(rank_padj, desc(rank_abs_ratio), name) %>%
    distinct(name, .keep_all = TRUE) %>%
    slice_head(n = n_keep)
  
  if (nrow(tmp) > 0) {
    tmp$match_pattern <- pat
    pattern_matches_list[[pat]] <- tmp
  }
}

pattern_df <- bind_rows(pattern_matches_list)

# ---------------------------------------------------------
# 10. COMBINE, DEDUPLICATE, ORDER
# ---------------------------------------------------------
plot_df <- bind_rows(exact_df, pattern_df) %>%
  distinct(name, .keep_all = TRUE)

preferred_order <- c(
  "PC 34:1",
  "PC 36:7",
  "PC(18:0/Aze)",
  "PC (18:0/Aze)",
  "PC(16:0/8:(COOH))",
  "PC (16:0/8:(COOH))",
  "PC(Azelaoyl-PAF)",
  "DG 30:3",
  "DG 32:5",
  "DG 34:2",
  "PS 34:4"
)

other_names <- setdiff(plot_df$name, preferred_order)
ordered_names <- c(preferred_order[preferred_order %in% plot_df$name], sort(other_names))

plot_df <- plot_df %>%
  mutate(name = factor(name, levels = rev(ordered_names))) %>%
  arrange(name)

# ---------------------------------------------------------
# 11. PRINT WHAT WAS SELECTED
# ---------------------------------------------------------
message("\nLipids selected for the recovery targeted heatmap:\n")

selected_tbl <- plot_df %>%
  select(
    name,
    any_of(padj_cols_found),
    any_of(ratio_cols_found)
  )

print(as.data.frame(selected_tbl), row.names = FALSE)

# ---------------------------------------------------------
# 12. BUILD HEATMAP MATRIX
# ---------------------------------------------------------
heat_mat <- plot_df %>%
  select(
    name,
    all_of(c(
      rdh12_veh_col,
      rdh12_100_col,
      rdh12_200_col,
      wt_veh_col,
      wt_100_col,
      wt_200_col
    ))
  ) %>%
  as.data.frame()

rownames(heat_mat) <- as.character(heat_mat$name)
heat_mat$name <- NULL
heat_mat <- as.matrix(heat_mat)

colnames(heat_mat) <- c(
  "RDH12 Veh",
  "RDH12 100",
  "RDH12 200",
  "WT Veh",
  "WT 100",
  "WT 200"
)

# ---------------------------------------------------------
# 13. OPTIONAL ROW ANNOTATION BY CLASS
# ---------------------------------------------------------
# ---------------------------------------------------------
# 13. OPTIONAL ROW ANNOTATION BY CLASS
# ---------------------------------------------------------
get_lipid_class <- function(x) {
  case_when(
    str_detect(x, "^PC O-") ~ "PC ether",
    str_detect(x, "^PC P-") ~ "PC plasmalogen",
    str_detect(x, "^PC\\s") ~ "PC",
    str_detect(x, "^PE\\s") ~ "PE",
    str_detect(x, "^PI\\s") ~ "PI",
    str_detect(x, "^DG\\s") ~ "DG",
    str_detect(x, "^PS\\s") ~ "PS",
    str_detect(x, "^TG\\s") ~ "TG",
    str_detect(x, "^CE\\s") ~ "CE",
    str_detect(x, regex("Aze|COOH|Azelaoyl", ignore_case = TRUE)) ~ "Oxidized PC",
    TRUE ~ "Other"
  )
}

annotation_row <- data.frame(
  Class = get_lipid_class(rownames(heat_mat)),
  row.names = rownames(heat_mat),
  stringsAsFactors = FALSE
)

annotation_colors <- list(
  Class = c(
    "PC" = "#E9C46A",
    "PC ether" = "#DDB892",
    "PC plasmalogen" = "#C68B59",
    "Oxidized PC" = "#E76F51",
    "PE" = "#5DADE2",
    "PI" = "#4A90E2",
    "DG" = "#C77DFF",
    "PS" = "#A1C181",
    "TG" = "#90BE6D",
    "CE" = "#F4A261",
    "Other" = "grey70"
  )
)

# ---------------------------------------------------------
# 14. COLUMN ANNOTATION (GENOTYPE / DOSE)
# ---------------------------------------------------------
annotation_col <- data.frame(
  Genotype = c("RDH12", "RDH12", "RDH12", "WT", "WT", "WT"),
  Treatment = c("Veh", "100", "200", "Veh", "100", "200"),
  row.names = colnames(heat_mat),
  stringsAsFactors = FALSE
)

annotation_col_colors <- list(
  Genotype = c("RDH12" = "#6A3D9A", "WT" = "#1F78B4"),
  Treatment = c("Veh" = "grey70", "100" = "#FDBF6F", "200" = "#E31A1C")
)

# ---------------------------------------------------------
# 15. DRAW HEATMAP
# ---------------------------------------------------------
# Make color scale symmetric around zero so 0 = white
max_abs <- max(abs(heat_mat), na.rm = TRUE)
my_breaks <- seq(-max_abs, max_abs, length.out = 101)
my_colors <- colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100)

draw_heatmap <- function() {
  pheatmap(
    mat = heat_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    border_color = NA,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    annotation_colors = c(annotation_colors, annotation_col_colors),
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 45,
    main = "Recovery targeted lipid heatmap",
    color = my_colors,
    breaks = my_breaks
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 16. SAVE OUTPUTS
# ---------------------------------------------------------
pdf(out_pdf, width = 9, height = 9)
draw_heatmap()
dev.off()

png(out_png, width = 2600, height = 2600, res = 300)
draw_heatmap()
dev.off()

write.csv(
  plot_df %>%
    select(
      name,
      all_of(c(
        rdh12_veh_col,
        rdh12_100_col,
        rdh12_200_col,
        wt_veh_col,
        wt_100_col,
        wt_200_col
      )),
      any_of(padj_cols_found),
      any_of(ratio_cols_found)
    ),
  out_csv,
  row.names = FALSE
)

cat("\nSaved files:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
cat(out_csv, "\n")