# =========================================================
# Recovery targeted lipid heatmap for main-text figure
# Uses: Recovery_DEA_results_v1_corrected.csv
# Goal:
#   - compact heatmap of biologically relevant recovery lipids
#   - lipid classes grouped next to each other
#   - cleaner TG/DG labels without NL fragments
#   - zero centered at white in color scale
#   - shows WT/Control and RDH12 recovery conditions
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(pheatmap)
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

required_cols <- c(
  "name",
  
  # centered abundance columns
  "Control_Veh_centered",
  "Control_100_centered",
  "Control_200_centered",
  "RDH12_Veh_centered",
  "RDH12_100_centered",
  "RDH12_200_centered",
  
  # within WT/Control treatment comparisons
  "Control_100_vs_Control_Veh_p.adj",
  "Control_200_vs_Control_Veh_p.adj",
  "Control_100_vs_Control_Veh_ratio",
  "Control_200_vs_Control_Veh_ratio",
  
  # within RDH12 treatment comparisons
  "RDH12_100_vs_RDH12_Veh_p.adj",
  "RDH12_200_vs_RDH12_Veh_p.adj",
  "RDH12_100_vs_RDH12_Veh_ratio",
  "RDH12_200_vs_RDH12_Veh_ratio",
  
  # RDH12 vs WT/Control comparisons
  "RDH12_Veh_vs_Control_Veh_p.adj",
  "RDH12_100_vs_Control_100_p.adj",
  "RDH12_200_vs_Control_200_p.adj",
  "RDH12_Veh_vs_Control_Veh_ratio",
  "RDH12_100_vs_Control_100_ratio",
  "RDH12_200_vs_Control_200_ratio"
)

missing_cols <- setdiff(required_cols, colnames(lipid_df))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 3. USER-DEFINED TARGETS
# ---------------------------------------------------------
# These targets are based on the recovery lipid story:
# - RDH12 recovery has increased PC / PC O- species
# - oxidized PCs appear/persist in recovery
# - DG species are lower in RDH12 relative to WT
# - PE/PI depletion persists after recovery
# - PS 34:4 is enriched in WT relative to RDH12

exact_targets <- c(
  # PC / membrane restoration
  "PC 34:1",
  "PC 36:7",
  "PC 32:0",
  "PC 28:0",
  
  # oxidized PC species / secondary lipid oxidation products
  "PC(18:0/Aze)",
  "PC (18:0/Aze)",
  "PC(16:0/8:0(COOH))",
  "PC(16:0/8:(COOH))",
  "PC (16:0/8:0(COOH))",
  "PC (16:0/8:(COOH))",
  "PC(Azelaoyl-PAF)",
  "PC (Azelaoyl-PAF)",
  
  # DG remodeling / lower in RDH12 relative to WT
  "DG 30:3",
  "DG 32:5",
  "DG 34:2",
  "DG 32:0",
  "DG 34:1",
  "DG 38:0 NL 20:0",
  
  # PE / PI depletion or persistent remodeling
  "PE 34:1",
  "PE 36:4",
  "PE 40:5",
  "PI 38:5",
  
  # PS enriched in WT relative to RDH12
  "PS 34:4",
  
  # LPC / phospholipid turnover
  "LPC 14:0",
  "LPC 16:0",
  "LPC 16:1",
  "LPC 18:1",
  
  # ether-PC / plasmalogen examples
  "PC O-30:0",
  "PC O-34:1",
  "PC O-36:4",
  "PC O-36:5",
  "PC O-36:6",
  
  # sphingolipid stress examples
  "Cer(d18:0/20:0(2OH))",
  "Cer(d14:2/16:0)",
  "Cer(d14:1/16:0)"
)

# Pattern-based additions for representative species if exact names differ.
# These are ranked by strongest adjusted p-value across recovery contrasts.
pattern_targets <- c(
  "^PC O-",
  "^PC P-",
  "^PC\\s",
  "Aze|COOH|Azelaoyl",
  "^DG\\s",
  "^PE\\s",
  "^PI\\s",
  "^PS\\s",
  "^LPC\\s",
  "^Cer\\s",
  "^dhCer\\s",
  "^SM\\s",
  "^TG\\s",
  "^\\[TG"
)

n_keep_per_pattern <- list(
  "^PC O-"          = 5,
  "^PC P-"          = 2,
  "^PC\\s"          = 4,
  "Aze|COOH|Azelaoyl" = 3,
  "^DG\\s"          = 4,
  "^PE\\s"          = 3,
  "^PI\\s"          = 3,
  "^PS\\s"          = 2,
  "^LPC\\s"         = 2,
  "^Cer\\s"         = 2,
  "^dhCer\\s"       = 2,
  "^SM\\s"          = 2,
  "^TG\\s"          = 1,
  "^\\[TG"          = 1
)

# ---------------------------------------------------------
# 4. HELPER FUNCTIONS
# ---------------------------------------------------------

clean_lipid_label <- function(x) {
  x %>%
    str_trim() %>%
    # Convert bracketed TG labels like [TG44:5] NL 20:0 to TG 44:5
    str_replace("^\\[TG([0-9]+:[0-9]+)\\].*$", "TG \\1") %>%
    # Remove neutral-loss text from DG/TG labels like DG 38:0 NL 20:0
    str_replace("\\s+NL\\s+[0-9]+:[0-9]+.*$", "") %>%
    str_replace("\\s+", " ") %>%
    str_trim()
}

get_lipid_class <- function(x) {
  case_when(
    str_detect(x, regex("Aze|COOH|Azelaoyl", ignore_case = TRUE)) ~ "Oxidized PC",
    str_detect(x, "^\\[TG")    ~ "TG",
    str_detect(x, "^TG\\s")    ~ "TG",
    str_detect(x, "^LPC\\s")   ~ "LPC",
    str_detect(x, "^PC O-")    ~ "PC ether",
    str_detect(x, "^PC P-")    ~ "PC plasmalogen",
    str_detect(x, "^PC\\s")    ~ "PC",
    str_detect(x, "^DG\\s")    ~ "DG",
    str_detect(x, "^PE\\s")    ~ "PE",
    str_detect(x, "^PI\\s")    ~ "PI",
    str_detect(x, "^PS\\s")    ~ "PS",
    str_detect(x, "^Cer\\s")   ~ "Cer",
    str_detect(x, "^dhCer\\s") ~ "dhCer",
    str_detect(x, "^SM\\s")    ~ "SM",
    TRUE ~ "Other"
  )
}

# Class order controls how lipid classes appear together
class_order <- c(
  "LPC",
  "PC",
  "PC ether",
  "PC plasmalogen",
  "Oxidized PC",
  "DG",
  "PE",
  "PI",
  "PS",
  "Cer",
  "dhCer",
  "SM",
  "TG",
  "Other"
)

# ---------------------------------------------------------
# 5. PREP AND RANK LIPIDS
# ---------------------------------------------------------

padj_cols_for_rank <- c(
  "Control_100_vs_Control_Veh_p.adj",
  "Control_200_vs_Control_Veh_p.adj",
  "RDH12_100_vs_RDH12_Veh_p.adj",
  "RDH12_200_vs_RDH12_Veh_p.adj",
  "RDH12_Veh_vs_Control_Veh_p.adj",
  "RDH12_100_vs_Control_100_p.adj",
  "RDH12_200_vs_Control_200_p.adj"
)

ratio_cols_for_rank <- c(
  "Control_100_vs_Control_Veh_ratio",
  "Control_200_vs_Control_Veh_ratio",
  "RDH12_100_vs_RDH12_Veh_ratio",
  "RDH12_200_vs_RDH12_Veh_ratio",
  "RDH12_Veh_vs_Control_Veh_ratio",
  "RDH12_100_vs_Control_100_ratio",
  "RDH12_200_vs_Control_200_ratio"
)

lipid_df2 <- lipid_df %>%
  mutate(
    name = str_trim(name),
    DisplayName = clean_lipid_label(name),
    Class = factor(get_lipid_class(name), levels = class_order)
  )

padj_matrix <- as.matrix(lipid_df2[, padj_cols_for_rank, drop = FALSE])
ratio_matrix <- as.matrix(lipid_df2[, ratio_cols_for_rank, drop = FALSE])

lipid_df2$rank_padj <- apply(padj_matrix, 1, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
})

lipid_df2$rank_abs_ratio <- apply(ratio_matrix, 1, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else max(abs(x), na.rm = TRUE)
})

# ---------------------------------------------------------
# 6. GET EXACT TARGETS
# ---------------------------------------------------------
exact_df <- lipid_df2 %>%
  filter(name %in% exact_targets)

missing_exact <- setdiff(exact_targets, exact_df$name)

if (length(missing_exact) > 0) {
  message("Exact targets not found:")
  print(missing_exact)
}

# ---------------------------------------------------------
# 7. GET PATTERN-BASED REPRESENTATIVE TARGETS
# ---------------------------------------------------------
pattern_matches_list <- list()

for (pat in pattern_targets) {
  n_keep <- n_keep_per_pattern[[pat]]
  if (is.null(n_keep)) next
  
  tmp <- lipid_df2 %>%
    filter(str_detect(name, regex(pat, ignore_case = TRUE))) %>%
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
# 8. COMBINE, DEDUPLICATE, AND ORDER BY CLASS
# ---------------------------------------------------------
plot_df <- bind_rows(exact_df, pattern_df) %>%
  distinct(name, .keep_all = TRUE) %>%
  mutate(
    Class = factor(as.character(Class), levels = class_order)
  ) %>%
  arrange(Class, DisplayName)

# Make display names unique in case cleaned labels collide
plot_df$DisplayName <- make.unique(plot_df$DisplayName, sep = " ")

# ---------------------------------------------------------
# 9. PRINT WHAT WAS SELECTED
# ---------------------------------------------------------
message("\nLipids selected for the recovery targeted heatmap:\n")

selected_tbl <- plot_df %>%
  select(
    name,
    DisplayName,
    Class,
    all_of(padj_cols_for_rank),
    all_of(ratio_cols_for_rank)
  )

print(as.data.frame(selected_tbl), row.names = FALSE)

# ---------------------------------------------------------
# 10. BUILD HEATMAP MATRIX
# ---------------------------------------------------------
heat_mat <- plot_df %>%
  select(
    DisplayName,
    `RDH12_Veh_centered`,
    `RDH12_100_centered`,
    `RDH12_200_centered`,
    `Control_Veh_centered`,
    `Control_100_centered`,
    `Control_200_centered`
  ) %>%
  as.data.frame()

rownames(heat_mat) <- heat_mat$DisplayName
heat_mat$DisplayName <- NULL
heat_mat <- as.matrix(heat_mat)

colnames(heat_mat) <- c(
  "RDH12 Veh",
  "RDH12 100 µM",
  "RDH12 200 µM",
  "WT Veh",
  "WT 100 µM",
  "WT 200 µM"
)

# ---------------------------------------------------------
# 11. ROW ANNOTATION BY CLASS
# ---------------------------------------------------------
annotation_row <- data.frame(
  Class = plot_df$Class,
  row.names = plot_df$DisplayName,
  stringsAsFactors = FALSE
)

annotation_colors <- list(
  Class = c(
    "LPC" = "#F4A261",
    "PC" = "#E9C46A",
    "PC ether" = "#DDB892",
    "PC plasmalogen" = "#C68B59",
    "Oxidized PC" = "#E76F51",
    "DG" = "#C77DFF",
    "PE" = "#5DADE2",
    "PI" = "#4A90E2",
    "PS" = "#A1C181",
    "Cer" = "#E76F51",
    "dhCer" = "#F28482",
    "SM" = "#8E7DBE",
    "TG" = "#90BE6D",
    "Other" = "grey70"
  )
)

# ---------------------------------------------------------
# 12. GAPS BETWEEN LIPID CLASSES
# ---------------------------------------------------------
class_counts <- table(factor(annotation_row$Class, levels = class_order))
class_counts <- class_counts[class_counts > 0]

gaps_row <- cumsum(class_counts)
gaps_row <- gaps_row[-length(gaps_row)]

# ---------------------------------------------------------
# 13. COLOR SCALE CENTERED AT ZERO
# ---------------------------------------------------------
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
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    border_color = NA,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    gaps_row = gaps_row,
    fontsize_row = 10,
    fontsize_col = 11,
    angle_col = 45,
    main = "Recovery targeted lipid heatmap",
    color = my_colors,
    breaks = my_breaks,
    legend_breaks = legend_breaks,
    legend_labels = legend_labels
  )
}

draw_heatmap()

# ---------------------------------------------------------
# 15. SAVE OUTPUTS
# ---------------------------------------------------------
pdf(out_pdf, width = 8.5, height = 10)
draw_heatmap()
dev.off()

png(out_png, width = 2600, height = 3000, res = 900)
draw_heatmap()
dev.off()

write.csv(
  plot_df %>%
    select(
      name,
      DisplayName,
      Class,
      `RDH12_Veh_centered`,
      `RDH12_100_centered`,
      `RDH12_200_centered`,
      `Control_Veh_centered`,
      `Control_100_centered`,
      `Control_200_centered`,
      all_of(padj_cols_for_rank),
      all_of(ratio_cols_for_rank)
    ),
  out_csv,
  row.names = FALSE
)

cat("\nSaved files:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
cat(out_csv, "\n")