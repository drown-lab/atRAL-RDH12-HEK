# =========================================================
# Acute targeted lipid heatmap for main-text figure
# Uses: Acute_DEA_results_v1_corrected.csv
# Goal: compact heatmap of biologically relevant differential lipids
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
acute_file <- "Lipidomics/output_txts/Acute_DEA_results_v1_corrected.csv"

out_dir <- "Lipidomics/Figures/Acute/Acute_targeted_heatmap/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Acute_targeted_lipid_heatmap.pdf")
out_png <- file.path(out_dir, "Acute_targeted_lipid_heatmap.png")
out_csv <- file.path(out_dir, "Acute_targeted_lipid_heatmap_table.csv")

# ---------------------------------------------------------
# 2. READ DATA
# ---------------------------------------------------------
lipid_df <- read.csv(
  acute_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "name",
  "RDH12_Veh_centered",
  "RDH12_100_centered",
  "RDH12_200_centered",
  "RDH12_100_vs_RDH12_Veh_p.adj",
  "RDH12_200_vs_RDH12_Veh_p.adj",
  "RDH12_100_vs_RDH12_Veh_ratio",
  "RDH12_200_vs_RDH12_Veh_ratio"
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
# These are the lipids/classes you said are most relevant for the acute figure.
# Keep exact names where possible. Use regex patterns for representative classes.

exact_targets <- c(
  "LPC 18:1",
  "PC 32:0",
  "PC 28:0",
  "DG 32:0",
  "DG 34:1",
  "PE 34:1",
  "PE 36:4",
  "PE 40:5"
)

# Pattern-based additions for representative species if exact names differ
# You can edit these after seeing the printed matches.
pattern_targets <- c(
  "^PI\\s",             # any PI species
  "^Cer\\s",            # ceramide species
  "^dhCer\\s",          # dihydroceramide species
  "^TG\\s",             # TG species
  "^PC O-",             # ether-PC / plasmalogen-like entries
  "^PC P-",             # plasmalogen naming if present
  "^LPC\\s",            # LPC entries
  "^DG\\s"              # DG entries
)

# How many additional representatives to keep from each pattern group
# ranked by smallest 200 vs Veh adjusted p-value, then smallest 100 vs Veh adjusted p-value
n_keep_per_pattern <- list(
  "^PI\\s"   = 2,
  "^Cer\\s"  = 2,
  "^dhCer\\s"= 2,
  "^TG\\s"   = 2,
  "^PC O-"   = 1,
  "^PC P-"   = 1,
  "^LPC\\s"  = 1,
  "^DG\\s"   = 2
)

# ---------------------------------------------------------
# 4. PREP AND RANK LIPIDS
# ---------------------------------------------------------
lipid_df2 <- lipid_df %>%
  mutate(
    name = str_trim(name),
    rank_padj = pmin(
      `RDH12_100_vs_RDH12_Veh_p.adj`,
      `RDH12_200_vs_RDH12_Veh_p.adj`,
      na.rm = TRUE
    ),
    rank_abs_ratio = pmax(
      abs(`RDH12_100_vs_RDH12_Veh_ratio`),
      abs(`RDH12_200_vs_RDH12_Veh_ratio`),
      na.rm = TRUE
    )
  )

# ---------------------------------------------------------
# 5. GET EXACT TARGETS
# ---------------------------------------------------------
exact_df <- lipid_df2 %>%
  filter(name %in% exact_targets)

# Warn if some exact targets were not found
missing_exact <- setdiff(exact_targets, exact_df$name)
if (length(missing_exact) > 0) {
  message("Exact targets not found:")
  print(missing_exact)
}

# ---------------------------------------------------------
# 6. GET PATTERN-BASED REPRESENTATIVE TARGETS
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
# 7. COMBINE, DEDUPLICATE, AND OPTIONALLY PRIORITIZE
# ---------------------------------------------------------
plot_df <- bind_rows(exact_df, pattern_df) %>%
  distinct(name, .keep_all = TRUE)

# Optional manual ordering for a cleaner pathway-aware figure
preferred_order <- c(
  "LPC 18:1",
  "PC 32:0",
  "PC 28:0",
  "DG 32:0",
  "DG 34:1",
  "PE 34:1",
  "PE 36:4",
  "PE 40:5"
)

# add any automatically selected representatives after the preferred ones
other_names <- setdiff(plot_df$name, preferred_order)
ordered_names <- c(preferred_order[preferred_order %in% plot_df$name], sort(other_names))

plot_df <- plot_df %>%
  mutate(name = factor(name, levels = rev(ordered_names))) %>%
  arrange(name)

# ---------------------------------------------------------
# 8. PRINT WHAT WAS SELECTED
# ---------------------------------------------------------
message("\nLipids selected for the acute targeted heatmap:\n")

selected_tbl <- plot_df %>%
  select(
    name,
    `RDH12_100_vs_RDH12_Veh_p.adj`,
    `RDH12_200_vs_RDH12_Veh_p.adj`,
    `RDH12_100_vs_RDH12_Veh_ratio`,
    `RDH12_200_vs_RDH12_Veh_ratio`
  )

print(as.data.frame(selected_tbl), row.names = FALSE)

# ---------------------------------------------------------
# 9. BUILD HEATMAP MATRIX
# ---------------------------------------------------------
heat_mat <- plot_df %>%
  select(
    name,
    `RDH12_Veh_centered`,
    `RDH12_100_centered`,
    `RDH12_200_centered`
  ) %>%
  as.data.frame()

rownames(heat_mat) <- as.character(heat_mat$name)
heat_mat$name <- NULL
heat_mat <- as.matrix(heat_mat)

colnames(heat_mat) <- c("Veh", "100 µM", "200 µM")

# ---------------------------------------------------------
# 10. OPTIONAL ROW ANNOTATION BY CLASS
# ---------------------------------------------------------
get_lipid_class <- function(x) {
  case_when(
    str_detect(x, "^LPC\\s")   ~ "LPC",
    str_detect(x, "^PC O-")    ~ "PC ether",
    str_detect(x, "^PC P-")    ~ "PC plasmalogen",
    str_detect(x, "^PC\\s")    ~ "PC",
    str_detect(x, "^DG\\s")    ~ "DG",
    str_detect(x, "^PE\\s")    ~ "PE",
    str_detect(x, "^PI\\s")    ~ "PI",
    str_detect(x, "^Cer\\s")   ~ "Cer",
    str_detect(x, "^dhCer\\s") ~ "dhCer",
    str_detect(x, "^TG\\s")    ~ "TG",
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
    "LPC" = "#F4A261",
    "PC" = "#E9C46A",
    "PC ether" = "#DDB892",
    "PC plasmalogen" = "#C68B59",
    "DG" = "#C77DFF",
    "PE" = "#5DADE2",
    "PI" = "#4A90E2",
    "Cer" = "#E76F51",
    "dhCer" = "#F28482",
    "TG" = "#90BE6D",
    "Other" = "grey70"
  )
)

# ---------------------------------------------------------
# 11. DRAW HEATMAP
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
    annotation_colors = annotation_colors,
    fontsize_row = 11,
    fontsize_col = 11,
    angle_col = 0,
    main = "Acute targeted lipid heatmap",
    color = my_colors,
    breaks = my_breaks
  )
}

draw_heatmap()
# ---------------------------------------------------------
# 12. SAVE OUTPUTS
# ---------------------------------------------------------
pdf(out_pdf, width = 7.5, height = 8.5)
draw_heatmap()
dev.off()

png(out_png, width = 2200, height = 2500, res = 900)
draw_heatmap()
dev.off()

write.csv(
  plot_df %>%
    select(
      name,
      `RDH12_Veh_centered`,
      `RDH12_100_centered`,
      `RDH12_200_centered`,
      `RDH12_100_vs_RDH12_Veh_p.adj`,
      `RDH12_200_vs_RDH12_Veh_p.adj`,
      `RDH12_100_vs_RDH12_Veh_ratio`,
      `RDH12_200_vs_RDH12_Veh_ratio`
    ),
  out_csv,
  row.names = FALSE
)

cat("\nSaved files:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
cat(out_csv, "\n")