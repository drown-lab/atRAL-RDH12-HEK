# BiocManager::install("ComplexHeatmap")
# BiocManager::install("circlize")

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

## Get PPFtable_v4 from the PPF prep script

input_dir <- if (dir.exists("PPF_datasets")) "PPF_datasets" else "."
prep_script_path <- file.path(input_dir, "01_Prep_data.R")
path_figures <- file.path(input_dir, "output_missingness")

if (!file.exists(prep_script_path)) {
  prep_script_path <- "01_Prep_data.R"
}

if (!file.exists(prep_script_path)) {
  stop("Could not find 01_Prep_data.R. Run this script from the project root or PPF_datasets folder.")
}

dir.create(path_figures, showWarnings = FALSE, recursive = TRUE)

source(prep_script_path)

exists("PPFtable_v4")
colnames(PPFtable_v4)

# ---- Check required columns exist
needed <- c("sample_name", "Run", "Protein.Group", "Genes", "PG.MaxLFQ")
stopifnot(all(needed %in% colnames(PPFtable_v4)))

# ---- 1) Build abundance matrix and sample metadata
PPF_missingness_input <- PPFtable_v4 |>
  mutate(
    sample = as.character(coalesce(sample_name, Run)),
    ProteinLabel = if_else(
      is.na(Genes) | Genes == "",
      Protein.Group,
      paste(Genes, Protein.Group, sep = "_")
    ),
    Condition = str_remove(sample, "\\d+$"),
    replicate = str_extract(sample, "\\d+$")
  ) |>
  mutate(
    Condition = if_else(is.na(Condition) | Condition == "", sample, Condition),
    replicate = if_else(is.na(replicate) | replicate == "", sample, replicate)
  ) |>
  group_by(Protein.Group, Genes, ProteinLabel, sample, Condition, replicate) |>
  summarise(
    PG.MaxLFQ = if_else(all(is.na(PG.MaxLFQ)), NA_real_, max(PG.MaxLFQ, na.rm = TRUE)),
    .groups = "drop"
  )

meta <- PPF_missingness_input |>
  distinct(sample, Condition, replicate) |>
  mutate(
    sample = as.character(sample),
    replicate = as.character(replicate)
  )

# ---- Check required metadata columns exist
needed_meta <- c("Condition", "replicate")
stopifnot(all(needed_meta %in% colnames(meta)))

# ---- Order columns by Condition -> replicate
meta_ord <- meta |>
  arrange(Condition, replicate)

abundance_table <- PPF_missingness_input |>
  select(ProteinLabel, sample, PG.MaxLFQ) |>
  distinct() |>
  pivot_wider(
    names_from = sample,
    values_from = PG.MaxLFQ
  )

mat <- abundance_table |>
  tibble::column_to_rownames("ProteinLabel") |>
  as.matrix()

mode(mat) <- "numeric"

# Reorder matrix columns to match
mat_ord <- mat[, meta_ord$sample, drop = FALSE]

# ---- 2) Missingness matrix (0=detected, 1=missing)
miss_num <- ifelse(is.na(mat_ord) | mat_ord <= 0, 1, 0)

# Optional: keep top N most-missing proteins to keep it readable
top_n <- 500
miss_frac <- rowMeans(miss_num == 1)
keep_rows <- order(miss_frac, decreasing = TRUE)[seq_len(min(top_n, nrow(miss_num)))]
miss_num <- miss_num[keep_rows, , drop = FALSE]

# ---- 3) Colors for annotations
condition_cols <- setNames(scales::hue_pal()(length(unique(meta_ord$Condition))),
                           unique(meta_ord$Condition))
rep_cols <- setNames(scales::hue_pal()(length(unique(meta_ord$replicate))),
                     unique(meta_ord$replicate))

# ---- 4) Build annotation layers
ha <- HeatmapAnnotation(
  Condition = meta_ord$Condition,
  Rep       = meta_ord$replicate,
  col = list(
    Condition = condition_cols,
    Rep       = rep_cols
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)

# ---- 5) Short x labels: ONLY show replicate values
# If you want NO labels, set show_column_names = FALSE below.
col_labels <- meta_ord$replicate

# ---- 6) Heatmap with clustering on rows (proteins)
ht <- Heatmap(
  miss_num,
  name = "Missing Proteins",
  col = c("0" = "#FFFFF0", "1" = "deepskyblue4"),
  top_annotation = ha,
  
  cluster_rows = TRUE,        # dendrogram "bracket" on left
  cluster_columns = FALSE,    # keep your Condition -> replicate order
  show_row_dend = TRUE,
  show_row_names = FALSE,     # do NOT print protein names
  show_column_names = TRUE,   # keep x labels short
  column_labels = col_labels,
  column_names_rot = 0,
  
  column_title = paste0("Top ", min(top_n, nrow(miss_num)), " Missing Proteins Across PPF Samples"),
  heatmap_legend_param = list(at = c(0, 1), labels = c("Detected", "Missing"))
)

draw(ht)

png(
  filename = file.path(path_figures, "HeatMap_missingness_PPF.png"),
  width = 4,
  height = 6,
  units = "in",
  res = 900
)

draw(ht, merge_legend = TRUE)

dev.off()

write.csv(
  as.data.frame(miss_num) |>
    tibble::rownames_to_column("ProteinLabel"),
  file.path(path_figures, "HeatMap_missingness_PPF_matrix.csv"),
  row.names = FALSE
)
