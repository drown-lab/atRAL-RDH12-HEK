# BiocManager::install("ComplexHeatmap")
# BiocManager::install("circlize")

library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

se <- data_filt   # or data_se

mat  <- assay(se)
meta <- as.data.frame(colData(se)) |> tibble::rownames_to_column("sample")

# ---- Check required columns exist (edit these if your names differ)
needed <- c("Genetype", "ExpType", "Treatment", "replicate")
stopifnot(all(needed %in% colnames(meta)))

# ---- 1) Order columns by Genetype -> ExpType -> Treatment -> replicate
meta_ord <- meta |>
  mutate(
    replicate = as.character(replicate),
    sample    = as.character(sample)
  ) |>
  arrange(Genetype, ExpType, Treatment, replicate)

# Reorder matrix columns to match
mat_ord <- mat[, meta_ord$sample, drop = FALSE]

# ---- 2) Missingness matrix (0=detected, 1=missing)
miss_num <- ifelse(is.na(mat_ord), 1, 0)

# Optional: keep top N most-missing proteins to keep it readable
top_n <- 750
miss_frac <- rowMeans(miss_num == 1)
keep_rows <- order(miss_frac, decreasing = TRUE)[seq_len(min(top_n, nrow(miss_num)))]
miss_num  <- miss_num[keep_rows, , drop = FALSE]

# ---- 3) Colors for annotations
gen_cols <- setNames(scales::hue_pal()(length(unique(meta_ord$Genetype))),
                     unique(meta_ord$Genetype))
exp_cols <- setNames(scales::hue_pal()(length(unique(meta_ord$ExpType))),
                     unique(meta_ord$ExpType))
trt_cols <- setNames(scales::hue_pal()(length(unique(meta_ord$Treatment))),
                     unique(meta_ord$Treatment))
rep_cols <- c(a = "#1f78b4", b = "#33a02c", c = "#e31a1c")  # adjust if needed

# ---- 4) Build 4-layer annotation (this is your “3-way layering” + rep)
ha <- HeatmapAnnotation(
  Genetype  = meta_ord$Genetype,
  ExpType   = meta_ord$ExpType,
  Treatment = meta_ord$Treatment,
  Rep       = meta_ord$replicate,
  col = list(
    Genetype  = gen_cols,
    ExpType   = exp_cols,
    Treatment = trt_cols,
    Rep       = rep_cols
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)

# ---- 5) Short x labels: ONLY show replicate letters (a/b/c)
# If you want NO labels, set show_column_names = FALSE below.
col_labels <- meta_ord$replicate

# ---- 6) Heatmap with clustering on rows (proteins)
ht <- Heatmap(
  miss_num,
  name = "Missing",
  col = c("0" = "#FFFFF0", "1" = "deepskyblue4"),
  top_annotation = ha,
  
  cluster_rows = TRUE,        # dendrogram "bracket" on left
  cluster_columns = FALSE,    # keep your Genetype->ExpType->Treatment->rep order
  show_row_dend = FALSE,
  show_row_names = FALSE,     # do NOT print protein names
  show_column_names = TRUE,   # keep x labels (short!)
  column_labels = col_labels,
  column_names_rot = 0,
  
  column_title = "Missingness (Detected=white, Missing=grey)",
  heatmap_legend_param = list(at = c(0, 1), labels = c("Detected", "Missing"))
)

draw(ht)
