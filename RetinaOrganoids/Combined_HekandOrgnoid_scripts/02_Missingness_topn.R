# --- Complete code: build + plot TOP N most-missing proteins heatmap ---

library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

top_n <- 2000  # how many proteins (rows) to display

# 0 = detected, 1 = missing
miss_df <- combined_proteome %>%
  group_by(Protein.Group, SampleName) %>%
  summarise(
    # detected if at least one positive, non-missing value exists in that group
    detected = any(!is.na(PG.MaxLFQ) & PG.MaxLFQ > 0),
    .groups = "drop"
  ) %>%
  mutate(Missing = as.integer(!detected)) %>%
  select(Protein.Group, SampleName, Missing)

miss_wide <- miss_df %>%
  pivot_wider(
    names_from = SampleName,
    values_from = Missing,
    values_fill = 1L   # if a protein has no row for a sample, treat as missing
  )

# ---- compute TOP N most-missing proteins ----
top_missing_tbl <- miss_wide %>%
  mutate(
    n_missing = rowSums(across(-Protein.Group)),
    n_detected = (ncol(miss_wide) - 1) - n_missing,
    missing_pct = n_missing / (ncol(miss_wide) - 1) * 100
  ) %>%
  arrange(desc(n_missing), Protein.Group) %>%
  slice_head(n = top_n)

top_ids <- top_missing_tbl$Protein.Group

# ---- build numeric matrix (full), then subset to top ----
miss_num <- miss_wide %>%
  select(-Protein.Group) %>%
  as.matrix()

rownames(miss_num) <- miss_wide$Protein.Group
storage.mode(miss_num) <- "numeric"

# ---- column order (keep only those that exist) ----
sample_order <- c(
  "EP_F1", "EP_F2", "EP_F3", "EP_F4",
  "WO_F1", "WO_F2", "WO_F3", "WO_F4",
  "GFP_1", "GFP_2", "GFP_3",
  "RDH12-2", "RDH12_1", "RDH12_3"
)

sample_order <- sample_order[sample_order %in% colnames(miss_num)]
miss_num <- miss_num[, sample_order, drop = FALSE]

# labels after ordering
col_labels <- colnames(miss_num)

# ---- subset to TOP N proteins ----
miss_num_top <- miss_num[rownames(miss_num) %in% top_ids, , drop = FALSE]

# Optional: keep rows in the same order as top_missing_tbl
miss_num_top <- miss_num_top[top_ids[top_ids %in% rownames(miss_num_top)], , drop = FALSE]

# ---- annotation data (match the plotted columns) ----
sample_info <- combined_proteome %>%
  distinct(SampleName, CellType, Fraction) %>%
  filter(SampleName %in% colnames(miss_num_top)) %>%
  mutate(SampleName = factor(SampleName, levels = colnames(miss_num_top))) %>%
  arrange(SampleName)

ha <- HeatmapAnnotation(
  CellType = sample_info$CellType,
  Fraction = sample_info$Fraction,
  annotation_name_side = "left"
)

# ---- plot heatmap (TOP N) ----
ht <- Heatmap(
  miss_num_top,
  name = "Missing",
  col = c("0" = "#FFFFF0", "1" = "deepskyblue4"),
  top_annotation = ha,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  
  # With TOP N, it can be useful to show names; set FALSE if too crowded
  show_row_names = FALSE,
  
  show_column_names = TRUE,
  column_labels = col_labels,
  column_names_rot = 0,
  
  column_title = paste0(
    "Top ", top_n, " most-missing proteins (Detected=white, Missing=blue)"
  ),
  heatmap_legend_param = list(at = c(0, 1), labels = c("Detected", "Missing"))
)

draw(ht)

# ---- Optional: inspect the top-missing table ----
# View(top_missing_tbl)
# print(top_missing_tbl %>% select(Protein.Group, n_missing, n_detected, missing_pct))
