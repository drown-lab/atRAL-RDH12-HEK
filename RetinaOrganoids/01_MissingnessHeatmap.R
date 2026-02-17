library(dplyr)
library(tidyr)
library(stringr)

# 0 = detected, 1 = missing
miss_df <- Orgnoid_v5 %>%
  group_by(Protein.Group, SampleName) %>%
  summarise(
    detected = any(!is.na(PG.MaxLFQ)),
   # detected = any(!is.na(PG.MaxLFQ) | PG.MaxLFQ > 0)
    
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

miss_num <- miss_wide %>%
  select(-Protein.Group) %>%
  as.matrix()

rownames(miss_num) <- miss_wide$Protein.Group
storage.mode(miss_num) <- "numeric"

sample_order <- c("EP_F1", "EP_F2", "EP_F3", "EP_F4", "WO_F1", "WO_F2", "WO_F3", "WO_F4")

sample_order <- sample_order[sample_order %in% colnames(miss_num)]
miss_num <- miss_num[, sample_order, drop = FALSE]


col_labels <- colnames(miss_num)

###Build Annotation
library(ComplexHeatmap)
library(circlize)

sample_info <- Orgnoid_v5 %>%
  distinct(SampleName, CellType, Fraction) %>%
  filter(SampleName %in% colnames(miss_num)) %>%
  mutate(SampleName = factor(SampleName, levels = colnames(miss_num))) %>%
  arrange(SampleName)

ha <- HeatmapAnnotation(
  CellType  = sample_info$CellType,
  Fraction   = sample_info$Fraction,
    annotation_name_side = "left"
)

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
  
  column_title = "Missingness (Detected=white, Missing=blue)",
  heatmap_legend_param = list(at = c(0, 1), labels = c("Detected", "Missing"))
)

draw(ht)