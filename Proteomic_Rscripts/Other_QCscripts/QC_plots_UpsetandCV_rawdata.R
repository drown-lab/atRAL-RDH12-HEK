# =========================================================
# UpSet plot + coefficient of variation plots
# From RawData_filtered_results.csv
#
# Input:
#   RawData_filtered_results.csv
#
# Outputs:
#   1. UpSet plots of protein detection overlap
#   2. CV plots across replicate LFQ intensities
#   3. CSV tables used for plotting
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ---------------------------------------------------------
# 1. INPUT / OUTPUT
# ---------------------------------------------------------

raw_file <- "Proteomic_output_txts/RawData_filtered_results.csv"

out_dir <- "Proteomic_Figs/QC_upset_CV/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Detection rule for UpSet plot
# A protein is considered detected in a condition if it is present in at least
# this many biological replicates.
min_reps_detected <- 2

# CV calculation rule
# CV is calculated only when at least this many replicate values are present.
min_reps_for_cv <- 2

# Optional cap for CV visualization only.
# This prevents a few extreme CVs from compressing the whole figure.
cv_plot_cap <- 200


# ---------------------------------------------------------
# 2. READ DATA
# ---------------------------------------------------------

raw_df <- read.csv(
  raw_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# 2B. CLEAN COLUMN NAMES
# ---------------------------------------------------------

bad_name_cols <- which(is.na(colnames(raw_df)) | colnames(raw_df) == "")

if (length(bad_name_cols) > 0) {
  raw_df <- raw_df[, -bad_name_cols, drop = FALSE]
}

raw_df <- raw_df %>%
  select(
    -any_of(c("Unnamed: 0", "...1", "X", "X.1"))
  )

colnames(raw_df) <- make.unique(colnames(raw_df))

cat("\nCleaned column names:\n")
print(colnames(raw_df))

cat("\nInput dimensions after column cleanup:\n")
print(dim(raw_df))

# ---------------------------------------------------------
# 3. IDENTIFY PROTEIN ID COLUMNS AND LFQ COLUMNS
# ---------------------------------------------------------

protein_id_col <- dplyr::case_when(
  "Protein.Group" %in% colnames(raw_df) ~ "Protein.Group",
  "Protein.Group.IDs" %in% colnames(raw_df) ~ "Protein.Group.IDs",
  "ID" %in% colnames(raw_df) ~ "ID",
  TRUE ~ NA_character_
)

gene_col <- dplyr::case_when(
  "Genes" %in% colnames(raw_df) ~ "Genes",
  "Gene" %in% colnames(raw_df) ~ "Gene",
  TRUE ~ NA_character_
)

if (is.na(protein_id_col)) {
  stop("Could not find a protein ID column. Expected Protein.Group, Protein.Group.IDs, or ID.")
}

lfq_cols <- grep("^LFQ\\.", colnames(raw_df), value = TRUE)

if (length(lfq_cols) == 0) {
  stop("No LFQ columns were found. Expected columns beginning with 'LFQ.'.")
}

cat("\nNumber of LFQ columns found:", length(lfq_cols), "\n")

# ---------------------------------------------------------
# 4. BUILD SAMPLE METADATA FROM LFQ COLUMN NAMES
# ---------------------------------------------------------

sample_meta <- tibble(
  SampleCol = lfq_cols,
  SampleName = str_remove(lfq_cols, "^LFQ\\.")
) %>%
  mutate(
    Treatment_raw = str_extract(SampleName, "^(100|200|control)"),
    Replicate = str_match(SampleName, "^(100|200|control)_([abc])_")[, 3],
    CellType_raw = str_match(SampleName, "^(100|200|control)_[abc]_(RDH12|GFP)_")[, 3],
    Phase = if_else(str_detect(SampleName, "\\+24h_recvr"), "Recovery", "Acute"),
    Treatment = case_when(
      Treatment_raw == "control" ~ "Veh",
      Treatment_raw == "100" ~ "100 µM",
      Treatment_raw == "200" ~ "200 µM",
      TRUE ~ Treatment_raw
    ),
    CellType = case_when(
      CellType_raw == "GFP" ~ "WT",
      CellType_raw == "RDH12" ~ "RDH12",
      TRUE ~ CellType_raw
    ),
    Condition = paste(Phase, CellType, Treatment)
  )

if (any(is.na(sample_meta$Treatment_raw)) |
    any(is.na(sample_meta$Replicate)) |
    any(is.na(sample_meta$CellType_raw))) {
  warning(
    "Some LFQ column names were not parsed cleanly. Check sample_meta output."
  )
}

cat("\nSample metadata:\n")
print(sample_meta)

write.csv(
  sample_meta,
  file.path(out_dir, "sample_metadata_from_LFQ_columns.csv"),
  row.names = FALSE
)

condition_order <- c(
  "Acute RDH12 Veh",
  "Acute RDH12 100 µM",
  "Acute RDH12 200 µM",
  "Recovery RDH12 Veh",
  "Recovery RDH12 100 µM",
  "Recovery RDH12 200 µM",
  "Recovery WT Veh",
  "Recovery WT 100 µM",
  "Recovery WT 200 µM"
)

condition_order <- condition_order[condition_order %in% sample_meta$Condition]

# ---------------------------------------------------------
# 5. LONG FORMAT TABLE
# ---------------------------------------------------------

long_df <- raw_df %>%
  mutate(
    ProteinID = .data[[protein_id_col]],
    Gene = if (!is.na(gene_col)) .data[[gene_col]] else .data[[protein_id_col]]
  ) %>%
  select(
    ProteinID,
    Gene,
    all_of(lfq_cols)
  ) %>%
  pivot_longer(
    cols = all_of(lfq_cols),
    names_to = "SampleCol",
    values_to = "LFQ"
  ) %>%
  left_join(sample_meta, by = "SampleCol") %>%
  mutate(
    LFQ = as.numeric(LFQ),
    Detected = !is.na(LFQ) & LFQ > 0,
    Condition = factor(Condition, levels = condition_order),
    Phase = factor(Phase, levels = c("Acute", "Recovery")),
    CellType = factor(CellType, levels = c("RDH12", "WT")),
    Treatment = factor(Treatment, levels = c("Veh", "100 µM", "200 µM"))
  )

write.csv(
  long_df,
  file.path(out_dir, "RawData_long_LFQ_table.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 6. DETECTION TABLE FOR UPSET PLOT
# ---------------------------------------------------------

detected_by_condition <- long_df %>%
  group_by(ProteinID, Gene, Condition) %>%
  summarise(
    n_reps = n(),
    n_detected = sum(Detected, na.rm = TRUE),
    detected_condition = n_detected >= min_reps_detected,
    .groups = "drop"
  )

detection_wide <- detected_by_condition %>%
  select(ProteinID, Gene, Condition, detected_condition) %>%
  pivot_wider(
    names_from = Condition,
    values_from = detected_condition,
    values_fill = FALSE
  )

write.csv(
  detection_wide,
  file.path(out_dir, "protein_detection_by_condition_for_upset.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 7. UPSET PLOT
# ---------------------------------------------------------

# Uses ComplexUpset if available. If not installed, install from CRAN.
if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  install.packages("ComplexUpset")
}

library(ComplexUpset)

upset_conditions <- condition_order

upset_df <- detection_wide %>%
  mutate(across(all_of(upset_conditions), as.logical))
p_upset_all <- ComplexUpset::upset(
  upset_df,
  intersect = upset_conditions,
  name = "Detected proteins",
  min_size = 1,
  n_intersections = 10,
  width_ratio = 0.18,
  
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      counts = TRUE
    )
  ),
  set_sizes = ComplexUpset::upset_set_size()
) +
  ggtitle(
    paste0(
      "Protein detection overlap across conditions\n"
      
    )
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_upset_all)

ggsave(
  file.path(out_dir, "UpSet_all_conditions_detected_proteins.pdf"),
  p_upset_all,
  width = 13,
  height = 7
)

ggsave(
  file.path(out_dir, "UpSet_all_conditions_detected_proteins.png"),
  p_upset_all,
  width = 13,
  height = 7,
  dpi = 300
)

# ---------------------------------------------------------
# 8. OPTIONAL SMALLER UPSET PLOTS
# ---------------------------------------------------------
# These are often easier to read than one large 9-condition plot.

acute_conditions <- condition_order[str_detect(condition_order, "^Acute")]
recovery_rdh12_conditions <- condition_order[str_detect(condition_order, "^Recovery RDH12")]
recovery_wt_conditions <- condition_order[str_detect(condition_order, "^Recovery WT")]

make_upset_subset <- function(set_cols, plot_title, file_stub) {
  
  if (length(set_cols) < 2) {
    message("Skipping ", file_stub, ": fewer than 2 sets.")
    return(NULL)
  }
  
  p <- ComplexUpset::upset(
    upset_df,
    intersect = set_cols,
    name = "Detected proteins",
    min_size = 1,
    width_ratio = 0.18,
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(
        counts = TRUE
      )
    ),
    set_sizes = ComplexUpset::upset_set_size()
  ) +
    ggtitle(plot_title) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  ggsave(
    file.path(out_dir, paste0(file_stub, ".pdf")),
    p,
    width = 9,
    height = 5.5
  )
  
  ggsave(
    file.path(out_dir, paste0(file_stub, ".png")),
    p,
    width = 9,
    height = 5.5,
    dpi = 300
  )
  
  return(p)
}

p_upset_acute <- make_upset_subset(
  acute_conditions,
  paste0("Acute RDH12 protein detection overlap\nDetected in ≥", min_reps_detected, " replicates"),
  "UpSet_acute_RDH12_detected_proteins"
)

p_upset_rec_rdh12 <- make_upset_subset(
  recovery_rdh12_conditions,
  paste0("Recovery RDH12 protein detection overlap\nDetected in ≥", min_reps_detected, " replicates"),
  "UpSet_recovery_RDH12_detected_proteins"
)

p_upset_rec_wt <- make_upset_subset(
  recovery_wt_conditions,
  paste0("Recovery WT protein detection overlap\nDetected in ≥", min_reps_detected, " replicates"),
  "UpSet_recovery_WT_detected_proteins"
)

# ---------------------------------------------------------
# 9. COEFFICIENT OF VARIATION TABLE
# ---------------------------------------------------------

cv_table <- long_df %>%
  group_by(ProteinID, Gene, Phase, CellType, Treatment, Condition) %>%
  summarise(
    n_total = n(),
    n_detected = sum(Detected, na.rm = TRUE),
    mean_LFQ = mean(LFQ[Detected], na.rm = TRUE),
    sd_LFQ = sd(LFQ[Detected], na.rm = TRUE),
    CV_percent = 100 * sd_LFQ / mean_LFQ,
    .groups = "drop"
  ) %>%
  mutate(
    CV_percent = if_else(
      n_detected >= min_reps_for_cv &
        is.finite(CV_percent) &
        mean_LFQ > 0,
      CV_percent,
      NA_real_
    ),
    CV_percent_capped = pmin(CV_percent, cv_plot_cap),
    Condition = factor(as.character(Condition), levels = condition_order),
    Phase = factor(as.character(Phase), levels = c("Acute", "Recovery")),
    CellType = factor(as.character(CellType), levels = c("RDH12", "WT")),
    Treatment = factor(as.character(Treatment), levels = c("Veh", "100 µM", "200 µM"))
  )

write.csv(
  cv_table,
  file.path(out_dir, "protein_CV_by_condition.csv"),
  row.names = FALSE
)

cv_summary <- cv_table %>%
  group_by(Condition) %>%
  summarise(
    n_proteins_with_CV = sum(!is.na(CV_percent)),
    median_CV = median(CV_percent, na.rm = TRUE),
    mean_CV = mean(CV_percent, na.rm = TRUE),
    q25_CV = quantile(CV_percent, 0.25, na.rm = TRUE),
    q75_CV = quantile(CV_percent, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  cv_summary,
  file.path(out_dir, "protein_CV_summary_by_condition.csv"),
  row.names = FALSE
)

cat("\nCV summary:\n")
print(cv_summary)

# ---------------------------------------------------------
# 10. CV BOXPLOT
# ---------------------------------------------------------

p_cv_box <- ggplot(
  cv_table %>% filter(!is.na(CV_percent)),
  aes(x = Condition, y = CV_percent_capped)
) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(
    aes(color = Phase),
    width = 0.15,
    alpha = 0.18,
    size = 0.5
  ) +
  geom_hline(
    yintercept = c(20),
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  theme_bw(base_size=16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Protein-level coefficient of variation by condition",
    subtitle = paste0(
      "CV calculated across biological replicates; dashed line at 20% "
    ),
    x = NULL,
    y = "CV (%)",
    color = "Phase"
  )

print(p_cv_box)

ggsave(
  file.path(out_dir, "CV_boxplot_by_condition.pdf"),
  p_cv_box,
  width = 11,
  height = 6
)

ggsave(
  file.path(out_dir, "CV_boxplot_by_condition.png"),
  p_cv_box,
  width = 11,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------
# 11. CV VIOLIN + BOXPLOT
# ---------------------------------------------------------

p_cv_violin <- ggplot(
  cv_table %>% filter(!is.na(CV_percent)),
  aes(x = Condition, y = CV_percent_capped, fill = Phase)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.7,
    color = "grey30"
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_hline(
    yintercept = c(20),
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  theme_bw(base_size=16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Distribution of protein CVs by condition",
    subtitle = paste0(
      "CV calculated across biological replicates; dashed line at 20% "
   
    ),
    x = NULL,
    y = "CV (%)",
    fill = "Phase"
  )

print(p_cv_violin)

ggsave(
  file.path(out_dir, "CV_violin_by_condition.pdf"),
  p_cv_violin,
  width = 11,
  height = 6
)

ggsave(
  file.path(out_dir, "CV_violin_by_condition.png"),
  p_cv_violin,
  width = 11,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------
# 12. CV DENSITY PLOT
# ---------------------------------------------------------

p_cv_density <- ggplot(
  cv_table %>% filter(!is.na(CV_percent), CV_percent <= cv_plot_cap),
  aes(x = CV_percent, color = Condition)
) +
  geom_density(linewidth = 0.8) +
  geom_vline(
    xintercept = c(20, 30),
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Protein CV density by condition",
    x = "CV (%)",
    y = "Density",
    color = "Condition"
  )

print(p_cv_density)

ggsave(
  file.path(out_dir, "CV_density_by_condition.pdf"),
  p_cv_density,
  width = 10,
  height = 6
)

ggsave(
  file.path(out_dir, "CV_density_by_condition.png"),
  p_cv_density,
  width = 10,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------
# 13. DETECTION COUNTS BY CONDITION
# ---------------------------------------------------------

detection_summary <- detected_by_condition %>%
  group_by(Condition) %>%
  summarise(
    n_proteins_detected = sum(detected_condition, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(as.character(Condition), levels = condition_order)
  )

write.csv(
  detection_summary,
  file.path(out_dir, "protein_detection_counts_by_condition.csv"),
  row.names = FALSE
)

p_detection_counts <- ggplot(
  detection_summary,
  aes(x = Condition, y = n_proteins_detected)
) +
  geom_col(width = 0.7) +
  theme_bw(base_size=16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste0(
      "Number of detected proteins by condition\n"
    ),
    x = NULL,
    y = "Detected proteins"
  )

print(p_detection_counts)

ggsave(
  file.path(out_dir, "Detected_protein_counts_by_condition.pdf"),
  p_detection_counts,
  width = 10,
  height = 5.5
)

ggsave(
  file.path(out_dir, "Detected_protein_counts_by_condition.png"),
  p_detection_counts,
  width = 10,
  height = 5.5,
  dpi = 300
)

# ---------------------------------------------------------
# 14. FINAL MESSAGE
# ---------------------------------------------------------

cat("\nSaved output directory:\n")
cat(out_dir, "\n")

cat("\nMain output files:\n")
cat(file.path(out_dir, "UpSet_all_conditions_detected_proteins.pdf"), "\n")
cat(file.path(out_dir, "CV_boxplot_by_condition.pdf"), "\n")
cat(file.path(out_dir, "CV_violin_by_condition.pdf"), "\n")
cat(file.path(out_dir, "CV_density_by_condition.pdf"), "\n")
cat(file.path(out_dir, "Detected_protein_counts_by_condition.pdf"), "\n")
cat(file.path(out_dir, "protein_CV_by_condition.csv"), "\n")
cat(file.path(out_dir, "protein_CV_summary_by_condition.csv"), "\n")
cat(file.path(out_dir, "protein_detection_by_condition_for_upset.csv"), "\n")