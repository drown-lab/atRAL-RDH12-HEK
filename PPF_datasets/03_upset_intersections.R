# =========================================================
# UpSet plot for PPF protein detection overlap
# From PPFtable_v4 created by 01_Prep_data.R
#
# Input:
#   PPFtable_v4
#
# Outputs:
#   1. UpSet plot of protein detection overlap across PPF samples
#   2. CSV tables used for plotting
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

input_dir <- if (dir.exists("PPF_datasets")) "PPF_datasets" else "."
prep_script_path <- file.path(input_dir, "01_Prep_data.R")
out_dir <- file.path(input_dir, "output_upset")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Detection rule for UpSet plot
# A protein is considered detected in a sample if PG.MaxLFQ is above this value.
min_intensity_for_detected <- 0

# Number of intersections to show in the UpSet plot.
top_n_intersections <- 9


# ---------------------------------------------------------
# 2. READ PREPARED PPF DATA
# ---------------------------------------------------------

if (!file.exists(prep_script_path)) {
  prep_script_path <- "01_Prep_data.R"
}

if (!file.exists(prep_script_path)) {
  stop("Could not find 01_Prep_data.R. Run this script from the project root or PPF_datasets folder.")
}

source(prep_script_path)

if (!exists("PPFtable_v4")) {
  stop("PPFtable_v4 was not created by 01_Prep_data.R.")
}

required_cols <- c("sample_name", "Run", "Protein.Group", "Genes", "PG.MaxLFQ")
missing_cols <- setdiff(required_cols, colnames(PPFtable_v4))

if (length(missing_cols) > 0) {
  stop("PPFtable_v4 is missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!is.numeric(PPFtable_v4$PG.MaxLFQ)) {
  stop("PG.MaxLFQ must be numeric before calculating sample intersections.")
}

cat("\nPPFtable_v4 dimensions:\n")
print(dim(PPFtable_v4))

cat("\nPPFtable_v4 columns:\n")
print(colnames(PPFtable_v4))


# ---------------------------------------------------------
# 3. BUILD SAMPLE METADATA
# ---------------------------------------------------------

sample_meta <- PPFtable_v4 %>%
  transmute(
    SampleName = as.character(coalesce(sample_name, Run))
  ) %>%
  distinct() %>%
  mutate(
    Condition = str_remove(SampleName, "\\d+$"),
    Replicate = str_extract(SampleName, "\\d+$"),
    Condition = if_else(is.na(Condition) | Condition == "", SampleName, Condition),
    Replicate = if_else(is.na(Replicate) | Replicate == "", SampleName, Replicate)
  ) %>%
  arrange(Condition, Replicate, SampleName)

cat("\nSample metadata:\n")
print(sample_meta)

write.csv(
  sample_meta,
  file.path(out_dir, "sample_metadata_for_upset.csv"),
  row.names = FALSE
)

sample_order <- sample_meta$SampleName


# ---------------------------------------------------------
# 4. LONG FORMAT DETECTION TABLE
# ---------------------------------------------------------

long_df <- PPFtable_v4 %>%
  mutate(
    SampleName = as.character(coalesce(sample_name, Run)),
    ProteinID = Protein.Group,
    Gene = if_else(is.na(Genes) | Genes == "", Protein.Group, Genes),
    PG.MaxLFQ = as.numeric(PG.MaxLFQ)
  ) %>%
  group_by(ProteinID, Gene, SampleName) %>%
  summarise(
    PG.MaxLFQ = if_else(all(is.na(PG.MaxLFQ)), NA_real_, max(PG.MaxLFQ, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(sample_meta, by = "SampleName") %>%
  mutate(
    Detected = !is.na(PG.MaxLFQ) & PG.MaxLFQ > min_intensity_for_detected,
    SampleName = factor(SampleName, levels = sample_order),
    Condition = factor(Condition, levels = unique(sample_meta$Condition))
  )

write.csv(
  long_df,
  file.path(out_dir, "protein_detection_long_table.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# 5. DETECTION TABLE FOR UPSET PLOT
# ---------------------------------------------------------

detection_wide <- long_df %>%
  select(ProteinID, Gene, SampleName, Detected) %>%
  pivot_wider(
    names_from = SampleName,
    values_from = Detected,
    values_fill = FALSE
  )

write.csv(
  detection_wide,
  file.path(out_dir, "protein_detection_by_sample_for_upset.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# 6. UPSET PLOT
# ---------------------------------------------------------

# Uses ComplexUpset if available. If not installed, install from CRAN.
if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  install.packages("ComplexUpset")
}

library(ComplexUpset)

upset_df <- detection_wide %>%
  mutate(across(all_of(sample_order), as.logical))

p_upset_all <- ComplexUpset::upset(
  upset_df,
  intersect = sample_order,
  name = "Detected proteins",
  min_size = 1,
  n_intersections = top_n_intersections,
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
      "intersection"
    )
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_upset_all)

ggsave(
  file.path(out_dir, "UpSet_all_samples_detected_proteins.pdf"),
  p_upset_all,
  width = 8,
  height = 7
)

ggsave(
  file.path(out_dir, "UpSet_all_samples_detected_proteins.png"),
  p_upset_all,
  width = 8,
  height = 7,
  dpi = 300
)


# ---------------------------------------------------------
# 7. FINAL MESSAGE
# ---------------------------------------------------------

cat("\nSaved output directory:\n")
cat(out_dir, "\n")

cat("\nMain output files:\n")
cat(file.path(out_dir, "UpSet_all_samples_detected_proteins.pdf"), "\n")
cat(file.path(out_dir, "UpSet_all_samples_detected_proteins.png"), "\n")
cat(file.path(out_dir, "protein_detection_by_sample_for_upset.csv"), "\n")
cat(file.path(out_dir, "protein_detection_long_table.csv"), "\n")
cat(file.path(out_dir, "sample_metadata_for_upset.csv"), "\n")

cat("\nAnalyzed proteins:\n")
cat(dplyr::n_distinct(long_df$ProteinID), "\n")

cat("\nAnalyzed samples:\n")
cat(length(sample_order), "\n")
