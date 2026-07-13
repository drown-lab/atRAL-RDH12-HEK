# =========================================================
# Lipidomics CV% analysis for acute and recovery datasets
#
# Fixes:
#   - robustly detects RDH12 and WT/GFP/Control genotype labels
#   - collapses lipid classes to main classes only: PC, PE, PI, DG, TG, etc.
#   - removes ether/plasmalogen subclasses from CV class plot
#   - uses only one dashed CV reference line at 30%
#
# CV (%) = standard deviation / mean * 100
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(forcats)
})

# ---------------------------------------------------------
# 1. INPUT FILES
# ---------------------------------------------------------

acute_lipid_file <- "Lipidomics/output_txts/Acute_filteredlipidtable.csv"
recovery_lipid_file <- "Lipidomics/output_txts/Recovery_filteredlipidtable2.csv"

acute_sample_file <- "Lipidomics/metadata/Sample_description guide_5hratral.csv"
recovery_sample_file <- "Lipidomics/metadata/Sample_description guide_recovery.csv"

out_dir <- "Lipidomics/Figures/QC_CV/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_reps_for_cv <- 2
cv_plot_cap <- 200
cv_reference_line <- 30

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------

standardize_genotype <- function(x) {
  x_clean <- x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_to_upper()
  
  case_when(
    str_detect(x_clean, "RDH12") ~ "RDH12",
    str_detect(x_clean, "RHD12") ~ "RDH12",  # catches common typo
    str_detect(x_clean, "GFP") ~ "WT",
    str_detect(x_clean, "WT") ~ "WT",
    str_detect(x_clean, "CONTROL") ~ "WT",
    TRUE ~ x_clean
  )
}

standardize_treatment <- function(x) {
  x_clean <- x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_to_upper()
  
  case_when(
    x_clean %in% c("VEH", "CONTROL", "CTRL", "0") ~ "Veh",
    x_clean %in% c("100", "100UM", "100µM", "100ΜM") ~ "100 µM",
    x_clean %in% c("200", "200UM", "200µM", "200ΜM") ~ "200 µM",
    TRUE ~ as.character(x)
  )
}

standardize_sample_metadata <- function(sample_df, experiment_label) {
  
  colnames(sample_df) <- str_trim(colnames(sample_df))
  
  # Fix possible column name variants
  if ("Treament" %in% colnames(sample_df)) {
    sample_df <- sample_df %>% rename(Treatment = Treament)
  }
  
  if ("Treatment " %in% colnames(sample_df)) {
    sample_df <- sample_df %>% rename(Treatment = `Treatment `)
  }
  
  if ("GeneType" %in% colnames(sample_df)) {
    sample_df <- sample_df %>% rename(Genotype = GeneType)
  }
  
  if ("Genetype" %in% colnames(sample_df)) {
    sample_df <- sample_df %>% rename(Genotype = Genetype)
  }
  
  if ("genoType" %in% colnames(sample_df)) {
    sample_df <- sample_df %>% rename(Genotype = genoType)
  }
  
  if ("Sample" %in% colnames(sample_df)) {
    sample_df <- sample_df %>% rename(sample = Sample)
  }
  
  required_meta_cols <- c("sample", "Treatment", "Rep", "Genotype")
  missing_meta_cols <- setdiff(required_meta_cols, colnames(sample_df))
  
  if (length(missing_meta_cols) > 0) {
    stop(
      paste0(
        "Missing required sample metadata columns in ", experiment_label, " sample guide:\n",
        paste(missing_meta_cols, collapse = "\n")
      )
    )
  }
  
  sample_df %>%
    mutate(
      sample = as.character(sample),
      Treatment_raw = as.character(Treatment),
      Genotype_raw = as.character(Genotype),
      Rep = as.character(Rep),
      
      Genotype = standardize_genotype(Genotype_raw),
      Treatment = standardize_treatment(Treatment_raw),
      
      Experiment = experiment_label,
      Condition = paste(Experiment, Genotype, Treatment, sep = " | ")
    )
}

clean_lipid_table <- function(lipid_df) {
  
  bad_name_cols <- which(is.na(colnames(lipid_df)) | colnames(lipid_df) == "")
  
  if (length(bad_name_cols) > 0) {
    lipid_df <- lipid_df[, -bad_name_cols, drop = FALSE]
  }
  
  lipid_df <- lipid_df %>%
    select(-any_of(c("Unnamed: 0", "...1", "X", "X.1")))
  
  colnames(lipid_df) <- make.unique(colnames(lipid_df))
  
  lipid_df
}

get_main_lipid_class_from_name <- function(x) {
  x <- as.character(x) %>%
    str_trim()
  
  case_when(
    str_detect(x, regex("^PC\\(|^PC\\s|^PCO-|^PC O-|^PCP-|^PC P-|Aze|COOH|Azelaoyl", ignore_case = TRUE)) ~ "PC",
    str_detect(x, regex("^LPC\\s|^LPC\\(", ignore_case = TRUE)) ~ "LPC",
    
    str_detect(x, regex("^PE\\(|^PE\\s|^PEO-|^PE O-|^PEP-|^PE P-", ignore_case = TRUE)) ~ "PE",
    str_detect(x, regex("^LPE\\s|^LPE\\(", ignore_case = TRUE)) ~ "LPE",
    
    str_detect(x, regex("^PI\\(|^PI\\s|^PIO-|^PI O-", ignore_case = TRUE)) ~ "PI",
    str_detect(x, regex("^PS\\(|^PS\\s|^PSO-|^PS O-", ignore_case = TRUE)) ~ "PS",
    str_detect(x, regex("^PG\\(|^PG\\s", ignore_case = TRUE)) ~ "PG",
    
    str_detect(x, regex("^DG\\(|^DG\\s|^DGO-|^DG O-|^\\[?DG", ignore_case = TRUE)) ~ "DG",
    str_detect(x, regex("^TG\\(|^TG\\s|^\\[?TG", ignore_case = TRUE)) ~ "TG",
    
    str_detect(x, regex("^Cer\\(|^Cer\\s", ignore_case = TRUE)) ~ "Cer",
    str_detect(x, regex("^dhCer\\(|^dhCer\\s", ignore_case = TRUE)) ~ "dhCer",
    str_detect(x, regex("^SM\\(|^SM\\s", ignore_case = TRUE)) ~ "SM",
    
    str_detect(x, regex("^CE\\s|^\\[?CE", ignore_case = TRUE)) ~ "CE",
    str_detect(x, regex("^CAR\\s|^CAR\\(", ignore_case = TRUE)) ~ "CAR",
    
    TRUE ~ "Other"
  )
}

get_main_lipid_class <- function(df) {
  
  # Prefer explicit lipid name and classify from that.
  if ("lipid_name" %in% colnames(df)) {
    return(get_main_lipid_class_from_name(df$lipid_name))
  }
  
  if ("name" %in% colnames(df)) {
    return(get_main_lipid_class_from_name(df$name))
  }
  
  # Fallback to existing class columns, but collapse ether/plasmalogen labels.
  class_col <- NULL
  
  if ("lipid_class" %in% colnames(df)) {
    class_col <- df$lipid_class
  } else if ("lipid_class1" %in% colnames(df)) {
    class_col <- df$lipid_class1
  } else if ("lipid_class_ether" %in% colnames(df)) {
    class_col <- df$lipid_class_ether
  } else if ("lipid_class_ether2" %in% colnames(df)) {
    class_col <- df$lipid_class_ether2
  } else {
    return(rep("Other", nrow(df)))
  }
  
  class_col <- as.character(class_col) %>%
    str_trim() %>%
    str_replace("^NA_", "") %>%
    str_replace("^O-_", "") %>%
    str_replace("^P-_", "") %>%
    str_replace("_O-$", "") %>%
    str_replace("_P-$", "")
  
  case_when(
    str_detect(class_col, "PC") ~ "PC",
    str_detect(class_col, "LPC") ~ "LPC",
    str_detect(class_col, "PE") ~ "PE",
    str_detect(class_col, "LPE") ~ "LPE",
    str_detect(class_col, "PI") ~ "PI",
    str_detect(class_col, "PS") ~ "PS",
    str_detect(class_col, "PG") ~ "PG",
    str_detect(class_col, "DG") ~ "DG",
    str_detect(class_col, "TG") ~ "TG",
    str_detect(class_col, "Cer") ~ "Cer",
    str_detect(class_col, "dhCer") ~ "dhCer",
    str_detect(class_col, "SM") ~ "SM",
    str_detect(class_col, "CE") ~ "CE",
    str_detect(class_col, "CAR") ~ "CAR",
    TRUE ~ "Other"
  )
}

make_lipid_long <- function(lipid_df, sample_meta, experiment_label) {
  
  lipid_df <- clean_lipid_table(lipid_df)
  
  if (!"lipid_name" %in% colnames(lipid_df)) {
    stop(paste0("lipid_name column not found in ", experiment_label, " lipid table."))
  }
  
  sample_cols <- sample_meta$sample
  missing_sample_cols <- setdiff(sample_cols, colnames(lipid_df))
  
  if (length(missing_sample_cols) > 0) {
    stop(
      paste0(
        "These sample columns from the ", experiment_label, " sample guide are missing in the lipid table:\n",
        paste(missing_sample_cols, collapse = "\n")
      )
    )
  }
  
  lipid_df2 <- lipid_df %>%
    mutate(
      Lipid = lipid_name,
      LipidClass = get_main_lipid_class(lipid_df),
      LipidClass = if_else(is.na(LipidClass) | LipidClass == "", "Other", as.character(LipidClass))
    )
  
  long_df <- lipid_df2 %>%
    select(
      Lipid,
      LipidClass,
      all_of(sample_cols)
    ) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "Intensity"
    ) %>%
    left_join(sample_meta, by = "sample") %>%
    mutate(
      Intensity = suppressWarnings(as.numeric(Intensity)),
      Detected = !is.na(Intensity) & Intensity > 0,
      log2Intensity = if_else(Detected, log2(Intensity), NA_real_)
    )
  
  long_df
}

calculate_cv_table <- function(long_df) {
  
  long_df %>%
    group_by(
      Experiment,
      Genotype,
      Treatment,
      Condition,
      Lipid,
      LipidClass
    ) %>%
    summarise(
      n_total = n(),
      n_detected = sum(Detected, na.rm = TRUE),
      mean_intensity = mean(Intensity[Detected], na.rm = TRUE),
      sd_intensity = sd(Intensity[Detected], na.rm = TRUE),
      CV_percent = 100 * sd_intensity / mean_intensity,
      mean_log2Intensity = mean(log2Intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      CV_percent = if_else(
        n_detected >= min_reps_for_cv &
          is.finite(CV_percent) &
          mean_intensity > 0,
        CV_percent,
        NA_real_
      ),
      CV_percent_capped = pmin(CV_percent, cv_plot_cap)
    )
}

# ---------------------------------------------------------
# 3. READ INPUT DATA
# ---------------------------------------------------------

acute_lipid_df <- read.csv(
  acute_lipid_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

recovery_lipid_df <- read.csv(
  recovery_lipid_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

acute_sample_df <- read.csv(
  acute_sample_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

recovery_sample_df <- read.csv(
  recovery_sample_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

acute_sample_meta <- standardize_sample_metadata(
  sample_df = acute_sample_df,
  experiment_label = "Acute"
)

recovery_sample_meta <- standardize_sample_metadata(
  sample_df = recovery_sample_df,
  experiment_label = "Recovery"
)

sample_meta_all <- bind_rows(
  acute_sample_meta,
  recovery_sample_meta
)

cat("\nCleaned sample metadata:\n")
print(sample_meta_all)

cat("\nGenotype counts after cleaning:\n")
print(sample_meta_all %>% count(Experiment, Genotype, Treatment))

write.csv(
  sample_meta_all,
  file.path(out_dir, "lipid_sample_metadata_cleaned.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 4. MAKE LONG TABLES
# ---------------------------------------------------------

acute_long <- make_lipid_long(
  lipid_df = acute_lipid_df,
  sample_meta = acute_sample_meta,
  experiment_label = "Acute"
)

recovery_long <- make_lipid_long(
  lipid_df = recovery_lipid_df,
  sample_meta = recovery_sample_meta,
  experiment_label = "Recovery"
)

lipid_long_all <- bind_rows(
  acute_long,
  recovery_long
)

condition_order <- c(
  "Acute | RDH12 | Veh",
  "Acute | RDH12 | 100 µM",
  "Acute | RDH12 | 200 µM",
  "Recovery | RDH12 | Veh",
  "Recovery | RDH12 | 100 µM",
  "Recovery | RDH12 | 200 µM",
  "Recovery | WT | Veh",
  "Recovery | WT | 100 µM",
  "Recovery | WT | 200 µM"
)

condition_order <- condition_order[condition_order %in% unique(lipid_long_all$Condition)]

lipid_long_all <- lipid_long_all %>%
  mutate(
    Condition = factor(Condition, levels = condition_order),
    Experiment = factor(Experiment, levels = c("Acute", "Recovery")),
    Genotype = factor(Genotype, levels = c("RDH12", "WT")),
    Treatment = factor(Treatment, levels = c("Veh", "100 µM", "200 µM")),
    LipidClass = factor(LipidClass)
  )

write.csv(
  lipid_long_all,
  file.path(out_dir, "lipid_long_raw_intensity_table.csv"),
  row.names = FALSE
)

cat("\nCondition counts after long-format conversion:\n")
print(lipid_long_all %>% distinct(sample, Experiment, Genotype, Treatment, Condition) %>% count(Experiment, Genotype, Treatment, Condition))

cat("\nMain lipid class counts:\n")
print(lipid_long_all %>% distinct(Lipid, LipidClass, Experiment) %>% count(Experiment, LipidClass, sort = TRUE))

# ---------------------------------------------------------
# 5. CALCULATE CV TABLE
# ---------------------------------------------------------

cv_table <- calculate_cv_table(lipid_long_all) %>%
  mutate(
    Condition = factor(as.character(Condition), levels = condition_order),
    Experiment = factor(as.character(Experiment), levels = c("Acute", "Recovery")),
    Genotype = factor(as.character(Genotype), levels = c("RDH12", "WT")),
    Treatment = factor(as.character(Treatment), levels = c("Veh", "100 µM", "200 µM")),
    LipidClass = factor(LipidClass)
  )

write.csv(
  cv_table,
  file.path(out_dir, "lipid_CV_by_condition.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 6. SUMMARY TABLES
# ---------------------------------------------------------

cv_summary_condition <- cv_table %>%
  group_by(Experiment, Genotype, Treatment, Condition) %>%
  summarise(
    n_lipids_with_CV = sum(!is.na(CV_percent)),
    median_CV = median(CV_percent, na.rm = TRUE),
    mean_CV = mean(CV_percent, na.rm = TRUE),
    q25_CV = quantile(CV_percent, 0.25, na.rm = TRUE),
    q75_CV = quantile(CV_percent, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  cv_summary_condition,
  file.path(out_dir, "lipid_CV_summary_by_condition.csv"),
  row.names = FALSE
)

cv_summary_class <- cv_table %>%
  group_by(Experiment, Genotype, Treatment, Condition, LipidClass) %>%
  summarise(
    n_lipids_with_CV = sum(!is.na(CV_percent)),
    median_CV = median(CV_percent, na.rm = TRUE),
    mean_CV = mean(CV_percent, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  cv_summary_class,
  file.path(out_dir, "lipid_CV_summary_by_class.csv"),
  row.names = FALSE
)

detection_summary <- cv_table %>%
  group_by(Experiment, Genotype, Treatment, Condition) %>%
  summarise(
    n_lipids_total = n(),
    n_lipids_detected_min2 = sum(n_detected >= min_reps_for_cv, na.rm = TRUE),
    n_lipids_detected_all3 = sum(n_detected == 3, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  detection_summary,
  file.path(out_dir, "lipid_detection_summary_by_condition.csv"),
  row.names = FALSE
)

cat("\nCV summary by condition:\n")
print(cv_summary_condition)

cat("\nDetection summary by condition:\n")
print(detection_summary)

# ---------------------------------------------------------
# 7. PLOT: CV VIOLIN + BOXPLOT BY CONDITION
# ---------------------------------------------------------

p_cv_violin <- ggplot(
  cv_table %>% filter(!is.na(CV_percent)),
  aes(x = Condition, y = CV_percent_capped, fill = Experiment)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.7,
    color = "grey30"
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.9,
    color = "grey20"
  ) +
  geom_hline(
    yintercept = cv_reference_line,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "Acute" = "#E76F51",
      "Recovery" = "#457B9D"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Lipid coefficient of variation by condition",
    subtitle = "CV calculated from raw intensities across biological replicates; dashed line at 30% CV",
    x = NULL,
    y = "CV (%)",
    fill = "Experiment"
  )

print(p_cv_violin)

ggsave(
  file.path(out_dir, "Lipid_CV_violin_by_condition.pdf"),
  p_cv_violin,
  width = 11,
  height = 6
)

ggsave(
  file.path(out_dir, "Lipid_CV_violin_by_condition.png"),
  p_cv_violin,
  width = 11,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------
# 8. PLOT: CV BOXPLOT BY CONDITION
# ---------------------------------------------------------

p_cv_box <- ggplot(
  cv_table %>% filter(!is.na(CV_percent)),
  aes(x = Condition, y = CV_percent_capped, color = Experiment)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 0.6
  ) +
  geom_hline(
    yintercept = cv_reference_line,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  scale_color_manual(
    values = c(
      "Acute" = "#E76F51",
      "Recovery" = "#457B9D"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Lipid CV distribution by condition",
    subtitle = "CV calculated from raw intensities across biological replicates; dashed line at 30% CV",
    x = NULL,
    y = "CV (%)",
    color = "Experiment"
  )

print(p_cv_box)

ggsave(
  file.path(out_dir, "Lipid_CV_boxplot_by_condition.pdf"),
  p_cv_box,
  width = 11,
  height = 6
)

ggsave(
  file.path(out_dir, "Lipid_CV_boxplot_by_condition.png"),
  p_cv_box,
  width = 11,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------
# 9. PLOT: MEDIAN CV BY CONDITION
# ---------------------------------------------------------

p_median_cv <- ggplot(
  cv_summary_condition,
  aes(x = Condition, y = median_CV, fill = Experiment)
) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = q25_CV, ymax = q75_CV),
    width = 0.2,
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = cv_reference_line,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "Acute" = "#E76F51",
      "Recovery" = "#457B9D"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Median lipid CV by condition",
    subtitle = "Bars show median CV; error bars show interquartile range; dashed line at 30% CV",
    x = NULL,
    y = "Median CV (%)",
    fill = "Experiment"
  )

print(p_median_cv)

ggsave(
  file.path(out_dir, "Lipid_median_CV_by_condition.pdf"),
  p_median_cv,
  width = 11,
  height = 5.5
)

ggsave(
  file.path(out_dir, "Lipid_median_CV_by_condition.png"),
  p_median_cv,
  width = 11,
  height = 5.5,
  dpi = 300
)

# ---------------------------------------------------------
# 10. PLOT: CV BY MAIN LIPID CLASS
# ---------------------------------------------------------

class_keep <- cv_table %>%
  filter(!is.na(CV_percent)) %>%
  count(LipidClass, name = "n") %>%
  filter(n >= 5, !is.na(LipidClass), LipidClass != "Other") %>%
  arrange(desc(n)) %>%
  pull(LipidClass)

p_cv_class <- ggplot(
  cv_table %>%
    filter(!is.na(CV_percent), LipidClass %in% class_keep) %>%
    mutate(LipidClass = fct_infreq(LipidClass)),
  aes(x = LipidClass, y = CV_percent_capped, fill = Experiment)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = cv_reference_line,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "Acute" = "#E76F51",
      "Recovery" = "#457B9D"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Lipid CV by main lipid class",
    subtitle = "Showing main lipid classes with ≥5 lipids with CV values; dashed line at 30% CV",
    x = NULL,
    y = "CV (%)",
    fill = "Experiment"
  )

print(p_cv_class)

ggsave(
  file.path(out_dir, "Lipid_CV_by_main_class.pdf"),
  p_cv_class,
  width = 10,
  height = 6
)

ggsave(
  file.path(out_dir, "Lipid_CV_by_main_class.png"),
  p_cv_class,
  width = 10,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------
# 11. PLOT: DETECTED LIPID COUNTS BY CONDITION
# ---------------------------------------------------------

p_detection <- ggplot(
  detection_summary,
  aes(x = Condition, y = n_lipids_detected_min2, fill = Experiment)
) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Acute" = "#E76F51",
      "Recovery" = "#457B9D"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste0(
      "Detected lipid count by condition\n",
      "Detected in ≥", min_reps_for_cv, " biological replicates"
    ),
    x = NULL,
    y = "Detected lipids",
    fill = "Experiment"
  )

print(p_detection)

ggsave(
  file.path(out_dir, "Lipid_detected_counts_by_condition.pdf"),
  p_detection,
  width = 10,
  height = 5.5
)

ggsave(
  file.path(out_dir, "Lipid_detected_counts_by_condition.png"),
  p_detection,
  width = 10,
  height = 5.5,
  dpi = 300
)

# ---------------------------------------------------------
# 12. FINAL MESSAGE
# ---------------------------------------------------------

cat("\nSaved output directory:\n")
cat(out_dir, "\n")

cat("\nMain output files:\n")
cat(file.path(out_dir, "lipid_CV_by_condition.csv"), "\n")
cat(file.path(out_dir, "lipid_CV_summary_by_condition.csv"), "\n")
cat(file.path(out_dir, "lipid_CV_summary_by_class.csv"), "\n")
cat(file.path(out_dir, "lipid_detection_summary_by_condition.csv"), "\n")
cat(file.path(out_dir, "Lipid_CV_violin_by_condition.pdf"), "\n")
cat(file.path(out_dir, "Lipid_CV_boxplot_by_condition.pdf"), "\n")
cat(file.path(out_dir, "Lipid_median_CV_by_condition.pdf"), "\n")
cat(file.path(out_dir, "Lipid_CV_by_main_class.pdf"), "\n")
cat(file.path(out_dir, "Lipid_detected_counts_by_condition.pdf"), "\n")