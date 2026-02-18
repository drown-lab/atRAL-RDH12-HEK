library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# =========================================================
# Choose the SE object to summarize (pre-imputation!)
# =========================================================
se <- data_filt   # or data_se

mat  <- assay(se)
meta <- as.data.frame(colData(se))

stopifnot("condition" %in% names(meta))
stopifnot(all(colnames(mat) %in% rownames(meta)) || nrow(meta) == ncol(mat))

# Ensure we can join by sample name
meta <- meta |> tibble::rownames_to_column("sample")

# =========================================================
# (1) Categorize proteins per CONDITION by replicate detection
#     MNAR = 0/3
#     MAR_1of3 = 1/3
#     MAR_2of3 = 2/3
#     Present_all = 3/3
# =========================================================

det_by_cond <- as.data.frame(!is.na(mat)) |>
  tibble::rownames_to_column("protein") |>
  pivot_longer(-protein, names_to = "sample", values_to = "detected") |>
  left_join(meta |> select(sample, condition), by = "sample") |>
  group_by(protein, condition) |>
  summarise(
    n_detect = sum(detected),
    n_reps   = n(),
    .groups  = "drop"
  ) |>
  mutate(
    detect_class = case_when(
      n_detect == 0 ~ "MNAR_0of3",
      n_detect == 1 ~ "MAR_1of3",
      n_detect == 2 ~ "MAR_2of3",
      n_detect >= 3 ~ "Present_3of3",
      TRUE ~ paste0("Other_", n_detect, "of", n_reps)
    ),
    detect_class = factor(
      detect_class,
      levels = c("MNAR_0of3", "MAR_1of3", "MAR_2of3", "Present_3of3")
    )
  )

# =========================================================
# (1.1) Summaries as NEW TABLES
#   A) counts per condition × class
#   B) totals across dataset
#   C) how many conditions per protein are 0/1/2/3 (wide summary)
# =========================================================
library(dplyr)
# A) Condition-level counts
summary_by_condition <- det_by_cond |>
  dplyr::count(condition, detect_class, name = "n_proteins") |>
  tidyr::pivot_wider(
    names_from  = detect_class,
    values_from = n_proteins,
    values_fill = 0
  ) |>
  dplyr::arrange(condition)


# B) Overall totals across all protein×condition blocks
summary_overall_blocks <- det_by_cond |>
  dplyr::count(detect_class, name = "n_protein_condition_blocks") |>
  arrange(detect_class)|>
  unique()

# C) Protein-level: for each protein, how many conditions are 0/1/2/3
summary_per_protein <- det_by_cond |>
  mutate(detect_bucket = as.character(detect_class)) |>
  group_by(protein) |>
  summarise(
    n_conditions_total = n(),
    n_cond_MNAR_0of3   = sum(detect_bucket == "MNAR_0of3"),
    n_cond_MAR_1of3    = sum(detect_bucket == "MAR_1of3"),
    n_cond_MAR_2of3    = sum(detect_bucket == "MAR_2of3"),
    n_cond_Present_3of3= sum(detect_bucket == "Present_3of3"),
    .groups = "drop"
  )

# And a distribution table: how many proteins have >=k conditions with >=2/3, etc.
dist_conditions_ge2 <- summary_per_protein |>
  dplyr::count(n_cond_MAR_2of3, name = "n_proteins") |>
  arrange(n_cond_MAR_2of3)

# Print tables
summary_by_condition
summary_overall_blocks
summary_per_protein |> head(10)
dist_conditions_ge2

# If you want to save to CSV:
# write.csv(summary_by_condition, "summary_by_condition.csv", row.names = FALSE)
# write.csv(summary_overall_blocks, "summary_overall_blocks.csv", row.names = FALSE)
# write.csv(summary_per_protein, "summary_per_protein.csv", row.names = FALSE)

# =========================================================
# (1.2) Missing value plot (editable ggplot heatmap)
#     We'll plot proteins (rows) vs samples (columns) as detected/missing.
#     For readability, show the top N proteins with the most missingness.
# =========================================================

# Build long missingness table protein×sample
miss_long <- as.data.frame(is.na(mat)) |>
  tibble::rownames_to_column("protein") |>
  pivot_longer(-protein, names_to = "sample", values_to = "is_missing") |>
  left_join(meta |> select(sample, condition), by = "sample")

# Choose how many proteins to display (heatmap gets huge otherwise)
top_n <- 300

# pick proteins with highest missing fraction
prot_miss_rank <- miss_long |>
  group_by(protein) |>
  summarise(miss_frac = mean(is_missing), .groups = "drop") |>
  arrange(desc(miss_frac)) |>
  slice_head(n = top_n)

miss_plot_df <- miss_long |>
  semi_join(prot_miss_rank, by = "protein") |>
  mutate(
    protein = fct_reorder(protein, is_missing, .fun = mean),  # order by missingness
    sample  = factor(sample, levels = unique(meta$sample))
  )

p_miss <- ggplot(miss_plot_df, aes(x = sample, y = protein, fill = is_missing)) +
  geom_tile() +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c(`TRUE` = "grey30", `FALSE` = "white"),
    labels = c(`TRUE` = "Missing", `FALSE` = "Detected"),
    name   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 10)
  ) +
  labs(
    title = paste0("Missingness heatmap (top ", top_n, " proteins by missing fraction)"),
    x = NULL,
    y = NULL
  )

print(p_miss)
