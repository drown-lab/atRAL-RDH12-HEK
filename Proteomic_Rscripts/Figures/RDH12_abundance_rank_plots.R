# =========================================================
# RDH12 abundance and ranked proteome plots
# =========================================================

# -------------------------
# 1. Load packages
# -------------------------
packages_needed <- c("dplyr", "ggplot2", "ggrepel", "readr", "stringr", "tibble")

missing_packages <- packages_needed[!vapply(packages_needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(stringr)

# =========================================================
# 2. User-editable parameters
# =========================================================

input_file <- "Proteomic_output_txts/Data_imputed_results.csv"
out_dir <- "Proteomic_Figs/RDH12_abundance_rank_plots"

rdh12_gene <- "RDH12"
rdh12_uniprot <- "Q96NR8"
rdh_family_regex <- "^RDH"

condition_order <- c(
  "RDH12_control_atRAL5hr",
  "RDH12_100_atRAL5hr",
  "RDH12_200_atRAL5hr",
  "RDH12_control_atRAL5hr.24h_recvr",
  "RDH12_100_atRAL5hr.24h_recvr",
  "RDH12_200_atRAL5hr.24h_recvr",
  "GFP_control_atRAL5hr.24h_recvr",
  "GFP_100_atRAL5hr.24h_recvr",
  "GFP_200_atRAL5hr.24h_recvr"
)

condition_labels <- c(
  "RDH12_control_atRAL5hr" = "Acute RDH12 Veh",
  "RDH12_100_atRAL5hr" = "Acute RDH12 100uM",
  "RDH12_200_atRAL5hr" = "Acute RDH12 200uM",
  "RDH12_control_atRAL5hr.24h_recvr" = "Recvr RDH12 Veh",
  "RDH12_100_atRAL5hr.24h_recvr" = "Recvr RDH12 100uM",
  "RDH12_200_atRAL5hr.24h_recvr" = "Recvr RDH12 200uM",
  "GFP_control_atRAL5hr.24h_recvr" = "Recvr WT Veh",
  "GFP_100_atRAL5hr.24h_recvr" = "Recvr WT 100uM",
  "GFP_200_atRAL5hr.24h_recvr" = "Recvr WT 200uM"
)

condition_colors <- c(
  "RDH12_control_atRAL5hr" = "#6B7280",
  "RDH12_100_atRAL5hr" = "#C2410C",
  "RDH12_200_atRAL5hr" = "#EA580C",
  "RDH12_control_atRAL5hr.24h_recvr" = "#374151",
  "RDH12_100_atRAL5hr.24h_recvr" = "#A16207",
  "RDH12_200_atRAL5hr.24h_recvr" = "#CA8A04",
  "GFP_control_atRAL5hr.24h_recvr" = "#2563EB",
  "GFP_100_atRAL5hr.24h_recvr" = "#0891B2",
  "GFP_200_atRAL5hr.24h_recvr" = "#0F766E"
)

plot_dpi <- 300

abundance_plot_width <- 9.5
abundance_plot_height <- 5.8

rank_plot_width <- 12
rank_plot_height <- 8.5

recovery_rank_plot_width <- 10.5
recovery_rank_plot_height <- 7

# =========================================================
# 3. Import data and check required columns
# =========================================================

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

proteome_df <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c(
  "condition",
  "replicate",
  "Genetype",
  "ExpType",
  "Treatment",
  "name",
  "Genes",
  "Protein.Group",
  "ID",
  "intensity"
)

missing_cols <- setdiff(required_cols, colnames(proteome_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

proteome_df <- proteome_df %>%
  mutate(
    condition = as.character(.data$condition),
    condition_label = dplyr::recode(.data$condition, !!!condition_labels, .default = .data$condition),
    condition_factor = factor(.data$condition, levels = condition_order),
    condition_label_factor = factor(.data$condition_label, levels = condition_labels[condition_order]),
    replicate = as.character(.data$replicate),
    Genetype = as.character(.data$Genetype),
    ExpType = as.character(.data$ExpType),
    Treatment = as.character(.data$Treatment),
    name = as.character(.data$name),
    Genes = as.character(.data$Genes),
    Protein.Group = as.character(.data$Protein.Group),
    ID = as.character(.data$ID),
    intensity = suppressWarnings(as.numeric(.data$intensity))
  )

if (any(is.na(proteome_df$intensity))) {
  warning("There are ", sum(is.na(proteome_df$intensity)), " rows with missing or non-numeric intensity values.")
}

unknown_conditions <- setdiff(unique(proteome_df$condition), condition_order)
if (length(unknown_conditions) > 0) {
  warning("Conditions not included in condition_order: ", paste(unknown_conditions, collapse = ", "))
}

if (all(is.na(proteome_df$condition_factor))) {
  stop("No rows match the configured condition_order.")
}

# =========================================================
# 4. Helper functions
# =========================================================

export_plot <- function(plot_obj, file_stub, width, height, dpi = 300) {
  png_file <- file.path(out_dir, paste0(file_stub, ".png"))
  pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))

  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )

  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in"
  )

  c(png_file, pdf_file)
}

rdh12_flag <- function(name, genes, id, protein_group) {
  name == rdh12_gene |
    genes == rdh12_gene |
    id == rdh12_uniprot |
    protein_group == rdh12_uniprot
}

rdh_family_flag <- function(name, genes) {
  coalesce(str_detect(name, rdh_family_regex), FALSE) |
    coalesce(str_detect(genes, rdh_family_regex), FALSE)
}

make_rank_table <- function(df) {
  df %>%
    filter(!is.na(.data$condition_factor), !is.na(.data$intensity)) %>%
    group_by(
      .data$condition,
      .data$condition_label,
      .data$condition_label_factor,
      .data$Genetype,
      .data$ExpType,
      .data$Treatment,
      .data$name,
      .data$Genes,
      .data$Protein.Group,
      .data$ID
    ) %>%
    summarise(
      mean_log2_abundance = mean(.data$intensity, na.rm = TRUE),
      sd_log2_abundance = sd(.data$intensity, na.rm = TRUE),
      n_replicates = sum(!is.na(.data$intensity)),
      .groups = "drop"
    ) %>%
    group_by(.data$condition, .data$condition_label, .data$condition_label_factor) %>%
    arrange(desc(.data$mean_log2_abundance), .data$name, .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      highlight_rdh12 = rdh12_flag(
        name = .data$name,
        genes = .data$Genes,
        id = .data$ID,
        protein_group = .data$Protein.Group
      ),
      highlight_rdh_family = rdh_family_flag(
        name = .data$name,
        genes = .data$Genes
      )
    ) %>%
    ungroup()
}

base_rank_plot <- function(plot_df) {
  ggplot(plot_df, aes(x = .data$rank, y = .data$mean_log2_abundance)) +
    geom_line(color = "grey70", linewidth = 0.45) +
    geom_point(color = "grey78", size = 0.45, alpha = 0.65) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    labs(
      x = "Rank (highest abundance -> lowest)",
      y = "Mean log2 abundance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "grey45"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# =========================================================
# 5. RDH12 abundance across RDH12-expressing cells
# =========================================================

rdh12_abundance_df <- proteome_df %>%
  filter(.data$Genetype == "RDH12") %>%
  filter(
    rdh12_flag(
      name = .data$name,
      genes = .data$Genes,
      id = .data$ID,
      protein_group = .data$Protein.Group
    )
  ) %>%
  filter(!is.na(.data$condition_factor), !is.na(.data$intensity)) %>%
  arrange(.data$condition_factor, .data$replicate)

if (nrow(rdh12_abundance_df) == 0) {
  stop("No RDH12 rows were found in RDH12-expressing conditions.")
}

rdh12_abundance_summary <- rdh12_abundance_df %>%
  group_by(.data$condition, .data$condition_label, .data$condition_label_factor) %>%
  summarise(
    mean_log2_abundance = mean(.data$intensity, na.rm = TRUE),
    sd_log2_abundance = sd(.data$intensity, na.rm = TRUE),
    n_replicates = sum(!is.na(.data$intensity)),
    .groups = "drop"
  ) %>%
  arrange(.data$condition_label_factor)

readr::write_csv(
  rdh12_abundance_df,
  file.path(out_dir, "rdh12_abundance_replicate_values.csv")
)

readr::write_csv(
  rdh12_abundance_summary,
  file.path(out_dir, "rdh12_abundance_condition_summary.csv")
)

rdh12_abundance_plot <- ggplot(
  rdh12_abundance_df,
  aes(x = .data$condition_label_factor, y = .data$intensity)
) +
  geom_line(
    data = rdh12_abundance_summary,
    aes(y = .data$mean_log2_abundance, group = 1),
    color = "#111827",
    linewidth = 0.7
  ) +
  geom_point(
    data = rdh12_abundance_summary,
    aes(y = .data$mean_log2_abundance),
    shape = 21,
    fill = "white",
    color = "#111827",
    size = 3.1,
    stroke = 0.9
  ) +
  geom_point(
    aes(fill = .data$condition),
    shape = 21,
    color = "black",
    size = 2.5,
    alpha = 0.85,
    position = position_jitter(width = 0.08, height = 0, seed = 1)
  ) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  labs(
    title = "RDH12 abundance across RDH12-expressing proteome datasets",
    x = "Condition",
    y = "Log2 abundance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

rdh12_abundance_files <- export_plot(
  plot_obj = rdh12_abundance_plot,
  file_stub = "rdh12_abundance_across_rdh12_conditions",
  width = abundance_plot_width,
  height = abundance_plot_height,
  dpi = plot_dpi
)

# =========================================================
# 6. RDH12 ranked abundance across all conditions
# =========================================================

rank_table <- make_rank_table(proteome_df)

if (nrow(rank_table) == 0) {
  stop("No proteins remain for ranked abundance plotting.")
}

rdh12_rank_table <- rank_table %>%
  filter(.data$highlight_rdh12)

if (nrow(rdh12_rank_table) == 0) {
  stop("RDH12 was not found in the ranked abundance table.")
}

readr::write_csv(
  rank_table,
  file.path(out_dir, "ranked_abundance_all_conditions.csv")
)

readr::write_csv(
  rdh12_rank_table,
  file.path(out_dir, "rdh12_ranked_abundance_positions.csv")
)

rdh12_rank_plot <- base_rank_plot(rank_table) +
  geom_point(
    data = rdh12_rank_table,
    color = "#E11D48",
    size = 2.4
  ) +
  ggrepel::geom_text_repel(
    data = rdh12_rank_table,
    aes(label = .data$name),
    color = "#E11D48",
    size = 3,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    seed = 1,
    max.overlaps = Inf
  ) +
  facet_wrap(~condition_label_factor, scales = "free_x", ncol = 3) +
  labs(
    title = "RDH12 ranked abundance within each proteome condition"
  )

rdh12_rank_files <- export_plot(
  plot_obj = rdh12_rank_plot,
  file_stub = "rdh12_ranked_abundance_all_conditions",
  width = rank_plot_width,
  height = rank_plot_height,
  dpi = plot_dpi
)

# =========================================================
# 7. Recovery RDH-family ranked abundance plots
# =========================================================

recovery_rank_table <- make_rank_table(
  proteome_df %>%
    filter(.data$ExpType == "atRAL5hr+24h_recvr")
)

if (nrow(recovery_rank_table) == 0) {
  stop("No recovery rows remain for ranked abundance plotting.")
}

recovery_rdh_family_table <- recovery_rank_table %>%
  filter(.data$highlight_rdh_family) %>%
  mutate(rdh_label = if_else(!is.na(.data$name) & .data$name != "", .data$name, .data$Genes))

if (nrow(recovery_rdh_family_table) == 0) {
  stop("No RDH-family proteins matched regex: ", rdh_family_regex)
}

readr::write_csv(
  recovery_rank_table,
  file.path(out_dir, "ranked_abundance_recovery_conditions.csv")
)

readr::write_csv(
  recovery_rdh_family_table,
  file.path(out_dir, "recovery_rdh_family_ranked_abundance_positions.csv")
)

recovery_rdh_rank_plot <- base_rank_plot(recovery_rank_table) +
  geom_point(
    data = recovery_rdh_family_table,
    aes(color = .data$rdh_label),
    size = 2.5
  ) +
  ggrepel::geom_text_repel(
    data = recovery_rdh_family_table,
    aes(label = .data$rdh_label, color = .data$rdh_label),
    size = 3,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    seed = 1,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(~condition_label_factor, scales = "free_x", ncol = 3) +
  scale_color_brewer(palette = "Dark2", name = "RDH protein") +
  labs(
    title = "RDH-family ranked abundance in recovery proteome datasets"
  )

recovery_rdh_rank_files <- export_plot(
  plot_obj = recovery_rdh_rank_plot,
  file_stub = "recovery_rdh_family_ranked_abundance",
  width = recovery_rank_plot_width,
  height = recovery_rank_plot_height,
  dpi = plot_dpi
)

# =========================================================
# 8. Completion summary
# =========================================================

rank_counts <- rank_table %>%
  count(.data$condition_label, name = "n_ranked_proteins") %>%
  arrange(match(.data$condition_label, condition_labels[condition_order]))

recovery_rdh_proteins <- recovery_rdh_family_table %>%
  distinct(.data$rdh_label) %>%
  arrange(.data$rdh_label) %>%
  pull(.data$rdh_label)

message("\nRDH12 abundance and ranked proteome plots completed.")
message("Rows read: ", nrow(proteome_df))
message("Unique proteins/IDs read: ", n_distinct(proteome_df$ID))
message("RDH12 replicate rows in RDH12-expressing conditions: ", nrow(rdh12_abundance_df))
message("RDH12 ranked rows: ", nrow(rdh12_rank_table))
message("RDH-family proteins highlighted in recovery plots: ", paste(recovery_rdh_proteins, collapse = ", "))
message("\nRanked proteins per condition:")
print(rank_counts)

message("\nFigure files written:")
for (file in c(rdh12_abundance_files, rdh12_rank_files, recovery_rdh_rank_files)) {
  message("  ", file)
}

message("\nCSV files written to: ", out_dir)
