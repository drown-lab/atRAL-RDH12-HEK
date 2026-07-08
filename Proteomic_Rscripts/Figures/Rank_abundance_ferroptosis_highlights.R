# =========================================================
# Rank-abundance plots with ferroptosis/redox protein highlights
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
library(tibble)

# =========================================================
# 2. User-editable parameters
# =========================================================

input_file <- "Proteomic_output_txts/Data_imputed_results.csv"
out_dir <- "Proteomic_Figs/rank_abundance_ferroptosis_highlights"

highlight_genes <- c(
  "CYBA", "GSTA4", "NFE2L2", "GSTZ1", "TXNRD1", "GSS", "GSR", "FTL",
  "GPX1", "TXNRD2", "HMOX1", "FTH1", "GCLM", "PRDX6", "ACSL4", "ACSL3",
  "GCLC", "MGST2", "KEAP1", "MGST3", "SLC7A11", "TFRC", "SLC11A2"
)

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

plot_width <- 13
plot_height <- 9
plot_dpi <- 300

out_rank_table <- file.path(out_dir, "ranked_abundance_all_conditions.csv")
out_highlight_table <- file.path(out_dir, "highlighted_gene_rank_positions.csv")
out_missing_table <- file.path(out_dir, "missing_highlight_genes.csv")
out_png <- file.path(out_dir, "rank_abundance_ferroptosis_highlights.png")
out_pdf <- file.path(out_dir, "rank_abundance_ferroptosis_highlights.pdf")

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

split_gene_tokens <- function(x) {
  x <- toupper(trimws(x))
  x <- x[!is.na(x) & x != ""]

  if (length(x) == 0) {
    return(character(0))
  }

  tokens <- unlist(strsplit(x, "[;,|[:space:]]+"), use.names = FALSE)
  unique(tokens[tokens != ""])
}

find_highlight_gene_one <- function(name, genes) {
  tokens <- unique(c(split_gene_tokens(name), split_gene_tokens(genes)))
  hit <- highlight_genes[highlight_genes %in% tokens]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[[1]]
}

find_highlight_gene <- function(name, genes) {
  mapply(
    find_highlight_gene_one,
    name = name,
    genes = genes,
    USE.NAMES = FALSE
  )
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
    mutate(
      highlight_gene = find_highlight_gene(.data$name, .data$Genes),
      is_highlight = !is.na(.data$highlight_gene)
    ) %>%
    group_by(.data$condition, .data$condition_label, .data$condition_label_factor) %>%
    arrange(desc(.data$mean_log2_abundance), .data$name, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

# =========================================================
# 5. Rank proteins and identify highlighted genes
# =========================================================

rank_table <- make_rank_table(proteome_df)

if (nrow(rank_table) == 0) {
  stop("No proteins remain for ranked abundance plotting.")
}

highlight_table <- rank_table %>%
  filter(.data$is_highlight) %>%
  arrange(.data$condition_label_factor, .data$highlight_gene, .data$rank)

if (nrow(highlight_table) == 0) {
  stop("None of the configured highlight genes were found in name or Genes.")
}

missing_highlight_genes <- tibble(highlight_gene = highlight_genes) %>%
  anti_join(
    highlight_table %>% distinct(.data$highlight_gene),
    by = "highlight_gene"
  )

readr::write_csv(rank_table, out_rank_table)
readr::write_csv(highlight_table, out_highlight_table)
readr::write_csv(missing_highlight_genes, out_missing_table)

# =========================================================
# 6. Plot ranked abundance distributions
# =========================================================

highlight_colors <- setNames(
  grDevices::hcl.colors(length(highlight_genes), palette = "Dark 3"),
  highlight_genes
)

rank_plot <- ggplot(rank_table, aes(x = .data$rank, y = .data$mean_log2_abundance)) +
  geom_line(color = "grey70", linewidth = 0.45) +
  geom_point(color = "grey78", size = 0.45, alpha = 0.65) +
  geom_point(
    data = highlight_table,
    aes(color = .data$highlight_gene),
    size = 2.2,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = highlight_table,
    aes(label = .data$highlight_gene, color = .data$highlight_gene),
    size = 2.7,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.15,
    seed = 1,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(~condition_label_factor, scales = "free_x", ncol = 3) +
  scale_color_manual(values = highlight_colors, breaks = highlight_genes, name = "Highlighted protein") +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(
    title = "Rank-abundance distributions with ferroptosis/redox protein highlights",
    x = "Rank (highest abundance -> lowest)",
    y = "Mean log2 abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey45"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.key.width = grid::unit(0.8, "lines")
  ) +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 2.8, alpha = 1)))

ggsave(
  filename = out_png,
  plot = rank_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = plot_dpi
)

ggsave(
  filename = out_pdf,
  plot = rank_plot,
  width = plot_width,
  height = plot_height,
  units = "in"
)

# =========================================================
# 7. Completion summary
# =========================================================

rank_counts <- rank_table %>%
  count(.data$condition_label, name = "n_ranked_proteins") %>%
  arrange(match(.data$condition_label, condition_labels[condition_order]))

found_genes <- highlight_table %>%
  distinct(.data$highlight_gene) %>%
  arrange(.data$highlight_gene) %>%
  pull(.data$highlight_gene)

message("\nRank-abundance ferroptosis/redox highlight plots completed.")
message("Rows read: ", nrow(proteome_df))
message("Unique proteins/IDs read: ", n_distinct(proteome_df$ID))
message("Ranked protein-condition rows: ", nrow(rank_table))
message("Highlighted protein-condition rows: ", nrow(highlight_table))
message("Highlighted genes found: ", paste(found_genes, collapse = ", "))

if (nrow(missing_highlight_genes) > 0) {
  message("Highlighted genes not found: ", paste(missing_highlight_genes$highlight_gene, collapse = ", "))
} else {
  message("All configured highlighted genes were found.")
}

message("\nRanked proteins per condition:")
print(rank_counts)

message("\nFiles written:")
for (file in c(out_png, out_pdf, out_rank_table, out_highlight_table, out_missing_table)) {
  message("  ", file)
}
