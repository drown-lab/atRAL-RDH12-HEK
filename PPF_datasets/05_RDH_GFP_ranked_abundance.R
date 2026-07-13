# =========================================================
# RDH/GFP ranked abundance plot with protein highlights
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

input_dir <- if (dir.exists("PPF_datasets")) "PPF_datasets" else "."
prep_script_path <- file.path(input_dir, "01_Prep_data.R")
out_dir <- file.path(input_dir, "output_ranked_abundance")

highlight_genes <- c("RDH12", "TBCD", "GAPDH", "TBP")
highlight_colors <- c(
  "RDH12" = "#E11D48",
  "TBCD" = "#2563EB",
  "GAPDH" = "#059669",
  "TBP" = "#EE1289"
)

group_order <- c("RDH", "GFP")
min_intensity_for_detected <- 0

plot_width <- 7
plot_height <- 5
plot_dpi <- 600

out_rank_table <- file.path(out_dir, "RDH_GFP_ranked_abundance_table.csv")
out_highlight_table <- file.path(out_dir, "RDH_GFP_highlight_rank_positions.csv")
out_missing_table <- file.path(out_dir, "RDH_GFP_missing_highlight_genes.csv")
out_png <- file.path(out_dir, "RDH_GFP_ranked_abundance_highlights.png")
out_pdf <- file.path(out_dir, "RDH_GFP_ranked_abundance_highlights.pdf")

# =========================================================
# 3. Import data and check required columns
# =========================================================

if (!file.exists(prep_script_path)) {
  prep_script_path <- "01_Prep_data.R"
}

if (!file.exists(prep_script_path)) {
  stop("Could not find 01_Prep_data.R. Run this script from the project root or PPF_datasets folder.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(prep_script_path)

if (!exists("PPFtable_v4")) {
  stop("PPFtable_v4 was not created by 01_Prep_data.R.")
}

required_cols <- c("sample_name", "Protein.Group", "Genes", "PG.MaxLFQ")
missing_cols <- setdiff(required_cols, colnames(PPFtable_v4))
if (length(missing_cols) > 0) {
  stop("PPFtable_v4 is missing required columns: ", paste(missing_cols, collapse = ", "))
}

PPF_rank_input <- PPFtable_v4 %>%
  mutate(
    sample_name = as.character(.data$sample_name),
    Group = str_extract(.data$sample_name, "^(RDH|GFP)"),
    Protein.Group = as.character(.data[["Protein.Group"]]),
    Genes = as.character(.data$Genes),
    PG.MaxLFQ = suppressWarnings(as.numeric(.data$PG.MaxLFQ))
  ) %>%
  filter(.data$Group %in% group_order)

if (any(is.na(PPF_rank_input$PG.MaxLFQ))) {
  warning("There are ", sum(is.na(PPF_rank_input$PG.MaxLFQ)), " rows with missing or non-numeric PG.MaxLFQ values.")
}

missing_groups <- setdiff(group_order, unique(PPF_rank_input$Group))
if (length(missing_groups) > 0) {
  stop("Missing expected sample groups: ", paste(missing_groups, collapse = ", "))
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

find_highlight_gene_one <- function(genes) {
  tokens <- split_gene_tokens(genes)
  hit <- highlight_genes[highlight_genes %in% tokens]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[[1]]
}

find_highlight_gene <- function(genes) {
  vapply(genes, find_highlight_gene_one, character(1), USE.NAMES = FALSE)
}

make_rank_table <- function(df) {
  df %>%
    filter(!is.na(.data$PG.MaxLFQ), .data$PG.MaxLFQ > min_intensity_for_detected) %>%
    mutate(
      log2_abundance = log2(.data$PG.MaxLFQ),
      Group = factor(.data$Group, levels = group_order)
    ) %>%
    group_by(.data$Group, .data[["Protein.Group"]], .data$Genes) %>%
    summarise(
      mean_log2_abundance = mean(.data$log2_abundance, na.rm = TRUE),
      sd_log2_abundance = sd(.data$log2_abundance, na.rm = TRUE),
      n_replicates = n_distinct(.data$sample_name),
      .groups = "drop"
    ) %>%
    mutate(
      highlight_gene = find_highlight_gene(.data$Genes),
      is_highlight = !is.na(.data$highlight_gene)
    ) %>%
    group_by(.data$Group) %>%
    arrange(desc(.data$mean_log2_abundance), .data$Genes, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

# =========================================================
# 5. Rank proteins and identify highlighted genes
# =========================================================

rank_table <- make_rank_table(PPF_rank_input)

if (nrow(rank_table) == 0) {
  stop("No proteins remain for ranked abundance plotting.")
}

highlight_table <- rank_table %>%
  filter(.data$is_highlight) %>%
  arrange(.data$Group, .data$highlight_gene, .data$rank)

missing_highlight_genes <- tidyr::expand_grid(
  Group = factor(group_order, levels = group_order),
  highlight_gene = highlight_genes
) %>%
  anti_join(
    highlight_table %>% distinct(.data$Group, .data$highlight_gene),
    by = c("Group", "highlight_gene")
  )

readr::write_csv(rank_table, out_rank_table)
readr::write_csv(highlight_table, out_highlight_table)
readr::write_csv(missing_highlight_genes, out_missing_table)

# =========================================================
# 6. Plot ranked abundance distributions
# =========================================================

rank_plot <- ggplot(rank_table, aes(x = .data$rank, y = .data$mean_log2_abundance)) +
  geom_line(color = "grey70", size = 0.45) +
  geom_point(color = "grey78", size = 0.45, alpha = 0.65) +
  geom_point(
    data = highlight_table,
    aes(color = .data$highlight_gene),
    size = 2.8,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = highlight_table,
    aes(label = .data$highlight_gene, color = .data$highlight_gene),
    size = 3.2,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    seed = 1,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(~Group, scales = "free_x", ncol = 2) +
  scale_color_manual(values = highlight_colors, breaks = highlight_genes, name = "Highlighted protein") +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(
    title = "RDH vs GFP ranked protein abundance",
    subtitle = paste0("Replicates averaged within group; detected proteins have PG.MaxLFQ > ", min_intensity_for_detected),
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
  count(.data$Group, name = "n_ranked_proteins") %>%
  arrange(.data$Group)

found_genes <- highlight_table %>%
  distinct(.data$highlight_gene) %>%
  arrange(.data$highlight_gene) %>%
  pull(.data$highlight_gene)

message("\nRDH/GFP ranked abundance plot completed.")
message("Rows read from PPFtable_v4: ", nrow(PPFtable_v4))
message("Detected rows used for ranking: ", nrow(PPF_rank_input %>% filter(!is.na(.data$PG.MaxLFQ), .data$PG.MaxLFQ > min_intensity_for_detected)))
message("Ranked protein-group rows: ", nrow(rank_table))
message("Highlighted protein-group rows: ", nrow(highlight_table))
message("Highlighted genes found: ", paste(found_genes, collapse = ", "))

if (nrow(missing_highlight_genes) > 0) {
  message("Highlighted genes missing from at least one group:")
  print(missing_highlight_genes)
} else {
  message("All configured highlighted genes were found in both groups.")
}

message("\nRanked proteins per group:")
print(rank_counts)

message("\nFiles written:")
for (file in c(out_png, out_pdf, out_rank_table, out_highlight_table, out_missing_table)) {
  message("  ", file)
}
