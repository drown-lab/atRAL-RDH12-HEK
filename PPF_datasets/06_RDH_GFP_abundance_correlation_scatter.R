# =========================================================
# RDH vs GFP abundance correlation scatter plot
# =========================================================

# -------------------------
# 1. Load packages
# -------------------------
packages_needed <- c("dplyr", "ggplot2", "ggrepel", "readr", "stringr", "tibble", "tidyr")

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
library(tidyr)

# =========================================================
# 2. User-editable parameters
# =========================================================

input_dir <- if (dir.exists("PPF_datasets")) "PPF_datasets" else "."
prep_script_path <- file.path(input_dir, "01_Prep_data.R")
out_dir <- file.path(input_dir, "output_correlation")

highlight_genes <- c( "TBP")
highlight_colors <- c(

  "TBP" = "#EE1289"
)

group_order <- c("RDH", "GFP")
min_intensity_for_detected <- 0

plot_width <- 6.5
plot_height <- 5.8
plot_dpi <- 600

out_correlation_table <- file.path(out_dir, "RDH_vs_GFP_abundance_correlation_table.csv")
out_summary_table <- file.path(out_dir, "RDH_vs_GFP_correlation_summary.csv")
out_highlight_table <- file.path(out_dir, "RDH_vs_GFP_highlighted_proteins.csv")
out_png <- file.path(out_dir, "RDH_vs_GFP_abundance_correlation_scatter.png")
out_pdf <- file.path(out_dir, "RDH_vs_GFP_abundance_correlation_scatter.pdf")

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

PPF_correlation_input <- PPFtable_v4 %>%
  mutate(
    sample_name = as.character(.data$sample_name),
    Group = str_extract(.data$sample_name, "^(RDH|GFP)"),
    Protein.Group = as.character(.data[["Protein.Group"]]),
    Genes = as.character(.data$Genes),
    PG.MaxLFQ = suppressWarnings(as.numeric(.data$PG.MaxLFQ))
  ) %>%
  filter(.data$Group %in% group_order)

if (any(is.na(PPF_correlation_input$PG.MaxLFQ))) {
  warning("There are ", sum(is.na(PPF_correlation_input$PG.MaxLFQ)), " rows with missing or non-numeric PG.MaxLFQ values.")
}

missing_groups <- setdiff(group_order, unique(PPF_correlation_input$Group))
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

# =========================================================
# 5. Average abundance by group and calculate correlations
# =========================================================

group_abundance <- PPF_correlation_input %>%
  filter(!is.na(.data$PG.MaxLFQ), .data$PG.MaxLFQ > min_intensity_for_detected) %>%
  mutate(
    log2_abundance = log2(.data$PG.MaxLFQ),
    Group = factor(.data$Group, levels = group_order)
  ) %>%
  group_by(.data[["Protein.Group"]], .data$Genes, .data$Group) %>%
  summarise(
    mean_log2_abundance = mean(.data$log2_abundance, na.rm = TRUE),
    sd_log2_abundance = sd(.data$log2_abundance, na.rm = TRUE),
    n_replicates = n_distinct(.data$sample_name),
    .groups = "drop"
  )

correlation_table <- group_abundance %>%
  select(Protein.Group, Genes, Group, mean_log2_abundance, n_replicates) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(mean_log2_abundance, n_replicates)
  ) %>%
  filter(!is.na(.data$mean_log2_abundance_RDH), !is.na(.data$mean_log2_abundance_GFP)) %>%
  mutate(
    highlight_gene = find_highlight_gene(.data$Genes),
    is_highlight = !is.na(.data$highlight_gene),
    abundance_difference_RDH_minus_GFP = .data$mean_log2_abundance_RDH - .data$mean_log2_abundance_GFP
  )

if (nrow(correlation_table) < 3) {
  stop("Fewer than 3 shared proteins were available for correlation.")
}

pearson_test <- cor.test(
  correlation_table$mean_log2_abundance_GFP,
  correlation_table$mean_log2_abundance_RDH,
  method = "pearson"
)

spearman_test <- cor.test(
  correlation_table$mean_log2_abundance_GFP,
  correlation_table$mean_log2_abundance_RDH,
  method = "spearman",
  exact = FALSE
)

correlation_summary <- tibble(
  comparison = "RDH_vs_GFP",
  n_shared_proteins = nrow(correlation_table),
  pearson_r = unname(pearson_test$estimate),
  pearson_p_value = pearson_test$p.value,
  spearman_rho = unname(spearman_test$estimate),
  spearman_p_value = spearman_test$p.value
)

highlight_table <- correlation_table %>%
  filter(.data$is_highlight) %>%
  arrange(.data$highlight_gene)

readr::write_csv(correlation_table, out_correlation_table)
readr::write_csv(correlation_summary, out_summary_table)
readr::write_csv(highlight_table, out_highlight_table)

# =========================================================
# 6. Plot abundance correlation
# =========================================================

axis_limits <- range(
  c(
    correlation_table$mean_log2_abundance_GFP,
    correlation_table$mean_log2_abundance_RDH
  ),
  na.rm = TRUE
)

axis_padding <- diff(axis_limits) * 0.04
axis_limits <- c(axis_limits[1] - axis_padding, axis_limits[2] + axis_padding)

subtitle_text <- paste0(
  "Pearson r = ", round(correlation_summary$pearson_r, 3),
  "; Spearman rho = ", round(correlation_summary$spearman_rho, 3),
  "; shared proteins = ", correlation_summary$n_shared_proteins
)

correlation_plot <- ggplot(
  correlation_table,
  aes(x = .data$mean_log2_abundance_GFP, y = .data$mean_log2_abundance_RDH)
) +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", size = 0.45) +
  geom_point(color = "grey70", size = 0.7, alpha = 0.65) +
  geom_smooth(method = "lm", se = FALSE, color = "grey20", size = 0.7) +
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
  scale_color_manual(values = highlight_colors, breaks = highlight_genes, name = "Highlighted protein") +
  coord_equal(xlim = axis_limits, ylim = axis_limits) +
  labs(
    title = "RDH vs GFP protein abundance correlation",
    subtitle = subtitle_text,
    x = "WT  log2 abundance",
    y = "RDH12  log2 abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = out_png,
  plot = correlation_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = plot_dpi
)

ggsave(
  filename = out_pdf,
  plot = correlation_plot,
  width = plot_width,
  height = plot_height,
  units = "in"
)

# =========================================================
# 7. Completion summary
# =========================================================

message("\nRDH/GFP abundance correlation scatter plot completed.")
message("Rows read from PPFtable_v4: ", nrow(PPFtable_v4))
message("Shared proteins used for correlation: ", nrow(correlation_table))
message("Pearson r: ", round(correlation_summary$pearson_r, 4), " (p = ", signif(correlation_summary$pearson_p_value, 4), ")")
message("Spearman rho: ", round(correlation_summary$spearman_rho, 4), " (p = ", signif(correlation_summary$spearman_p_value, 4), ")")

if (nrow(highlight_table) > 0) {
  message("Highlighted proteins found: ", paste(highlight_table$highlight_gene, collapse = ", "))
} else {
  message("No configured highlighted proteins were found in both RDH and GFP.")
}

message("\nFiles written:")
for (file in c(out_png, out_pdf, out_correlation_table, out_summary_table, out_highlight_table)) {
  message("  ", file)
}
