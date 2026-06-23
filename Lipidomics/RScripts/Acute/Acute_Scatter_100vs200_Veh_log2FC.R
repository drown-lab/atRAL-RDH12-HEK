# =========================================================
# Acute lipid scatterplot:
# 100uM vs vehicle log2FC compared with 200uM vs vehicle log2FC
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------
acute_file <- "Lipidomics/output_txts/Acute_DEA_results_v1_corrected.csv"

out_dir <- "Lipidomics/Figures/Acute/Scatter_100vs200/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.1
fc_cutoff <- 1.3
lfc_cutoff <- log2(fc_cutoff)
n_top_labels <- 10

contrast_100 <- "RDH12_100_vs_RDH12_Veh"
contrast_200 <- "RDH12_200_vs_RDH12_Veh"

out_table <- file.path(out_dir, "Acute_scatter_100vsVeh_vs_200vsVeh_table.csv")
out_pdf <- file.path(out_dir, "Acute_scatter_100vsVeh_vs_200vsVeh_main_text.pdf")
out_png <- file.path(out_dir, "Acute_scatter_100vsVeh_vs_200vsVeh_main_text.png")

# ---------------------------------------------------------
# 2. READ DATA
# ---------------------------------------------------------
lipid_df <- read.csv(
  acute_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ratio_100_col <- paste0(contrast_100, "_ratio")
ratio_200_col <- paste0(contrast_200, "_ratio")
padj_100_col <- paste0(contrast_100, "_p.val")
padj_200_col <- paste0(contrast_200, "_p.val")

required_cols <- c("name", "ID", ratio_100_col, ratio_200_col, padj_100_col, padj_200_col)
missing_cols <- setdiff(required_cols, colnames(lipid_df))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 3. PREPARE PLOTTING TABLE
# ---------------------------------------------------------
call_direction <- function(ratio, padj, padj_cutoff, lfc_cutoff) {
  case_when(
    !is.na(padj) & !is.na(ratio) & padj <= padj_cutoff & ratio >= lfc_cutoff ~ "Up",
    !is.na(padj) & !is.na(ratio) & padj <= padj_cutoff & ratio <= -lfc_cutoff ~ "Down",
    TRUE ~ "NS"
  )
}

plot_df <- lipid_df %>%
  transmute(
    name = str_trim(as.character(name)),
    ID = as.character(ID),
    log2fc_100_vs_vehicle = suppressWarnings(as.numeric(.data[[ratio_100_col]])),
    log2fc_200_vs_vehicle = suppressWarnings(as.numeric(.data[[ratio_200_col]])),
    pvalue_100_vs_vehicle = suppressWarnings(as.numeric(.data[[padj_100_col]])),
    pvalue_200_vs_vehicle = suppressWarnings(as.numeric(.data[[padj_200_col]]))
  ) %>%
  filter(
    !is.na(log2fc_100_vs_vehicle),
    !is.na(log2fc_200_vs_vehicle)
  ) %>%
  mutate(
    direction_100 = call_direction(
      ratio = log2fc_100_vs_vehicle,
      padj = pvalue_100_vs_vehicle,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    ),
    direction_200 = call_direction(
      ratio = log2fc_200_vs_vehicle,
      padj = pvalue_200_vs_vehicle,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    ),
    sig_100 = direction_100 != "NS",
    sig_200 = direction_200 != "NS",
    significance_group = case_when(
      sig_100 & sig_200 ~ "Significant in both",
      sig_100 & !sig_200 ~ "100uM only",
      !sig_100 & sig_200 ~ "200uM only",
      TRUE ~ "Not significant"
    ),
    same_direction = case_when(
      direction_100 == direction_200 & direction_100 != "NS" ~ "same",
      sig_100 | sig_200 ~ "dose-specific or opposite",
      TRUE ~ "not significant"
    ),
    combined_abs_log2fc = abs(log2fc_100_vs_vehicle) + abs(log2fc_200_vs_vehicle),
    min_pvalue = pmin(pvalue_100_vs_vehicle, pvalue_200_vs_vehicle, na.rm = TRUE)
  )

if (nrow(plot_df) == 0) {
  stop("No rows remain after removing missing log2FC values.")
}

significance_levels <- c(
  "Significant in both",
  "100uM only",
  "200uM only",
  "Not significant"
)

plot_df <- plot_df %>%
  mutate(significance_group = factor(significance_group, levels = significance_levels))

label_df <- plot_df %>%
  filter(sig_100 | sig_200) %>%
  arrange(desc(combined_abs_log2fc), min_pvalue) %>%
  slice_head(n = n_top_labels)

write.csv(plot_df, out_table, row.names = FALSE)

# ---------------------------------------------------------
# 4. MAKE MAIN TEXT SCATTERPLOT
# ---------------------------------------------------------
axis_limit <- max(
  abs(c(plot_df$log2fc_100_vs_vehicle, plot_df$log2fc_200_vs_vehicle, lfc_cutoff)),
  na.rm = TRUE
)
axis_limit <- ceiling(axis_limit * 10) / 10

cor_pearson <- cor(
  plot_df$log2fc_100_vs_vehicle,
  plot_df$log2fc_200_vs_vehicle,
  method = "pearson",
  use = "complete.obs"
)

group_colors <- c(
  "Significant in both" = "#6A3D9A",
  "100uM only" = "#E08214",
  "200uM only" = "#0072B2",
  "Not significant" = "grey78"
)

p_scatter <- ggplot(
  plot_df,
  aes(x = log2fc_100_vs_vehicle, y = log2fc_200_vs_vehicle)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_hline(
    yintercept = c(-lfc_cutoff, lfc_cutoff),
    color = "grey35",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    color = "grey35",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dotted", linewidth = 0.45) +
  geom_point(
    data = plot_df %>% filter(significance_group == "Not significant"),
    aes(color = significance_group),
    alpha = 0.45,
    size = 1.8
  ) +
  geom_point(
    data = plot_df %>% filter(significance_group != "Not significant"),
    aes(color = significance_group),
    alpha = 0.9,
    size = 2.6
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = name),
    size = 2.8,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey45",
    segment.size = 0.25,
    seed = 1,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -axis_limit,
    y = axis_limit,
    hjust = 0,
    vjust = 1,
    label = paste0("Pearson r = ", sprintf("%.2f", cor_pearson)),
    size = 3.2,
    color = "grey20"
  ) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  coord_equal(
    xlim = c(-axis_limit, axis_limit),
    ylim = c(-axis_limit, axis_limit),
    clip = "off"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color = "grey20"),
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey25"),
    plot.margin = margin(6, 14, 6, 6)
  ) +
  labs(
    title = "Acute lipid response by atRAL dose",
    subtitle = paste0("Dashed lines: p <= ", padj_cutoff, " and |log2FC| >= log2(", fc_cutoff, ")"),
    x = "log2FC: 100uM vs vehicle",
    y = "log2FC: 200uM vs vehicle",
    color = NULL
  )

print(p_scatter)

ggsave(out_pdf, p_scatter, width = 5.8, height = 5.1, useDingbats = FALSE)
ggsave(out_png, p_scatter, width = 5.8, height = 5.1, dpi = 600)

# ---------------------------------------------------------
# 5. CONSOLE SUMMARY
# ---------------------------------------------------------
summary_counts <- plot_df %>%
  count(significance_group, name = "n") %>%
  arrange(significance_group)

print(summary_counts)

cat("\nSaved files:\n")
cat(out_table, "\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
