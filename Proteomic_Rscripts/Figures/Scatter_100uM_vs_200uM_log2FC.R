# =========================================================
# 100uM vs 200uM log2FC scatterplot for acute RDH12 atRAL
# - x-axis: RDH12 100uM vs vehicle log2FC
# - y-axis: RDH12 200uM vs vehicle log2FC
# - colors show significant direction classes in each dose
# =========================================================

# -------------------------
# 1. Load packages
# -------------------------
packages_needed <- c("dplyr", "ggplot2", "ggrepel", "readr", "stringr", "svglite")

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

input_file <- "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv"
out_dir <- "Proteomic_Figs/Overlap_100uM_vs_200uM"

padj_cutoff <- 0.01
lfc_cutoff <- log2(2)
n_top_labels <- 12

contrast_100 <- "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr"
contrast_200 <- "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr"

out_table <- file.path(out_dir, "scatter_100uM_vs_200uM_log2FC_table.csv")
out_svg <- file.path(out_dir, "scatter_100uM_vs_200uM_log2FC.svg")
out_pdf <- file.path(out_dir, "scatter_100uM_vs_200uM_log2FC.pdf")

# =========================================================
# 3. Import data and check required columns
# =========================================================

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

DEPresults_v2 <- readr::read_csv(input_file, show_col_types = FALSE)

ratio_100_col <- paste0(contrast_100, "_ratio")
ratio_200_col <- paste0(contrast_200, "_ratio")
padj_100_col <- paste0(contrast_100, "_p.adj")
padj_200_col <- paste0(contrast_200, "_p.adj")

required_cols <- c("name", "ID", ratio_100_col, ratio_200_col, padj_100_col, padj_200_col)
missing_cols <- setdiff(required_cols, colnames(DEPresults_v2))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# =========================================================
# 4. Prepare plotting table
# =========================================================

call_direction <- function(ratio, padj, padj_cutoff, lfc_cutoff) {
  dplyr::case_when(
    !is.na(padj) & !is.na(ratio) & padj <= padj_cutoff & ratio >= lfc_cutoff ~ "up",
    !is.na(padj) & !is.na(ratio) & padj <= padj_cutoff & ratio <= -lfc_cutoff ~ "down",
    TRUE ~ "ns"
  )
}

plot_df <- DEPresults_v2 %>%
  transmute(
    name = as.character(.data$name),
    ID = as.character(.data$ID),
    log2fc_100uM_vs_vehicle = suppressWarnings(as.numeric(.data[[ratio_100_col]])),
    log2fc_200uM_vs_vehicle = suppressWarnings(as.numeric(.data[[ratio_200_col]])),
    padj_100uM_vs_vehicle = suppressWarnings(as.numeric(.data[[padj_100_col]])),
    padj_200uM_vs_vehicle = suppressWarnings(as.numeric(.data[[padj_200_col]]))
  ) %>%
  filter(
    !is.na(log2fc_100uM_vs_vehicle),
    !is.na(log2fc_200uM_vs_vehicle)
  ) %>%
  mutate(
    call_100uM = call_direction(
      ratio = log2fc_100uM_vs_vehicle,
      padj = padj_100uM_vs_vehicle,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    ),
    call_200uM = call_direction(
      ratio = log2fc_200uM_vs_vehicle,
      padj = padj_200uM_vs_vehicle,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    ),
    direction_class = case_when(
      call_100uM == "up" & call_200uM == "up" ~ "Up in both",
      call_100uM == "down" & call_200uM == "down" ~ "Down in both",
      call_100uM == "up" & call_200uM == "down" ~ "Up 100 / Down 200",
      call_100uM == "down" & call_200uM == "up" ~ "Down 100 / Up 200",
      call_100uM == "up" & call_200uM == "ns" ~ "100uM up only",
      call_100uM == "down" & call_200uM == "ns" ~ "100uM down only",
      call_100uM == "ns" & call_200uM == "up" ~ "200uM up only",
      call_100uM == "ns" & call_200uM == "down" ~ "200uM down only",
      TRUE ~ "Not significant"
    ),
    significant_either_dose = call_100uM != "ns" | call_200uM != "ns",
    combined_abs_log2fc = abs(log2fc_100uM_vs_vehicle) + abs(log2fc_200uM_vs_vehicle)
  )

if (nrow(plot_df) == 0) {
  stop("No rows remain after removing missing log2FC values.")
}

direction_levels <- c(
  "Up in both",
  "Down in both",
  "Up 100 / Down 200",
  "Down 100 / Up 200",
  "100uM up only",
  "100uM down only",
  "200uM up only",
  "200uM down only",
  "Not significant"
)

plot_df <- plot_df %>%
  mutate(direction_class = factor(direction_class, levels = direction_levels))

label_df <- plot_df %>%
  filter(significant_either_dose) %>%
  arrange(desc(combined_abs_log2fc)) %>%
  slice_head(n = n_top_labels)

# =========================================================
# 5. Export plotting table
# =========================================================

readr::write_csv(plot_df, out_table)

# =========================================================
# 6. Make scatterplot
# =========================================================

direction_colors <- c(
  "Up in both" = "#D73027",
  "Down in both" = "#4575B4",
  "Up 100 / Down 200" = "#E08214",
  "Down 100 / Up 200" = "#5E3C99",
  "100uM up only" = "#F46D43",
  "100uM down only" = "#74ADD1",
  "200uM up only" = "#A50026",
  "200uM down only" = "#313695",
  "Not significant" = "grey78"
)

axis_limit <- max(
  abs(c(plot_df$log2fc_100uM_vs_vehicle, plot_df$log2fc_200uM_vs_vehicle, lfc_cutoff)),
  na.rm = TRUE
)
axis_limit <- ceiling(axis_limit * 10) / 10

p <- ggplot(
  plot_df,
  aes(
    x = log2fc_100uM_vs_vehicle,
    y = log2fc_200uM_vs_vehicle
  )
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), color = "grey35", linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), color = "grey35", linetype = "dashed", linewidth = 0.4) +
  geom_point(
    data = plot_df %>% filter(direction_class == "Not significant"),
    aes(color = direction_class),
    alpha = 0.55,
    size = 2.2
  ) +
  geom_point(
    data = plot_df %>% filter(direction_class != "Not significant"),
    aes(color = direction_class),
    alpha = 0.88,
    size = 2.6
  ) +
  geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dotted", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = name),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.25,
    seed = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = direction_colors, drop = FALSE) +
  coord_equal(
    xlim = c(-axis_limit, axis_limit),
    ylim = c(-axis_limit, axis_limit)
  ) +
  labs(
    title = "Acute RDH12 atRAL: 100uM vs 200uM protein log2FC",
    subtitle = paste0("Colored if padj <= ", padj_cutoff, " and abs(log2FC) >= log2(2) in either dose"),
    x = "log2FC: RDH12 100uM vs vehicle",
    y = "log2FC: RDH12 200uM vs vehicle",
    color = "Direction class"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(out_svg, plot = p, width = 8, height = 7, units = "in")
ggsave(out_pdf, plot = p, width = 8, height = 7, units = "in")

# =========================================================
# 7. Completion summary
# =========================================================

summary_counts <- plot_df %>%
  count(direction_class, name = "n") %>%
  arrange(direction_class)

print(summary_counts)

message("\nScatterplot script completed.")
message("Rows plotted: ", nrow(plot_df))
message("Labeled proteins: ", nrow(label_df))
message("Table written to: ", out_table)
message("SVG written to: ", out_svg)
message("PDF written to: ", out_pdf)
