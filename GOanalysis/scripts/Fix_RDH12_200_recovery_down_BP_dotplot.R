# =========================================================
# Fix one GO BP simplified dotplot label/layout
# =========================================================

library(dplyr)
library(forcats)
library(ggplot2)
library(paletteer)
library(readr)

# =========================================================
# 1. USER-EDITABLE PARAMETERS
# =========================================================

contrast_name <- "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr_down"

out_dir <- file.path(
  "GOanalysis/output/clusterProfiler",
  contrast_name
)

input_csv <- file.path(out_dir, "enrichGO_BP_simplified.csv")
output_pdf <- file.path(out_dir, "enrichGO_BP_simplified_dotplot.pdf")
output_png <- file.path(out_dir, "enrichGO_BP_simplified_dotplot.png")

terms_to_plot <- 10
plot_width <- 11
plot_height <- 5.5
plot_dpi <- 900

long_label <- "positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction"
short_label <- "positive regulation of PI3K/PKBsignal transduction"

# =========================================================
# 2. READ AND CHECK INPUT
# =========================================================

if (!file.exists(input_csv)) {
  stop("Missing target simplified GO file: ", input_csv)
}

ego_bp_s_df <- read_csv(input_csv, show_col_types = FALSE)

required_cols <- c("Description", "FoldEnrichment", "Count", "p.adjust")
missing_cols <- setdiff(required_cols, colnames(ego_bp_s_df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (nrow(ego_bp_s_df) == 0) {
  stop("The simplified GO table is empty: ", input_csv)
}

if (any(is.na(ego_bp_s_df$FoldEnrichment))) {
  stop("FoldEnrichment contains missing values.")
}

if (any(is.na(ego_bp_s_df$Count))) {
  stop("Count contains missing values.")
}

if (any(is.na(ego_bp_s_df$p.adjust))) {
  stop("p.adjust contains missing values.")
}

# =========================================================
# 3. MAKE PLOT TABLE
# =========================================================

plot_df <- ego_bp_s_df %>%
  mutate(
    FoldEnrichment_num = as.numeric(FoldEnrichment),
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust),
    Description = if_else(Description == long_label, short_label, Description)
  ) %>%
  arrange(p.adjust, desc(Count), desc(FoldEnrichment_num)) %>%
  slice_head(n = terms_to_plot) %>%
  mutate(
    Description = fct_reorder(Description, FoldEnrichment_num)
  )

if (any(is.na(plot_df$FoldEnrichment_num))) {
  stop("FoldEnrichment could not be converted to numeric for at least one plotted term.")
}

# =========================================================
# 4. REBUILD AND SAVE ONLY THE TARGET DOTPLOT
# =========================================================

p_dot_simplified <- ggplot(
  plot_df,
  aes(x = FoldEnrichment_num, y = Description, size = Count, color = p.adjust)
) +
  geom_point() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 11)
  ) +
  labs(
    title = paste0(contrast_name, ": GO BP (simplified)"),
    x = "Fold enrichment",
    y = NULL,
    color = "p.adjust"
  ) +
  scale_color_paletteer_c("grDevices::Plasma", direction = -1)

ggsave(
  filename = output_pdf,
  plot = p_dot_simplified,
  width = plot_width,
  height = plot_height
)

ggsave(
  filename = output_png,
  plot = p_dot_simplified,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

cat("Rebuilt simplified GO BP dotplot for:\n")
cat(contrast_name, "\n\n")
cat("Saved PDF:", output_pdf, "\n")
cat("Saved PNG:", output_png, "\n")
