# =========================================================
# Bar plot of mean intensity +/- SD for one protein
# using df_long sample-level intensities
# =========================================================

library(dplyr)
library(ggplot2)
library(stringr)

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------
df_long<-read.csv("Proteomic_output_txts/Data_imputed_results.csv")

protein_gene   <- "TUBA4A"     # gene symbol to plot
protein_uniprot <- NULL      # optional, e.g. "P36969"
plot_points    <- TRUE       # show replicate points
use_sd         <- TRUE       # TRUE = SD, FALSE = SEM

# Set the order you want conditions displayed
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

# Pretty x-axis labels to match desired display
condition_labels <- c(
  "RDH12_control_atRAL5hr"            = "Acute RDH12 Veh",
  "RDH12_100_atRAL5hr"                = "Acute RDH12 100uM",
  "RDH12_200_atRAL5hr"                = "Acute RDH12 200uM",
  "RDH12_control_atRAL5hr.24h_recvr"  = "Recvr RDH12 Veh",
  "RDH12_100_atRAL5hr.24h_recvr"      = "Recvr RDH12 100uM",
  "RDH12_200_atRAL5hr.24h_recvr"      = "Recvr RDH12 200uM",
  "GFP_control_atRAL5hr.24h_recvr"    = "Recvr WT Veh",
  "GFP_100_atRAL5hr.24h_recvr"        = "Recvr WT 100uM",
  "GFP_200_atRAL5hr.24h_recvr"        = "Recvr WT 200uM"
)

# Optional fill colors by condition
condition_fill <- c(
  "RDH12_control_atRAL5hr"            = "#BDBDBD",
  "RDH12_100_atRAL5hr"                = "#E08A5F",
  "RDH12_200_atRAL5hr"                = "#E6CDB9",
  "RDH12_control_atRAL5hr.24h_recvr"  = "#BDBDBD",
  "RDH12_100_atRAL5hr.24h_recvr"      = "#E08A5F",
  "RDH12_200_atRAL5hr.24h_recvr"      = "#E6CDB9",
  "GFP_control_atRAL5hr.24h_recvr"    = "#74A2C6",
  "GFP_100_atRAL5hr.24h_recvr"        = "#74A2C6",
  "GFP_200_atRAL5hr.24h_recvr"        = "#74A2C6"
)

# output file
out_file <- paste0("Proteomic_Figs/Intensity_barplot_", protein_gene, ".svg")

# ---------------------------------------------------------
# 2. FILTER TO PROTEIN OF INTEREST
# ---------------------------------------------------------

protein_df <- df_long %>%
  filter(
    if (!is.null(protein_uniprot)) ID == protein_uniprot else TRUE,
    name == protein_gene | Genes == protein_gene
  )

if (nrow(protein_df) == 0) {
  stop("No rows found for the requested protein. Check 'protein_gene' or 'protein_uniprot'.")
}

# If multiple IDs map to same gene symbol, optionally inspect:
protein_df %>%
  distinct(name, Genes, ID, Protein.Group) %>%
  print()

# ---------------------------------------------------------
# 3. SUMMARIZE MEAN AND SD/SEM BY CONDITION
# ---------------------------------------------------------

summary_df <- protein_df %>%
  mutate(
    condition = factor(condition, levels = condition_order)
  ) %>%
  filter(!is.na(condition)) %>%
  group_by(condition) %>%
  summarise(
    mean_intensity = mean(intensity, na.rm = TRUE),
    sd_intensity   = sd(intensity, na.rm = TRUE),
    n              = sum(!is.na(intensity)),
    sem_intensity  = sd_intensity / sqrt(n),
    .groups = "drop"
  )

# FIX HERE
summary_df <- summary_df %>%
  mutate(
    error_value = if (use_sd) sd_intensity else sem_intensity
  )

print(summary_df)

# ---------------------------------------------------------
# 4. OPTIONAL: REPLICATE-LEVEL TABLE FOR PLOTTING POINTS
# ---------------------------------------------------------

plot_points_df <- protein_df %>%
  mutate(
    condition = factor(condition, levels = condition_order)
  ) %>%
  filter(!is.na(condition))

# ---------------------------------------------------------
# 5. BUILD TITLES
# ---------------------------------------------------------

gene_for_title <- protein_df %>%
  pull(name) %>%
  na.omit() %>%
  unique()

gene_for_title <- if (length(gene_for_title) > 0) gene_for_title[1] else protein_gene

id_for_title <- protein_df %>%
  pull(ID) %>%
  na.omit() %>%
  unique()

id_for_title <- if (length(id_for_title) > 0) id_for_title[1] else NA_character_

subtitle_text <- if (!is.na(id_for_title)) {
  paste("Mean intensity across proteomics samples | UniProt:", id_for_title)
} else {
  "Mean intensity across proteomics samples"
}

# ---------------------------------------------------------
# 6. PLOT
# ---------------------------------------------------------

baseline <- 20

summary_df2 <- summary_df %>%
  mutate(
    x = seq_len(n()),
    xmin = x - 0.35,
    xmax = x + 0.35
  )

plot_points_df2 <- plot_points_df %>%
  mutate(
    condition = factor(condition, levels = levels(summary_df$condition)),
    x = as.numeric(condition)
  )

p <- ggplot() +
  geom_rect(
    data = summary_df2,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = baseline,
      ymax = mean_intensity,
      fill = condition
    ),
    color = "black"
  ) +
  geom_errorbar(
    data = summary_df2,
    aes(
      x = x,
      ymin = pmax(mean_intensity - error_value, baseline),
      ymax = mean_intensity + error_value
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  {
    if (plot_points) geom_point(
      data = plot_points_df2,
      aes(x = x, y = intensity),
      size = 2.3,
      position = position_jitter(width = 0.08, height = 0)
    )
  } +
  scale_x_continuous(
    breaks = summary_df2$x,
    labels = condition_labels[as.character(summary_df2$condition)]
  ) +
  scale_y_continuous(limits = c(baseline, NA)) +
  scale_fill_manual(values = condition_fill, guide = "none") +
  labs(
    title = paste(gene_for_title, "Expression by Condition"),
    subtitle = subtitle_text,
    x = "Condition",
    y = "Mean log2 intensity"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 16)
  )

print(p)
# ---------------------------------------------------------
# 7. SAVE
# ---------------------------------------------------------

ggsave(
  filename = out_file,
  plot = p,
  width = 12,
  height = 7
)