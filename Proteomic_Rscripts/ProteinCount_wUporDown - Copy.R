# =========================================================
# Stacked DE counts by contrast
# - for NON-channel DEP output
# - uses existing *_ratio and *_p.adj columns
# - counts Upregulated / Downregulated proteins
# - DOES NOT round before filtering
# - rounds only for display if needed later
# =========================================================

# -------------------------
# 1. PACKAGES
# -------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# -------------------------
# 2. INPUTS
# -------------------------
DEPresults_v2 <- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")

padj_cutoff <- 0.05
lfc_cutoff  <- log2(1.5)

pretty_contrast_names <- c(
  "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr" =
    "A:RDH12 100uM vs Veh",
  "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr" =
    "A:RDH12 200uM vs 100uM",
  "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr" =
    "A:RDH12 200uM vs Veh",
  
  "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr" =
    "R:RDH12 100uM vs Veh",
  "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr" =
    "R:RDH12 200uM vs Veh",
  "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr.24h_recvr" =
    "R:RDH12 200uM vs 100uM",
  
  "GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" =
    "R:WT 100uM vs Veh",
  "GFP_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr" =
    "R:WT 200uM vs 100uM",
  "GFP_200_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" =
    "R:WT 200uM vs Veh",
  
  "RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr" =
    "R:RDH12 100uM vs WT 100uM",
  "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr" =
    "R:RDH12 200uM vs WT 200uM",
  
  "RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" =
    "R:RDH12 Veh vs WT Veh",
  "RDH12_control_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr" =
    "RvsA:RDH12 Veh vs Veh",
  "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr" =
    "RvsA:RDH12 100uM vs 100uM",
  "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_200_atRAL5hr" =
    "RvsA:RDH12 200uM vs 200uM"
)

out_file <- "Proteomic_Figs/stacked_DE_counts_by_contrast.svg"

# -------------------------
# 3. HELPER FUNCTIONS
# -------------------------

get_all_contrasts <- function(df) {
  colnames(df) %>%
    str_subset("_ratio$") %>%
    str_remove("_ratio$")
}

count_sig_by_contrast <- function(df, contrast, padj_cutoff = 0.05, lfc_cutoff = log2(1.5)) {
  
  ratio_col <- paste0(contrast, "_ratio")
  padj_col  <- paste0(contrast, "_p.adj")
  
  needed_cols <- c(ratio_col, padj_col)
  missing_cols <- setdiff(needed_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for contrast ", contrast, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  tmp <- df %>%
    transmute(
      ratio = as.numeric(.data[[ratio_col]]),
      padj  = as.numeric(.data[[padj_col]])
    ) %>%
    mutate(
      direction = case_when(
        !is.na(padj) & !is.na(ratio) &
          padj <= padj_cutoff & ratio >=  lfc_cutoff ~ "Upregulated",
        
        !is.na(padj) & !is.na(ratio) &
          padj <= padj_cutoff & ratio <= -lfc_cutoff ~ "Downregulated",
        
        TRUE ~ "NS"
      ),
      # display-only columns; not used for filtering
      ratio_display = round(ratio, 1),
      padj_display  = round(padj, 2)
    )
  
  tmp %>%
    filter(direction != "NS") %>%
    dplyr::count(direction, name = "n") %>%
    tidyr::complete(
      direction = c("Downregulated", "Upregulated"),
      fill = list(n = 0)
    ) %>%
    mutate(contrast = contrast)
}

# -------------------------
# 4. CHECK CONTRASTS
# -------------------------

all_contrasts <- get_all_contrasts(DEPresults_v2)

missing_from_mapping <- setdiff(all_contrasts, names(pretty_contrast_names))
if (length(missing_from_mapping) > 0) {
  message("These contrasts are in DEPresults_v2 but not in pretty_contrast_names:")
  print(missing_from_mapping)
}

missing_from_data <- setdiff(names(pretty_contrast_names), all_contrasts)
if (length(missing_from_data) > 0) {
  message("These contrasts are in pretty_contrast_names but not in DEPresults_v2:")
  print(missing_from_data)
}

# -------------------------
# 5. BUILD SUMMARY TABLE
# -------------------------

count_table <- lapply(
  all_contrasts,
  function(x) count_sig_by_contrast(
    df = DEPresults_v2,
    contrast = x,
    padj_cutoff = padj_cutoff,
    lfc_cutoff = lfc_cutoff
  )
) %>%
  bind_rows() %>%
  mutate(
    contrast_label = if_else(
      contrast %in% names(pretty_contrast_names),
      pretty_contrast_names[contrast],
      contrast
    )
  )

desired_order <- unname(pretty_contrast_names)
extra_labels <- setdiff(unique(count_table$contrast_label), desired_order)

count_table <- count_table %>%
  mutate(
    contrast_label = factor(
      contrast_label,
      levels = c(desired_order, extra_labels)
    ),
    direction = factor(direction, levels = c("Downregulated", "Upregulated"))
  )

# segment labels
count_table_labels <- count_table %>%
  mutate(label = if_else(n > 0, as.character(n), ""))

# totals above bars
count_totals <- count_table %>%
  group_by(contrast_label) %>%
  summarise(total = sum(n), .groups = "drop")

print(count_table_labels)
print(count_totals)

# -------------------------
# 6. PLOT
# -------------------------

p <- ggplot(count_table_labels, aes(x = contrast_label, y = n, fill = direction)) +
  geom_col(color = "black", width = 0.8) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4
  ) +
  geom_text(
    data = count_totals,
    aes(x = contrast_label, y = total, label = total),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Downregulated" = "#00B8B8",
      "Upregulated"   = "purple"
    )
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Proteins with Differential Expression",
    subtitle = "BH-adjusted p-value <= 0.05 and |log2FC| >= log2(1.5)",
    x = NULL,
    y = "Number of Proteins",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

print(p)

# -------------------------
# 7. SAVE
# -------------------------

ggsave(out_file, plot = p, width = 10, height = 7)

# -------------------------
# 8. OPTIONAL: SAVE SUMMARY TABLE
# -------------------------

write.csv(
  count_table_labels,
  "Proteomic_Figs/stacked_DE_counts_by_contrast_summary.csv",
  row.names = FALSE
)