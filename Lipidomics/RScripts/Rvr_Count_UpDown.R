# =========================================================
# Count up/down proteins per contrast and make stacked bar plot
# with user-defined plotting order
# =========================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------
#Recovery_DEA_results_v1_corrected.csv
DEPresults_v2 <- read.csv("Lipidomics/output_txts/Recovery_DEA_results_v1_corrected.csv")
colnames(DEPresults_v2)
padj_cutoff <- 0.1
lfc_cutoff  <- log2(1.3)

# Define labels in the EXACT order you want them plotted
pretty_contrast_names <- c(
  
  "RDH12_100_vs_RDH12_Veh" =
    "R:RDH12 100uM vs Veh",
  "RDH12_200_vs_RDH12_Veh" =
    "R:RDH12 200uM vs Veh",
  "RDH12_200_vs_RDH12_100" =
    "R:RDH12 200uM vs 100uM",
  
  "Control_100_vs_Control_Veh" =
    "R:WT 100uM vs Veh",
  "Control_200_vs_Control_100" =
    "R:WT 200uM vs 100uM",
  "Control_200_vs_Control_Veh" =
    "R:WT 200uM vs Veh",
  
  "RDH12_100_vs_Control_100" =
    "R:RDH12 100uM vs WT 100uM",
  "RDH12_200_vs_Control_200" =
    "R:RDH12 200uM vs WT 200uM",
  "RDH12_Veh_vs_Control_Veh" =
    "R:RDH12 Veh vs WT Veh"
 
)

# Output file
out_file <- "Lipidomics/Figures/Recovery/stacked_DE_counts_by_contrast.svg"

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------

get_all_contrasts <- function(df) {
  colnames(df) %>%
    .[str_detect(., "_ratio$")] %>%
    str_remove("_ratio$")
}

count_sig_by_contrast <- function(df, contrast, padj_cutoff = 0.1, lfc_cutoff = log2(1.5)) {
  ratio_col <- paste0(contrast, "_ratio")
  padj_col  <- paste0(contrast, "_p.val")
  
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
      ratio = .data[[ratio_col]],
      padj  = .data[[padj_col]]
    ) %>%
    mutate(
      direction = case_when(
        !is.na(padj) & !is.na(ratio) & padj <= padj_cutoff & ratio >=  lfc_cutoff ~ "Upregulated",
        !is.na(padj) & !is.na(ratio) & padj <= padj_cutoff & ratio <= -lfc_cutoff ~ "Downregulated",
        TRUE ~ "NS"
      )
    )
  
  tmp %>%
    filter(direction != "NS") %>%
    count(direction, name = "n") %>%
    complete(
      direction = c("Downregulated", "Upregulated"),
      fill = list(n = 0)
    ) %>%
    mutate(contrast = contrast)
}

# ---------------------------------------------------------
# 3. GET CONTRASTS AND CHECK MAPPINGS
# ---------------------------------------------------------

all_contrasts <- get_all_contrasts(DEPresults_v2)

# Check for contrasts in the data that are not in your desired order vector
missing_from_mapping <- setdiff(all_contrasts, names(pretty_contrast_names))
if (length(missing_from_mapping) > 0) {
  message("These contrasts are in DEPresults_v2 but not in pretty_contrast_names:")
  print(missing_from_mapping)
}

# Check for names in your mapping that are not in the data
missing_from_data <- setdiff(names(pretty_contrast_names), all_contrasts)
if (length(missing_from_data) > 0) {
  message("These contrasts are in pretty_contrast_names but not in DEPresults_v2:")
  print(missing_from_data)
}

# ---------------------------------------------------------
# 4. BUILD SUMMARY TABLE
# ---------------------------------------------------------

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

# Force the plotting order to follow pretty_contrast_names exactly
desired_order <- unname(pretty_contrast_names)

# If any contrasts were not mapped, append their raw names to the end
extra_labels <- setdiff(unique(count_table$contrast_label), desired_order)

count_table <- count_table %>%
  mutate(
    contrast_label = factor(
      contrast_label,
      levels = c(desired_order, extra_labels)
    ),
    direction = factor(direction, levels = c("Upregulated", "Downregulated"))
  )

print(count_table)

# ---------------------------------------------------------
# 5. PLOT
# ---------------------------------------------------------
p <- ggplot(count_table, aes(x = contrast_label, y = n, fill = direction)) +
  geom_col(color = "black", width = 0.8) +
  scale_fill_manual(
    values = c(
      "Downregulated" = "#00B8B8",
      "Upregulated"   = "purple"
    )
  ) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10)
  )+
  labs(
    title = "Lipids with Differential Abundance",
    subtitle =  "Fold Change >=1.3 (pval<=0.1)",
    x = NULL,
    y = "Number of Lipids",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top"
  )
print(p)

# ---------------------------------------------------------
# 6. OPTIONAL SAVE
# ---------------------------------------------------------

ggsave(
  filename = "Lipidomics/Figures/Recovery/stacked_DEP_counts_by_contrast.svg",
  plot = p,
  width = 8,
  height = 7
)