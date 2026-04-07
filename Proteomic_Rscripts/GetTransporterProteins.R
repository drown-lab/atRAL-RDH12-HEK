library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

#-----------------------------
# Load data
#-----------------------------
df <- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")

padj_cutoff <- 0.05
lfc_cutoff  <- log2(1.5)

#-----------------------------
# 1) Identify contrasts
#-----------------------------
pval_cols <- colnames(df)[str_detect(colnames(df), "_p\\.val$")]
contrasts <- str_remove(pval_cols, "_p\\.val$")

#-----------------------------
# 2) Define transporter genes
#-----------------------------
df <- df %>%
  mutate(
    is_transporter = str_detect(Gene, "^(SLC|ABC)") |
      str_detect(Gene, "ATP") |
      Gene %in% c("TFRC", "SLC39A14", "SLC7A11")
  )

#-----------------------------
# 3) Extract significant transporters per contrast
#-----------------------------
get_transporters <- function(ct) {
  
  pval_col  <- paste0(ct, "_p.val")
  padj_col  <- paste0(ct, "_p.adj")
  ratio_col <- paste0(ct, "_ratio")
  
  df %>%
    filter(is_transporter) %>%
    transmute(
      Gene,
      contrast = ct,
      ratio = .data[[ratio_col]],
      p.adj = .data[[padj_col]]
    ) %>%
    filter(!is.na(p.adj)) %>%
    filter(p.adj <= padj_cutoff) %>%
    mutate(
      direction = case_when(
        ratio >=  lfc_cutoff  ~ "up",
        ratio <= -lfc_cutoff  ~ "down",
        TRUE ~ "ns"
      )
    ) %>%
    filter(direction != "ns")
}

#-----------------------------
# 4) Combine across contrasts
#-----------------------------
transporters_all <- map_dfr(contrasts, get_transporters)

#-----------------------------
# 5) Summary table
#-----------------------------
transporters_summary <- transporters_all %>%
  arrange(contrast, direction, desc(abs(ratio)))

# View results
transporters_summary

write.csv(transporters_summary, "Proteomic_output_txts/Transporterproteins_signficant.csv")