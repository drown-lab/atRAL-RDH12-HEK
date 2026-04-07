library(dplyr)

# --------------------------------------------
# user-defined thresholds
# --------------------------------------------

data_results_w_missingclass<-read.csv("Lipidomics/output_txts/Acute_DEA_output_v1.csv")
alpha_cutoff <- 0.05
lfc_cutoff   <- log2(1.5)

# --------------------------------------------
# function to recompute BH-adjusted p-values
# and update significance columns
# --------------------------------------------
readjust_pvals_BH <- function(df, alpha = 0.05, lfc = log2(1.5)) {
  
  # find all raw p-value columns
  pval_cols <- grep("_p.val$", colnames(df), value = TRUE)
  
  for (p_col in pval_cols) {
    
    # matching column names
    contrast_base <- sub("_p.val$", "", p_col)
    padj_col  <- paste0(contrast_base, "_p.adj")
    ratio_col <- paste0(contrast_base, "_ratio")
    sig_col   <- paste0(contrast_base, "_significant")
    
    # extract p-values
    pvals <- df[[p_col]]
    
    # recompute BH, preserving NA positions
    new_padj <- rep(NA_real_, length(pvals))
    keep <- !is.na(pvals)
    new_padj[keep] <- p.adjust(pvals[keep], method = "BH")
    
    # overwrite or create p.adj column
    df[[padj_col]] <- new_padj
    
    # recompute significance if ratio column exists
    if (ratio_col %in% colnames(df)) {
      df[[sig_col]] <- !is.na(df[[padj_col]]) &
        !is.na(df[[ratio_col]]) &
        df[[padj_col]] <= alpha &
        abs(df[[ratio_col]]) >= lfc
    }
  }
  
  return(df)
}

# --------------------------------------------
# apply to your table
# --------------------------------------------
data_results_w_missingclass <- readjust_pvals_BH(
  df = data_results_w_missingclass,
  alpha = alpha_cutoff,
  lfc = lfc_cutoff
)

# --------------------------------------------
# optional: inspect updated p.adj columns
# --------------------------------------------
grep("_p.adj$", colnames(data_results_w_missingclass), value = TRUE)

# --------------------------------------------
# optional: save updated table
# --------------------------------------------
write.csv(
  data_results_w_missingclass,
  "Lipidomics/output_txts/Acute_DEA_results_v1_corrected.csv",
  row.names = FALSE
)