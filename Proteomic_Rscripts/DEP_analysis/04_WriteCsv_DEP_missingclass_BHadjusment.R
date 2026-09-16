library(dplyr)

# --------------------------------------------
# user-defined thresholds
#
# These drive only the `<contrast>_significant` columns below, which nothing reads and
# 05_make_SI_Table_Protein_DEA.R drops from the SI table. The Results-text counts use the
# stricter RESULTS_ALPHA / RESULTS_MIN_ABS_LOG2FC defined there, so these flags will not
# reproduce them.
# --------------------------------------------
alpha_cutoff <- 0.05
lfc_cutoff   <- log2(1.5)

# --------------------------------------------
# Quantification basis per contrast (added 12 Sep 2026)
#
# Each contrast "A_vs_B" is classified from the per-condition missingness classes
# (missclass_A / missclass_B; MNAR_0of3 = all three values QRILC-imputed):
#   measured      both conditions have >= 1 observed value
#   num_imputed   A entirely imputed (protein absent in A, present in B)
#   den_imputed   B entirely imputed (protein present in A, absent in B)
#   both_imputed  neither condition observed, so the test compares two QRILC draws.
#                 p.val and p.adj are set to NA before BH adjustment, removing them from
#                 testing and from the BH universe. The ratio is kept for transparency.
# One-side-imputed contrasts stay tested but are flagged: they are presence/absence calls,
# not fold changes.
# --------------------------------------------
basis_of <- function(num_class, den_class) {
  ni <- num_class == "MNAR_0of3"; di <- den_class == "MNAR_0of3"
  dplyr::case_when(is.na(num_class) | is.na(den_class) ~ NA_character_,
                   ni & di ~ "both_imputed", ni ~ "num_imputed", di ~ "den_imputed",
                   TRUE ~ "measured")
}

add_basis_and_mask <- function(df) {
  pval_cols <- grep("_p\\.val$", colnames(df), value = TRUE)
  summary_rows <- list()
  for (p_col in pval_cols) {
    contrast_base <- sub("_p\\.val$", "", p_col)
    parts <- strsplit(contrast_base, "_vs_", fixed = TRUE)[[1]]
    # A condition name containing "_vs_" would split into more than two parts and give
    # the wrong missclass lookup, leaving the contrast unmasked. Don't guess.
    if (length(parts) != 2) {
      warning("Cannot split contrast '", contrast_base, "' into exactly two conditions - basis not assigned"); next
    }
    num_col <- paste0("missclass_", parts[1]); den_col <- paste0("missclass_", parts[2])
    if (!all(c(num_col, den_col) %in% colnames(df))) {
      warning("No missclass columns for ", contrast_base, " - basis not assigned"); next
    }
    basis <- basis_of(df[[num_col]], df[[den_col]])
    df[[paste0(contrast_base, "_basis")]] <- basis
    both <- !is.na(basis) & basis == "both_imputed"
    n_removed_tested <- sum(both & !is.na(df[[p_col]]))
    df[[p_col]][both] <- NA_real_                 # removed from testing (p.adj recomputed below)
    summary_rows[[contrast_base]] <- data.frame(
      contrast = contrast_base,
      n_measured = sum(basis == "measured", na.rm = TRUE),
      n_num_imputed = sum(basis == "num_imputed", na.rm = TRUE),
      n_den_imputed = sum(basis == "den_imputed", na.rm = TRUE),
      n_both_imputed_removed = n_removed_tested)
  }
  attr(df, "basis_summary") <- do.call(rbind, summary_rows)
  df
}

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
data_results_w_missingclass <- add_basis_and_mask(data_results_w_missingclass)
basis_summary <- attr(data_results_w_missingclass, "basis_summary")
print(basis_summary)
write.csv(basis_summary, "Proteomic_output_txts/DEA_quantification_basis_summary.csv", row.names = FALSE)

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
  "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv",
  row.names = FALSE
)
