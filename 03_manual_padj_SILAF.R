families <- list(
  Color_effect1 = c("White_7_vs_Blue_7_p.val","White_10_vs_Blue_10_p.val"),
  Time_effect1= c("Blue_7_vs_Blue_10_p.val", "White_7_vs_White_10_p.val")
)


adjust_by_family <- function(df,
                             families,
                             pval_suffix = "_p.val",
                             out_suffix  = "_p.adj_family",
                             method      = "BH") {
  stopifnot(is.data.frame(df))
  for (fam_name in names(families)) {
    pcols <- families[[fam_name]]
    # keep only columns that exist
    present <- intersect(pcols, names(df))
    missing <- setdiff(pcols, names(df))
    if (length(missing) > 0) {
      warning(sprintf("[%s] missing columns: %s",
                      fam_name, paste(missing, collapse = ", ")))
    }
    if (length(present) == 0) next
    
    # pull p-values as a matrix (proteins x contrasts in this family)
    mat <- as.matrix(df[, present, drop = FALSE])
    vec <- as.numeric(unlist(mat))              # column-wise unlist
    keep <- is.finite(vec)
    
    adj_vec <- rep(NA_real_, length(vec))
    adj_vec[keep] <- p.adjust(vec[keep], method = method)
    
    # reshape back to matrix with same dims/colnames
    adj_mat <- matrix(adj_vec, nrow = nrow(mat), ncol = ncol(mat))
    colnames(adj_mat) <- sub(paste0(pval_suffix, "$"),
                             out_suffix,
                             present)
    
    # write adjusted columns back
    df[, colnames(adj_mat)] <- adj_mat
  }
  df
}


LghtStrs_Fly_sumHLLFQ_v1 <- adjust_by_family(
  df          = LghtStrs_Fly_sumHLLFQ,
  families    = families,
  pval_suffix = "_p.val",
  out_suffix  = "_p.adj_family",
  method      = "BH"
)
