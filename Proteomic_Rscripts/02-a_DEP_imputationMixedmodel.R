library(DEP)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
###Run code for mixed imputation

# ------------------------------------------------------------
# Mixed imputation:
#   - convert 0 -> NA
#   - MNAR_0of3 (within condition) -> QRILC
#   - all other missing -> KNN
# Returns:
#   $se_mixed   : imputed SummarizedExperiment
#   $status_tbl : protein × sample table with missing_class + imputation_used
# ------------------------------------------------------------
se<-data_filt

mixed_impute_se <- function(se, condition_col = "condition", seed = 1) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(colData(se)))
  
  set.seed(seed)
  
  mat0 <- assay(se)
  meta <- as.data.frame(colData(se))
  cond_vec <- as.character(meta[[condition_col]])
  
  # 0 -> NA so imputation treats zeros as missing
  mat <- mat0
  mat[is.finite(mat) & mat == 0] <- NA_real_
  
  se0 <- se
  assay(se0) <- mat
  
  # detected = non-NA after 0->NA
  detected_mat <- !is.na(mat)
  
  # ------------------------------------------------------------
  # MNAR_0of3 mask per condition: protein has 0 detections in that condition
  # ------------------------------------------------------------
  mnar0_mask <- matrix(FALSE, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  
  for (cond in unique(cond_vec)) {
    cols <- which(cond_vec == cond)
    det_counts <- rowSums(detected_mat[, cols, drop = FALSE])
    is_mnar0 <- det_counts == 0
    mnar0_mask[is_mnar0, cols] <- TRUE
  }
  
  # ------------------------------------------------------------
  # Two imputations
  # ------------------------------------------------------------
  se_knn   <- DEP::impute(se0, fun = "knn")
  se_qrilc <- DEP::impute(se0, fun = "QRILC")
  
  mat_knn   <- assay(se_knn)
  mat_qrilc <- assay(se_qrilc)
  
  # ------------------------------------------------------------
  # Merge: start from KNN, replace MNAR_0of3 missing positions with QRILC
  # ------------------------------------------------------------
  mat_mixed <- mat_knn
  replace_idx <- mnar0_mask & is.na(mat)   # only positions that were missing after 0->NA
  mat_mixed[replace_idx] <- mat_qrilc[replace_idx]
  
  se_mixed <- se0
  assay(se_mixed) <- mat_mixed
  
  # ------------------------------------------------------------
  # Status table: protein × sample (so you can merge back after DEP)
  # ------------------------------------------------------------
  status_tbl <- tibble(
    Protein.Group = rep(rownames(mat), times = ncol(mat)),
    sample_id     = rep(colnames(mat), each = nrow(mat)),
    condition     = rep(cond_vec, each = nrow(mat)),
    detected      = as.vector(detected_mat),
    MNAR_0of3     = as.vector(mnar0_mask)
  ) %>%
    mutate(
      missing_class = case_when(
        MNAR_0of3 ~ "MNAR_0of3",
        !detected ~ "MAR_partial",
        TRUE      ~ "Observed"
      ),
      imputation_used = case_when(
        MNAR_0of3 & !detected ~ "QRILC",
        !MNAR_0of3 & !detected ~ "KNN",
        TRUE ~ "None"
      )
    )
  
  list(se_mixed = se_mixed, status_tbl = status_tbl)
}

# ------------------------------------------------------------
# Run for Heavy and Light SEs
# ------------------------------------------------------------
out <- mixed_impute_se(data_filt, condition_col = "condition", seed = 1)
data_imp_mixed <- out$se_mixed
status <- out$status_tbl 


# Quick summaries
table(status$imputation_used)
table(status$missing_class)


#########===================================
#Check Specific Proteins
#MNAR RDH12
prot_idx <- which(rowData(data_filt)$Protein.Group == "Q96NR8")
MNAR_raw <- assay(data_filt)[prot_idx, ]
MNAR_raw

MNAR_imp<- assay(data_imp_mixed)[prot_idx, ]
MNAR_imp

#MAR1of3 RUSF1
prot_idx <- which(rowData(data_filt)$Protein.Group == "Q96GQ5")

mar1_raw <- assay(data_filt)[prot_idx, ]
mar1_raw

mar1_imp<- assay(data_imp_mixed)[prot_idx, ]
mar1_imp

#MAR2of3 RAB15
prot_idx <- which(rowData(data_filt)$Protein.Group == "P59190")

mar2_raw <- assay(data_filt)[prot_idx, ]
mar2_raw

mar2_imp<- assay(data_imp_mixed)[prot_idx, ]
mar2_imp




#===============================
## Plot intensity density distributions before and after imputation
plot_imputation(data_filt, data_imp_mixed) 
  labs(title = " Imputation")
#ggsave("figures_separate_DEP/HeavyChannel_Imputation.png", width = 6, height = 5)

#ggsave("figures_separate_DEP/LightChannel_Imputation.png", width = 6, height = 5)


# Visualize distribution of intensities by boxplots for all samples 
plot_normalization(data_filt, data_imp_mixed)

