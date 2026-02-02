library(DEP)
library(SummarizedExperiment)
library(dplyr)

# ------------------------------------------------------------
# 1) Start from your filtered SE (DO NOT mask 1/3!)
#    thr = 1 keeps proteins with >=2/3 in at least one condition
#    (you can keep thr=1 or change if you want)
# ------------------------------------------------------------
data_filt <- filter_missval(data_se, thr = 1)

mat  <- assay(data_filt)
meta <- as.data.frame(colData(data_filt))

stopifnot("condition" %in% colnames(meta))

cond_vec <- meta$condition

# ------------------------------------------------------------
# 2) Build a logical matrix of structural MNAR positions:
#    structural_na[row, col] = TRUE if that protein is 0/3 in
#    the sample's condition (i.e., all reps missing)
# ------------------------------------------------------------
structural_na <- matrix(FALSE, nrow(mat), ncol(mat), dimnames = dimnames(mat))

for (cond in unique(cond_vec)) {
  cols <- which(cond_vec == cond)
  
  # how many detections per protein inside this condition?
  det_counts <- rowSums(!is.na(mat[, cols, drop = FALSE]))
  
  # 0/3 in this condition => structural MNAR block
  is_zero_of_three <- det_counts == 0
  
  # mark ALL samples of this condition as structural NA for those proteins
  structural_na[is_zero_of_three, cols] <- TRUE
}

# ------------------------------------------------------------
# 3) Impute ALL other missing values (MAR) with kNN
#    This will fill missing values in 1/3 and 2/3 blocks
#    (and might also fill 0/3 blocks — we'll restore those next)
# ------------------------------------------------------------
data_imp_tmp <- impute(data_filt, fun = "knn")

# ------------------------------------------------------------
# 4) Restore 0/3 or 1/3 blocks back to NA forever
#replaces the 0/3 conditions back to NA
# ------------------------------------------------------------
mat_imp <- assay(data_imp_tmp)
mat_imp[structural_na] <- NA_real_

data_imp <- data_filt
assay(data_imp) <- mat_imp

# ------------------------------------------------------------
# 5) Quick sanity checks (pick any proteins you care about)
#    - 1/3 blocks should now be filled (except 0/3 blocks)
#    - 0/3 blocks should remain NA
# ------------------------------------------------------------
# BEFORE imputation
rdh12_raw <- assay(data_filt)[rownames(data_filt) == "RDH12", ]
rdh12_raw

rdh12_filt <- assay(data_filt)[rownames(data_filt) == "RDH12", ]
rdh12_filt
# AFTER imputation
rdh12_imp <- assay(data_imp)[rownames(data_imp) == "RDH12", ]
rdh12_imp

####Checking
# BEFORE imputation
RUSF1_raw <- assay(data_filt)[rownames(data_filt) == "RUSF1", ]
RUSF1_raw

# AFTER imputation
RUSF1_imp <- assay(data_imp)[rownames(data_imp) == "RUSF1", ]
RUSF1_imp

cat("Proteins kept for imputation:     ", nrow(data_imp), "\n")

# BEFORE imputation
RAB15_raw <- assay(data_filt)[rownames(data_filt) == "RAB15", ]
RAB15_raw

# AFTER imputation
RAB15_imp <- assay(data_imp)[rownames(data_imp) == "RAB15", ]
RAB15_imp
