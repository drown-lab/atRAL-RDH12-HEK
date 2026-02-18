## 0) Normalize labels: ensure exactly one "LFQ." prefix
Hek_norm <- Protmcs_HekatRALexps_v1


Hek_norm <- Hek_norm %>%
  mutate(
    label = if_else(
      str_starts(Rundscrp, fixed("LFQ.")),
      Rundscrp,
      paste0("LFQ.", Rundscrp)
    )
  )
## 2) Experimental design: make condition UNIQUE across biology
experimental_design <- Hek_norm %>%
  distinct(label, Treatment, Rep, Genetype, ExpType) %>%
  transmute(
    label,
    # make condition include at least Genetype + Treatment; add ExpType if needed
    condition  = paste(Genetype, Treatment,ExpType, sep = "_"),
    # if you also need time/state to avoid any collision, use:
    # condition  = paste(Genetype, ExpType, Treatment, sep = "_"),
    replicate  = Rep,           # "a/b/c" is fine
    Genetype   = Genetype,
    ExpType    = ExpType,
    Treatment  = Treatment
  )

## 1) Build abundance_table (wide LFQ matrix for DEP)

abundance_table <- Hek_norm %>%
  select(Protein.Group, Genes, label, PG.MaxLFQ) %>%
  distinct() %>%
  pivot_wider(
    names_from  = label,
    values_from = PG.MaxLFQ
  )


# sanity checks
stopifnot(all(experimental_design$label %in% names(abundance_table)))
stopifnot(anyDuplicated(experimental_design$label) == 0)

## 3) Make required 'name' and 'ID'
data_unique <- DEP::make_unique(
  abundance_table,
  names = "Genes",
  ids   = "Protein.Group",
  delim = ";"
)

## 4) Align LFQ columns and coerce to numeric
columns <- match(experimental_design$label, names(data_unique))
stopifnot(all(names(data_unique)[columns] == experimental_design$label))
data_unique[columns] <- lapply(data_unique[columns], \(x) suppressWarnings(as.numeric(x)))

## 5) Convert design to data.frame and ensure UNIQUE rownames that DEP derives internally
# DEP internally uses condition/replicate to make rownames.
# Make sure those pairs are unique:
stopifnot(anyDuplicated(paste(experimental_design$condition, experimental_design$replicate, sep = "_")) == 0)

experimental_design_df <- as.data.frame(experimental_design)
# Set rownames to label (DEP will still be happy as long as condition/replicate pairs are unique)
rownames(experimental_design_df) <- experimental_design_df$label

## 6) Build SE (positional args only)
data_se <- DEP::make_se(data_unique, columns, experimental_design_df)

data_se


plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
#thr = 1 means "allow 1 missing per condition"
data_filt <- filter_missval(data_se, thr = 1)

#data_filt <- filter_missval(data_se, thr = 2)
plot_detect(data_filt)

plot_missval(data_filt)

plot_normalization(data_filt)

#data_norm<- normalize_vsn(data_filt)
#data_norm <- normalize(data_filt, method = "median")
#plot_normalization(data_norm, data_filt, data_imp)


## 2) Drop proteins that are ALL-NA in ANY condition
##    (prevents imputing structural absence, e.g., RDH12 in GFP)

mat  <- assay(data_filt)
meta <- as.data.frame(colData(data_filt))

# Ensure condition column exists
stopifnot("condition" %in% colnames(meta))

conditions <- unique(meta$condition)

# For each condition, flag proteins that are all NA across its replicates
all_na_any_condition <- sapply(conditions, function(cond) {
  cols <- rownames(meta)[meta$condition == cond]
  rowSums(is.na(mat[, cols, drop = FALSE])) == length(cols)
})

# If a protein is all-NA in at least one condition -> drop it
drop_idx <- apply(all_na_any_condition, 1, any)
data_filt2 <- data_filt[!drop_idx, ]

cat("Proteins in initial file:     ", nrow(data_se), "\n")
cat("Proteins after thr=1 filter:      ", nrow(data_filt),  "\n")
cat("Proteins dropped (all-NA in a condition): ", sum(drop_idx), "\n")
cat("Proteins kept for imputation:     ", nrow(data_filt2), "\n")

## ---------------------------------------------------------
data_imp <- impute(data_filt2, fun = "QRILC")
plot_imputation(data_se,data_filt, data_filt2, data_imp)
data_imp <- impute(data_filt, fun = "knn")

##randana are missing at random
randna<- 

data_imp <- impute(naset, method = "mixed",
                   randna = randna,
                   mar = "knn", mnar = "none")

####Checking
# BEFORE imputation
rdh12_raw <- assay(data_filt)[rownames(data_filt) == "RDH12", ]
rdh12_raw

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

##################

data_diff <- test_diff(data_imp, type = "manual", 
                       test = c("RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
                                "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr",
                                "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
                                "RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr",
                                "RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr",
                                "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                                "GFP_200_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" ,
                                "GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr"
                       ))

dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1.5))
plot_volcano(dep, contrast = "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr", label_size = 2, add_names = TRUE)

data_results <- get_results(dep)

DEPallcontrasts_mixedimputation_v1<-data_results