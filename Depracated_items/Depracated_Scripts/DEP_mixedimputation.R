data_se <- DEP::make_se(data_unique, columns, experimental_design_df)

data_se


plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)

plot_missval(data_filt)

nrow(data_filt)

nrow(data_filt)
assay(data_filt)[rownames(data_filt) == "RDH12", ]


#===================================================
# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(data_filt) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()


# Get a logical vector
MNAR <- names(data_filt) %in% proteins_MNAR
test=MNAR=="RDH12"
data_imp <- impute(data_filt, fun = "mixed",
                   randna = !MNAR,# we have to define MAR which is the opposite of MNAR
                   mar = "knn", mnar = "none")# imputation function for MAR and MNAR


plot_imputation(data_filt, data_knn, data_QRILC)

plot_normalization(data_filt, data_imp)
plot_imputation(data_se,data_filt, data_imp)
plot_detect(data_filt)
plot_detect(data_imp)

#========================================================
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
