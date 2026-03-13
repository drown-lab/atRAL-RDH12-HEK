
#requires table made from 1_initial_processing_MRMdata_script that produces Recovery_lipidsall_v9

experimental_design<- read.csv("Lipidomics/metadata/Sample_description guide_5hratral.csv")
colnames(experimental_design)



data_abundance <- Acute_lipidsall_v9 %>%
  select(picked_candidate, mrm, sample, absint) %>%
  pivot_wider(
    names_from = sample,
    values_from = absint,
    names_prefix = "LFQ."
  )

data_unique <- make_unique(data_abundance, "picked_candidate", "mrm", delim = ";")

experimental_design2 <- experimental_design %>%
  transmute(
    label = paste0("LFQ.", sample),
    condition = paste(GeneType, Treament, sep = "_"),
    replicate = Rep
  )

LFQ_columns <- grep("LFQ\\.", colnames(data_unique))

data_se <- make_se(data_unique, LFQ_columns, experimental_design2)
data_se

setdiff(experimental_design2$label, colnames(data_unique))
setdiff(colnames(data_unique)[LFQ_columns], experimental_design2$label)


plot_frequency(data_se)

summary(data_unique[, LFQ_columns])
plot_normalization(data_se)

plot_imputation(data_se)
meanSdPlot(data_se)



data_normvsn<- normalize_vsn(data_se)
meanSdPlot(data_normvsn)
plot_normalization(data_se,data_normvsn)
plot_imputation(data_se, data_normvsn)

plot_missval(data_normvsn)
data_imp <- impute(data_normvsn, fun = "QRILC")

plot_imputation(data_se, data_normvsn, data_imp)
plot_normalization(data_se,data_normvsn, data_imp)
meanSdPlot(data_imp)


table(colData(data_normvsn)$condition)
colData(data_normvsn)$condition
#get correlation between mean and std dev
mat <- assay(data_se)

means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)

cor(means, sds)

mat <- assay(data_normvsn)

means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)

cor(means, sds)

mat <- assay(data_imp)

means <- rowMeans(mat, na.rm = TRUE)
sds <- apply(mat, 1, sd, na.rm = TRUE)

cor(means, sds)


test<-get_df_long(data_imp)
testnorm<-get_df_long(data_normvsn)