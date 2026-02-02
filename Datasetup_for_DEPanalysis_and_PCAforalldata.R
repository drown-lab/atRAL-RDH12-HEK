library(dplyr)
library(tidyr)
library(stringr)
library(DEP)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(paletteer)


##Get table from file
Protmcs_HekatRALexps_v1<-read_excel("Z:/data/Projects/Rams_Collab_RDH12experiments/atRAL_experiments/proteomics/quant_data_filtered_DIANNandExcels/HekRDH12andGFP_combinedDatasets_atRAL5hr_with24hrRecvry_rawdata.xlsx")
colnames(Protmcs_HekatRALexps_v1)

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


# what metadata columns exist?
colnames(as.data.frame(colData(data_se)))

# look at the first rows
head(as.data.frame(colData(data_imp)))

#print(as.data.frame(colData(data_imp)))

plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)

plot_detect(data_filt)

plot_missval(data_filt)

plot_normalization(data_filt)

data_imp <- impute(data_filt, fun = "QRILC")
plot_imputation(data_filt, data_imp)

################
#STEP: Plot PCA
##########
## =========================================================
## PCA from DEP SummarizedExperiment (data_imp)
## - Adds % variance to PC1/PC2 axis labels
## - Uses a discrete color palette
## =========================================================


## 1) Extract expression matrix (imputed SE)
mat <- assay(data_imp)   # if you have named assays, you can use: assay(data_imp, "your_assay_name")

## 2) Run PCA (samples must be rows, so transpose)
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

## 3) Compute variance explained (%)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_var <- round(100 * var_explained[1], 1)
pc2_var <- round(100 * var_explained[2], 1)

## 4) Merge PCA scores with sample metadata
meta <- as.data.frame(colData(data_imp))

scores <- as.data.frame(pca$x) |>
  rownames_to_column("sample") |>
  left_join(
    meta |> rownames_to_column("sample"),
    by = "sample"
  )

## 5) Ensure 'condition' is discrete (factor)
scores <- scores |>
  mutate(condition = as.factor(condition))

## 6) Plot PC1 vs PC2 with discrete palette
p <- ggplot(scores, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA of Proteomics Data",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Condition"
  ) +
  scale_color_paletteer_d("ggsci::category20_d3")

print(p)

# what metadata columns exist?
colnames(as.data.frame(colData(data_imp)))

# look at the first rows
head(as.data.frame(colData(data_imp)))

#######################
###Plot PCA of only acute AtRAL experiments
keep <- meta$ExpType == "atRAL5hr"
se_atRAL5hr <- data_imp[, keep]

mat <- assay(se_atRAL5hr)
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x) |>
  rownames_to_column("sample") |>
  left_join(
    as.data.frame(colData(se_atRAL5hr)) |> rownames_to_column("sample"),
    by = "sample"
  ) |>
  mutate(condition = factor(condition))

ggplot(scores, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCA: atRAL5hr only",
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)")
  )


#######

keep <- meta$ExpType == "atRAL5hr"
se_atRAL5hr <- data_imp[, keep]

mat <- assay(se_atRAL5hr)
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x) |>
  rownames_to_column("sample") |>
  left_join(
    as.data.frame(colData(se_atRAL5hr)) |> rownames_to_column("sample"),
    by = "sample"
  ) |>
  mutate(condition = factor(condition))

ggplot(scores, aes(PC1, PC2, color = condition, fill = condition)) +
  # shaded ellipses
  stat_ellipse(
    geom  = "polygon",
    type  = "t",
    level = 0.68,        # ~1 SD (good for overlap visualization)
    alpha = 0.2,
    color = NA
  ) +
  # ellipse outlines
  stat_ellipse(
    geom  = "path",
    type  = "t",
    linewidth = 1
  ) +
  # points
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA: atRAL5hr only",
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)"),
    color = "Condition",
    fill  = "Condition"
  )+
  theme(plot.title = element_text(hjust = 0.5))

