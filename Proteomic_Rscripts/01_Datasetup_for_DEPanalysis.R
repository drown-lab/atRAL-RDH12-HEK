library(dplyr)
library(tidyr)
library(stringr)
library(DEP)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(paletteer)


##Get table from file
#Protmcs_HekatRALexps_v1<-read_excel("Z:/data/Projects/Rams_Collab_RDH12experiments/atRAL_experiments/proteomics/quant_data_filtered_DIANNandExcels/HekRDH12andGFP_combinedDatasets_atRAL5hr_with24hrRecvry_rawdata.xlsx")

Protmcs_HekatRALexps_v1<-read_excel("quant_DIANN_outputs/HekRDH12andGFP_combinedDatasets_atRAL5hr_with24hrRecvry_rawdata.xlsx")
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

saveRDS(Hek_norm, file = "RData/Hekcells_Prep.rds")

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

saveRDS(data_se, file = "RData/DEP_SEtable_initial.rds")

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



#data_imp <- impute(data_filt, fun = "QRILC")
#plot_imputation(data_filt, data_imp)

head(as.data.frame(colData(meta)))

head(meta)
