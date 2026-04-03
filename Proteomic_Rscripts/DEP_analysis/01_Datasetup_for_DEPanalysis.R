library(dplyr)
library(tidyr)
library(stringr)
library(DEP)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(paletteer)


##Get table from file

#Protmcs_HekatRALexps_v1<-read_csv("quant_DIANN_outputs/HekRDH12andGFP_combinedDatasets_atRAL5hr_with24hrRecvry_rawdata.csv")
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


#remove extra contaminatns keratin protiens as found from doi.org/10.1021/acs.jproteome.2c00145
contam_protein_list <- c(
  #From Keratin
  "Q6E0U4", "P20930", "Q5D862", "Q86YZ3", "P04264", "P13645", "Q99456", "P13646",
  "P02533", "P19012", "P08779", "Q04695", "P05783", "P08727", "P35908", "P35900",
  "Q8N1A0", "Q9C075", "Q2M2I5", "Q7Z3Z0", "Q7Z3Y9", "Q7Z3Y8", "Q7Z3Y7", "P12035",
  "Q15323", "Q14532", "O76009", "Q14525", "O76011", "Q92764", "O76013", "O76014",
  "O76015", "Q6A163", "P19013", "Q6A162", "P13647", "P02538", "P04259", "P48668",
  "Q3KNV1", "P08729", "Q3SY84", "Q14CN4", "Q86Y46", "Q7RTS7", "O95678", "Q01546",
  "Q7Z794", "Q8N1N4", "Q5XKE5", "P05787", "Q6KB66", "Q14533", "Q9NSB4", "P78385",
  "Q9NSB2", "P78386", "O43790", "A6NCN2", "P35527", "P60331", "P60014", "P60412",
  "P60413", "P60368", "P60369", "P60372", "P60370", "P60371", "P60409", "P60410",
  "P60411", "Q07627", "Q8IUC1", "P59990", "P59991", "P60328", "P60329", "Q8IUG1",
  "Q8IUC0", "Q52LG2", "Q3SY46", "Q3LI77", "P0C5Y4", "Q9BYS1", "Q3LI76", "A8MUX0",
  "Q9BYP8", "Q8IUB9", "Q3LHN2", "Q7Z4W3", "Q3LI73", "Q3LI72", "Q3LI70", "Q3SYF9",
  "Q3LI54", "Q3LI63", "Q3LI61", "Q9BYU5", "Q3LI58", "Q3LI59", "Q3LHN1", "Q9BYT5",
  "Q3MIV0", "P0C7H8", "A1A580", "Q9BYR9", "Q3LI83", "Q3LHN0", "Q6PEX3", "Q3LI81",
  "A8MX34", "Q9BYR8", "Q9BYR7", "Q9BYR6", "Q9BYQ7", "Q9BYQ6", "Q9BQ66", "Q9BYR5",
  "Q9BYR4", "Q9BYR3", "Q9BYR2", "Q9BYQ5", "Q9BYR0", "Q9BYQ9", "Q9BYQ8", "Q6L8H4",
  "Q6L8G5", "Q6L8G4", "Q701N4", "Q6L8H2", "Q6L8H1", "Q701N2", "Q6L8G9", "Q6L8G8",
  "O75690", "P26371", "Q3LI64", "Q3LI66", "Q3LI67", "Q8IUC3", "Q8IUC2", "A8MXZ3",
  "Q9BYQ4", "Q9BYQ3", "Q9BYQ2", "A8MVA2", "A8MTY7", "Q9BYQ0", "Q9BYP9",
  #skin
  "P81605","Q08554","Q15517","P31151","Q9NZT1","Q96DR8",
  #From FBS
  "P01023","P02774","P00734","P01024","P02771","P04114","Q06033","P25311", "P00747"
  
)

abundance_table <- abundance_table %>%
  filter(!Protein.Group %in% contam_protein_list)
saveRDS(Hek_norm, file = "RData/Hekcells_Prep.rds")
write.csv(abundance_table, "Proteomic_output_txts/RawData_filtered_results.csv")

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

#saveRDS(data_se, file = "RData/DEP_SEtable_initial.rds")

# what metadata columns exist?
colnames(as.data.frame(colData(data_se)))

# look at the first rows
#head(as.data.frame(colData(data_imp)))

#print(as.data.frame(colData(data_imp)))

plot_frequency(data_se)
ggsave(filename = file.path(path_figures,"ProteinFrequency_acrossSamps.png"), width = 6, height = 5)

# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)

plot_detect(data_filt)

plot_missval(data_filt)

plot_normalization(data_filt)



