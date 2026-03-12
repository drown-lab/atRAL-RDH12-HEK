library(tidyverse)
library(DEP)
library(dplyr)
#This code provides the first steps for filtering data into separate channel and making summarized experiment tables after initial processing

# load prepared data or run 01_PrepareData.R
#SILAC_v1 <- readRDS("Cell_Perturbed/output_txts/SILACprbd_v7_Prep.rds")

SILAC_v1<-SILACprbd_v7v1
print(unique(SILACprbd_v7v1$ExperimentName))

#==============
#Made se tables for heavy and light

data_heavy <- SILAC_v1 |>
  select(c(Protein.Group, Genes, ExperimentName, H)) |>
  pivot_wider(names_from = ExperimentName, values_from = H, names_prefix = "LFQ.")

data_heavy$Genes |> duplicated() |> any()
data_heavy |> group_by(Genes) |> summarize(frequency = n()) |> 
  arrange(desc(frequency)) |> filter(frequency > 1)
data_unique <- make_unique(data_heavy, "Genes", "Protein.Group", delim = ";")

LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_H_se <- make_se(data_unique, LFQ_columns, experimental_design)

data_light <- SILAC_v1 |>
  select(c(Protein.Group, Genes, ExperimentName, L)) |>
  pivot_wider(names_from = ExperimentName, values_from = L, names_prefix = "LFQ.")

data_light$Genes |> duplicated() |> any()
data_light |> group_by(Genes) |> summarize(frequency = n()) |> 
  arrange(desc(frequency)) |> filter(frequency > 1)
data_unique <- make_unique(data_light, "Genes", "Protein.Group", delim = ";")

LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_L_se <- make_se(data_unique, LFQ_columns, experimental_design)

colnames(colData(data_H_se))
as.data.frame(colData(data_H_se))|>head()

# Inspect available conditions first (recommended)
unique(colData(data_H_se)$condition)

#======================================
plot_frequency(data_H_se) +
  labs(title = "Heavy Channel - Protein identifications overlap")
#ggsave("Cell_Perturbed/Paper1/Figures/HeavyChannel_DataCompleteness.png", width = 6, height = 5)
ggsave(filename = file.path(path_figures,"HeavyChannel_DataCompleteness.png"), width = 7, height = 5)

plot_frequency(data_L_se) +
  labs(title = "Light Channel - Protein identifications overlap")
#ggsave("Cell_Perturbed/Paper1/Figures/LightChannel_DataCompleteness.png", width = 6, height = 5)
ggsave(filename = file.path(path_figures,"LightChannel_DataCompleteness.png"), width = 7, height = 5)

# Filter for proteins that are identified in all replicates of at least one condition
#data_H_filt <- filter_missval(data_H_se, thr = 0)
#data_L_filt <- filter_missval(data_L_se, thr = 0)

# Filter for proteins that are identified in 2/3 replicates of at least one condition
data_H_filt <- filter_missval(data_H_se, thr = 1)
data_L_filt <- filter_missval(data_L_se, thr = 1)

plot_numbers(data_H_filt) +
  labs(title = "Heavy Channel - Proteins per sample")
#ggsave("Cell_Perturbed/Paper1/HeavyChannel_ProteinIdentifications.png", width = 6, height = 5)
ggsave(filename = file.path(path_figures,"HeavyChannel_ProteinIdentifications.png"), width = 7, height = 5)

plot_numbers(data_L_filt) +
  labs(title = "Light Channel - Proteins per sample")
#ggsave("Cell_Perturbed/Paper1/LightChannel_ProteinIdentifications.png", width = 6, height = 5)
ggsave(filename = file.path(path_figures,"LightChannel_ProteinIdentifications.png"), width = 7, height = 5)


# Visualize distribution of intensities by boxplots for all samples 
#plot_normalization(data_H_se, data_H_filt)
#ggsave("figures_separate_DEP/HeavyChannel_Normalization.png", width = 6, height = 5)
#plot_normalization(data_L_se, data_L_filt)
#ggsave("figures_separate_DEP/LightChannel_Normalization.png", width = 6, height = 5)
#===========================================================
as.data.frame(colData(data_H_filt))|>head()
as.data.frame(colData(data_L_filt))|>head()

#GO to imputation code


plot_detect(data_H_filt) 
plot_detect(data_L_filt) 

plot_missval(data_H_filt)
plot_missval(data_L_filt)


#write.csv(SILAC_v1, file = "Cell_Perturbed/Paper1/output_txts/Paper1_initialSILAC_v1.csv")

write.csv(SILAC_v1, path_data_prep_CSV)
