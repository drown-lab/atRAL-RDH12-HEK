
library(dplyr)
library(readxl)
library(tidyverse)
library(writexl)
library(ggplot2)
library(stringr)
path_figures<- "Lipidomics/Figures/Acute/"
#Load metadata
acute_expdesign<-read.csv("Lipidomics/metadata/Sample_description guide_5hratral.csv")

###Quick Processing Script for Finding Candidates of Lipids Analysis
##STEP 1: Take .csv ouputs from script, and change headers in each one so they are the same before running this script. 
######for doing the combining of all data folders from MRM QqQ MPF data
# Define the initial folder path
#need to define where the output of the Perl script that produces the .csv files of the MRM data files
#his script will look in the first folder for sub folders called "Absolute_intensity.csv" which is the compiled output of the Perl script
#This script will take that csv file and add each csv into one bit data table called "combined_data"
#==============
#Step 1:  Define the initial folder path
initial_folder <- "C:/Users/LabUser/Desktop/StemCells/lipidomic/20250702_HekCells_atRALtreated_rams/Life Sciences Native LIpids"

###run the next set of code togeter
# Create a list of all subfolders within the initial folder
subfolders <- dir(initial_folder, recursive = TRUE, full.names = TRUE)

# Filter out files and only keep directories
subfolders <- subfolders[file.info(subfolders)$isdir == FALSE]

# Remove the initial folder from the list
subfolders <- subfolders[subfolders != initial_folder]

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each subfolder
for (subfolder in subfolders) {
  # Check if the subfolder contains the "Absolute_intensity.csv" file
  file_path <- file.path(dirname(subfolder), "Absolute_intensity.csv")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    print(paste("File found:", file_path))
    # Read the CSV file into a data frame
    data_frame <- read.csv(file_path)
    
    # Check if combined_data is empty
    if (nrow(combined_data) == 0) {
      combined_data <- data_frame
    } else {
      # Bind the data frame to the combined data
      combined_data <- rbind(combined_data, data_frame)
    }
  } else {
    print(paste("File not found:", file_path))
  }
}

##############################################
##STEP 2: Clean up naming
Acute_lipidsall_v1 <- combined_data |> 
  janitor::clean_names() 

Acute_lipidsall_v2 <- Acute_lipidsall_v1 |>
  mutate(lipid_name = str_replace_all(lipid_name, ",", "_"))|> ##replace commas with underscore in Lipid name
  mutate(lipid_class1 = str_extract(lipid_name, "PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD")) |> # add all desired abbreviations
  mutate(lipid_class = str_extract(lipid_name, "LPC|LPE|PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD")) |>
  filter(!str_detect(lipid_class, "STD"))

Acute_lipidsall_v3 <- Acute_lipidsall_v2 |>
  filter(!str_detect(lipid_name, "_QUAL"))|>
  filter(!str_detect(lipid_class, "STD"))|>
  filter(!str_detect(lipid_class, "FA"))

Acute_lipidsall_v3<- Acute_lipidsall_v3|>
  unique()





##Mutate table further with TGs
Acute_lipidsall_v4 <-Acute_lipidsall_v3 |>
  mutate(DG_TG_acyl_chains = str_extract(lipid_name, "C\\d{1,2}:\\d")) ##adds column for acyl chain number for DGs and TGs

Acute_lipidsall_v4 <- Acute_lipidsall_v4 |>
  mutate(
    NL_chain = case_when(
      str_detect(lipid_name, "NL\\s*\\d{1,2}:\\d") ~ str_extract(lipid_name, "NL\\s*\\d{1,2}:\\d"),      # e.g. "NL 18:0"
      str_detect(lipid_name, "]_\\d{1,2}:\\d")     ~ str_extract(lipid_name, "(?<=]_)(\\d{1,2}:\\d)"),    # e.g. "]_22:6"
      TRUE ~ NA_character_
    )
  )



#######################################################
##Step 2: Calculate Intensity/Blank and Filter
#calc the Intensity /Blank

#do max value across all samples so with downstream analysis we retain intensity values
Acute_lipidsall_v5 <- Acute_lipidsall_v4 |>
  rowwise() |>
  mutate(maxvalue = max(c(s1, s2, s3, s4, s5, s6, s7, s8, s9)))

Acute_lipidsall_v5<-Acute_lipidsall_v5|>
  mutate(
    max_divid_blank = maxvalue/blank
  )


##Plot Distribution of Lipids
ggplot(Acute_lipidsall_v5, aes(x= max_divid_blank, fill = lipid_class1))+
  geom_histogram(bins = 40)+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Acute: Distribution of MRMs with Signal-to-Blank",
    subtitle = "Dashed line: 1.3x "
  )+
  scale_y_continuous(trans = "log10")
  
##Plot Distribution of Lipids
ggplot(Acute_lipidsall_v5, aes(x= max_divid_blank, y= lipid_class1, fill = lipid_class1))+
  geom_boxplot()+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Acute: Distribution of MRMs with Signal-to-Blank",
    subtitle = "Dashed line: 1.3x "
  )
  


#Filter by 30% higher than the blank
#Standard cut-off 
Acute_lipidsall_v6 <- Acute_lipidsall_v5 |>
  filter(max_divid_blank >= 1.3)

Acute_lipidsall_v6<- Acute_lipidsall_v6|>
  unique()

#============================================================
###STEP 3: Clean up MRM name
#this will add a column called mrm1 which will only keep 760->184 values and removes all decimal places 
Acute_lipidsall_v6 <- Acute_lipidsall_v6  |>
  mutate(mrm1 = str_replace(mrm, "(\\d+)\\.\\d+ -> (\\d+)\\.\\d+", "\\1 -> \\2"))



Acute_lipidsall_v7 <- Acute_lipidsall_v6 |>
  filter(!NL_chain %in% c("15:0"))

##Plot Distribution of Lipids
ggplot(Acute_lipidsall_v7, aes(x= max_divid_blank, y= lipid_class1, fill = lipid_class1))+
  geom_boxplot()+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Acute: Distribution of MRMs with Signal-to-Blank After filtering",
    subtitle = "Dashed line: 1.3x "
  )


# Find rows where 'mrm1' is duplicated
#should be all PC and PE because ran PEOx and PCOx methods
dups <- Acute_lipidsall_v6 %>%
  group_by(mrm1) %>%
  filter(n() > 1) %>%        # keep only duplicated groups
  arrange(mrm1)

dups2 <- Acute_lipidsall_v7 %>%
  group_by(mrm1) %>%
  filter(n() > 1) %>%        # keep only duplicated groups
  arrange(mrm1)

#remove duplicate values due to repeat technical inject but group by lipid_class to keep most info but with higher max/blank
Acute_lipidsall_v7v2 <- Acute_lipidsall_v7 |>
  group_by(lipid_class, mrm1) |>
  slice_max(max_divid_blank, n = 1, with_ties = FALSE) |>
  ungroup()
#========================================================================
###STEP 4: FIND # of IDs per lipid class or other filter
# Group the data by the desired variables and count the number of distinct MRMs

Acute_lipidsall_v7v2 <- Acute_lipidsall_v7v2 %>%
  mutate(
    precursor = as.numeric(str_split_fixed(mrm1, "\\s*->\\s*", 2)[,1]),
    product   = as.numeric(str_split_fixed(mrm1, "\\s*->\\s*", 2)[,2])
  )

dups3 <- Acute_lipidsall_v7v2 %>%
  group_by(mrm1) %>%
  filter(n() > 1) %>%        # keep only duplicated groups
  arrange(mrm1)

#============
#add more filtering to remove odd chain, and other non-sense identifications

df <- Acute_lipidsall_v7v2

#run script 1a_Filtering_Score_script.R
source("Lipidomics/RScripts/1a_Filtering_Score_script.R")

#this will output a parsed table: has a scoring value for lipid candidates when multiple options are present in lipid_name cell
#output-also separately break out cer_parsed table
#ouput-the "df_with_unsat_updated": final table output with the most likely candidate picked
#=============================
#Manual Filter

#Rename or manual filter out more based on score ouput
#First remove odd chains from chosen lipid classes
df_with_unsat_filtered <- df_with_unsat_updated %>%
  filter(!(lipid_class %in% c("CE", "Cer", "SM") & C_total %% 2 == 1))

##manually checked odd numbered lipids for false IDs
#found and remove or replace
#lipid_name
#remove DG 31:3 NL 20:0
#replace DG 33:0 NL 16:0 with TG 32:0 NL 16:0
#remove DG 33:2 NL 20:0
#replace DG 35:1 NL 16:1 with TG 34:1 NL 16:1
#remove 	DG 35:1 NL 20:0
#replace DG 37:7 NL 16:0 with DG 36:0 NL 16:0
#replace DG 37:7 NL 18:0 with DG 36:0 NL 18:0
#remove DG 39:6 NL 20:0
#remove DG 41:6 NL 16:1
#remove DG 41:6 NL 20:0
#remove DG 41:6 NL 22:6
#remove DG 41:5 NL 16:0
#remove DG 41:5 NL 22:5
#replace [TG37:0] NL 20:0 with DG 38:0 NL 20:0
#remove [TG37:0] NL 18:0
#remove [TG45:3] NL 16:0
#remove [TG50:9_TG49:2] NL 20:0
#remove [TG52:9_TG51:2] NL 16:0
#remove [TG52:9_TG51:2] NL 18:1
#remove DG 40:7 NL 20:0

df_with_unsat_filtered2 <- df_with_unsat_filtered %>%
  filter(
    !(lipid_name == "DG 31:3 NL 20:0"),
    !(lipid_name == "DG 33:1 NL 14:1"),
    
    !(lipid_name == "DG 33:2 NL 20:0"),
    !(lipid_name == "DG 35:1 NL 20:0"),
    !(lipid_name == "DG 39:6 NL 20:0"),
    !(lipid_name == "DG 41:6 NL 16:1"),
    !(lipid_name == "DG 41:6 NL 20:0"),
    !(lipid_name == "DG 41:6 NL 22:6"),
    !(lipid_name == "DG 41:5 NL 16:0"),
    !(lipid_name == "DG 41:5 NL 22:5"),
    !(lipid_name == "DG 40:7 NL 20:0"),
    !(lipid_name == "[TG37:0] NL 18:0"),
    !(lipid_name == "[TG45:3] NL 16:0"),
    !(lipid_name == "[TG50:9_TG49:2] NL 20:0"),
    !(lipid_name == "[TG52:9_TG51:2] NL 16:0"),
    !(lipid_name == "[TG52:9_TG51:2] NL 18:1")
  )

df_with_unsat_filtered2 <- df_with_unsat_filtered2 %>%
  mutate(
    lipid_name = case_when(
      lipid_name == "DG 33:0 NL 16:0" ~ "TG 32:0 NL 16:0",
      lipid_name == "DG 35:1 NL 16:1" ~ "TG 34:1 NL 16:1",
      lipid_name == "DG 37:7 NL 16:0" ~ "DG 36:0 NL 16:0",
      lipid_name == "DG 37:7 NL 18:0" ~ "DG 36:0 NL 18:0",
      lipid_name == "[TG37:0] NL 20:0" ~ "DG 38:0 NL 20:0",
      TRUE ~ lipid_name
    )
  )



#look aat SM and PCs identifiy isobars
SMandPCs <- df_with_unsat_updated %>%
  filter(lipid_class1 %in% c("PC", "SM"))
#Ids that need to be replace due to PC naming
#replace PC O-25:1;O2 with PC(16:0/8:0(COOH))
#replace PC(26:0)shift16 with PC 28:6_PC34:2OEP
#replace PC(27:0)_PC(O-28:0)shift16 with PC O-30:6_PC36:2OEP
#replace PC(29:0)_PC(O-30:0)shift16 with PC O-32:6_PC 38:2OEP
#replace PC(31:0)_PC(O-32:0)shift16 PC O-24:6_PC(22:0/8:0(COOH))
#replace PE(33:0)_PE(O-34:0)shift32 with PE 36:5
#replace PC 43:6 with PC O-44:6
#replace PC 43:2 with PC O-44:2_PC 44:9
#replace PC(30:0)_PC(O-31:0)shift32 with PC O-34:5
#replace PC(30:0)_PC(O-31:0)shift16 with PC 32:6
#replace LPC(16:0)_PC(O-16:0)_LPC(O-17:0)shift32 with LPC O-20:5
#replace PC(28:0)_PC(O-29:0)shift16 with PC(20:0/8:0(COOH))_PC(18:0/Aze)

df_with_unsat_filtered2 <- df_with_unsat_filtered2 %>%
  mutate(
    lipid_name = case_when(
      lipid_name == " PC O-25:1;O2 " ~ "PC(16:0/8:0(COOH))",
      lipid_name == "PC(26:0)shift16" ~ "PC 28:6_PC34:2OEP",
      lipid_name == "PC(27:0)_PC(O-28:0)shift16" ~ "PC O-30:6_PC36:2OEP",
      lipid_name == "PC(29:0)_PC(O-30:0)shift16" ~ "PC O-32:6_PC 38:2OEP",
      lipid_name == "PC(31:0)_PC(O-32:0)shift16" ~ "PC O-24:6_PC(22:0/8:0(COOH))",
      lipid_name == "PE(33:0)_PE(O-34:0)shift32" ~ "PE 36:5",
      lipid_name == "PC 43:6" ~ "PC O-44:6",
      lipid_name == "PC 43:2" ~ "PC O-44:2_PC 44:9",
      lipid_name == "PC(30:0)_PC(O-31:0)shift32" ~ "PC O-34:5",
      lipid_name == "PC(30:0)_PC(O-31:0)shift16" ~ "PC 32:6",
      lipid_name == "LPC(16:0)_PC(O-16:0)_LPC(O-17:0)shift32" ~ "LPC O-20:5",
      lipid_name == "PC(28:0)_PC(O-29:0)shift16" ~ "PC(20:0/8:0(COOH))_PC(18:0/Aze)",
     
      
        TRUE ~ lipid_name
    )
  )

df<-df_with_unsat_filtered2

df<-df|>
  select(-picked_candidate, -C_total, -DB_total, -NL_C, -NL_DB, -cand_score)

#run script 1a_Filtering_Score_script.R
source("Lipidomics/RScripts/1a_Filtering_Score_script.R")
#this reupdates the columns after replacement of likely IDs
#current table is titled df_with_unsat_updated

dups4 <- df_with_unsat_updated %>%
  group_by(mrm1) %>%
  filter(n() > 1) %>%        # keep only duplicated groups
  arrange(mrm1)

#still fix:
  #DG 29:1 NL 18:1 542-243 542.47847 remove
  #DG 37:7 NL 16:1 642*371 642.50977 remove 
  #DG 39:7 NL 20:0 670-341 670.54107 remove because duplicate
  #DG 39:0 NL 20:0 684-355  remove because duplicate

df_with_unsat_updated <- df_with_unsat_updated %>%
  filter(
    !(picked_candidate == "DG 29:1 NL 18:1"),
    !(picked_candidate == "DG 37:7 NL 16:1"),
    !(picked_candidate == "DG 39:7 NL 20:0"),
    !(picked_candidate == "DG 39:0 NL 20:0"),
    !(lipid_name == "LPS O-20:0;O"),
    
    
    )

df_with_unsat_updated <- df_with_unsat_updated %>%
  mutate(
    lipid_name = case_when(
  lipid_name == "PC O-40:7" ~ "PC 38:0",
  lipid_name == "PC O-38:7" ~ "PC 36:0",
  lipid_name == "PC O-42:7" ~ "PC 40:0",
  lipid_name == "DG O-38:9 NL 18:1" ~ "DG 36:2 NL 18:1",
  lipid_name == "DG O-38:8 NL 18:0" ~ "DG 36:1 NL 18:0",
  lipid_name == "DG O-38:8 NL 16:1" ~ "DG 36:1 NL 16:1",
  lipid_name == "DG O-38:8 NL 16:1" ~ "DG 36:1 NL 16:1",
  lipid_name == "DG O-40:9 NL 18:2" ~ "DG 38:2 NL 16:1",

  TRUE ~ lipid_name
    ))

df<-df_with_unsat_updated
df<-df|>
  select(-picked_candidate, -C_total, -DB_total, -NL_C, -NL_DB, -cand_score)
source("Lipidomics/RScripts/1a_Filtering_Score_script.R")


dups5 <- df_with_unsat_updated %>%
  group_by(picked_candidate) %>%
  filter(n() > 1) %>%        # keep only duplicated groups
  arrange(picked_candidate)


#===============
df_with_unsat_updated <- df_with_unsat_updated |>
  mutate(lipid_name = str_replace_all(picked_candidate, ",", "_"))|> ##replace commas with underscore in Lipid name
  mutate(lipid_class1 = str_extract(picked_candidate, "PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD")) |> # add all desired abbreviations
  mutate(lipid_class = str_extract(picked_candidate, "LPC|LPE|PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD"))
##Look at Summarization across data
Acute_filtered_summaryDGTG <- df_with_unsat_updated  %>%
  group_by( NL_chain,lipid_class) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

Acute_filtered_summary <- df_with_unsat_updated  %>%
  group_by( lipid_class) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

Acute_filtered_summary2 <- df_with_unsat_updated  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))
Acute_filtered_summary2 <- Acute_filtered_summary2 %>%
  mutate(
    percent = distinct_precursor / sum(distinct_precursor),
    label = paste0(lipid_class, " (", scales::percent(percent, accuracy = 1), ")")
  )

write.csv(df_with_unsat_updated, "Lipidomics/output_txts/Acute_filteredlipidtable.csv")
write.csv(Acute_filtered_summary2,"Lipidomics/output_txts/Acute_lipidclass_summary.csv")


#=====================
##Step 4B:
#pivot
Acute_lipidsall_v8<-df_with_unsat_updated|>
  select(c(picked_candidate, mrm,mrm1, precursor, product,lipid_class1, lipid_class, NL_chain,s1:s9))|>
  pivot_longer(
    cols = c(s1:s9),
    values_to = "absint",
    names_to = "sample"
  )


Acute_lipidsall_v9 <- Acute_lipidsall_v8 %>%
  left_join(acute_expdesign, by = "sample")

Acute_lipidsall_v9<-Acute_lipidsall_v9|>
  mutate(
    samplename = paste(GeneType, Treament, Rep, sep = "_")
  )

##Plot Distribution of Lipids
ggplot(Acute_lipidsall_v9, aes(x= log2(absint), y= lipid_class1, fill = lipid_class1))+
  geom_boxplot()+
  theme_bw(base_size=12)+
  labs(
    title = "Acute: Distribution of Intensity across Lipid Classes",
  )

ggplot(Acute_lipidsall_v9, aes(x= log2(absint), y= samplename, fill = samplename))+
  geom_boxplot()+
  theme_bw(base_size=12)+
  labs(
    title = "Acute: Distribution of Intensity across Samples",
  )


ggplot(Acute_lipidsall_v9, aes(x= log2(absint), y= samplename, fill = samplename))+
  geom_boxplot()+
  theme_bw(base_size=12)+
  labs(
    title = "Acute: Distribution of Intensity across Samples",
  )+
  facet_wrap(~lipid_class1, nrow= 5)


colnames(Acute_lipidsall_v9)

#====================================
#PCA
# optional, only if you use the paletteer scales below
# library(paletteer)

# 1) Wide matrix
dat_wide <- Acute_lipidsall_v9 %>%
  mutate(feature = paste(picked_candidate, mrm, sep = " | ")) %>%
  select(sample, GeneType, Treament, Rep, feature, absint) %>%
  summarise(absint = sum(absint, na.rm = TRUE),
            .by = c(sample, GeneType, Treament, Rep, feature)) %>%
  pivot_wider(names_from = feature, values_from = absint)

meta <- dat_wide %>% select(sample, GeneType, Treament, Rep)
X <- dat_wide %>% select(-sample, -GeneType, -Treament, -Rep)

X_log <- log1p(as.matrix(X))
X_log[is.na(X_log)] <- 0

# 2) PCA
pca <- prcomp(X_log, center = TRUE, scale. = TRUE)

scores <- as.data.frame(pca$x) %>%
  bind_cols(meta) %>%
  mutate(group = interaction(Treament, GeneType, drop = TRUE))

# variance explained
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

# 3) Helper: ellipse points
make_ellipse <- function(df, level = 0.68, n = 150, ridge = 1e-6) {
  if (nrow(df) < 3) return(NULL)
  
  mu <- colMeans(df[, c("PC1", "PC2")])
  S  <- stats::cov(df[, c("PC1", "PC2")]) + diag(ridge, 2)
  
  r2  <- stats::qchisq(level, df = 2)
  eig <- eigen(S)
  
  A <- eig$vectors %*% diag(sqrt(eig$values * r2)) %*% t(eig$vectors)
  
  tt <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(tt), sin(tt))
  
  ell <- sweep(circle %*% t(A), 2, mu, "+")
  tibble(PC1 = ell[, 1], PC2 = ell[, 2])
}

ellipse_df <- scores %>%
  group_by(group, Treament, GeneType) %>%
  group_modify(~ {
    e <- make_ellipse(.x, level = 0.68, n = 150, ridge = 1e-6)
    if (is.null(e)) return(tibble())
    e
  }) %>%
  ungroup()

# 4) Plot: filled ellipses + points
ggplot(scores, aes(PC1, PC2, color = Treament, shape = GeneType)) +
  geom_polygon(
    data = ellipse_df,
    aes(group = group, fill = Treament),
    alpha = 0.18,
    color = NA
  ) +
  geom_path(
    data = ellipse_df,
    aes(group = group, color = Treament),
    linewidth = 0.6
  ) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA",
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)"),
    color = "Treament",
    fill  = "Treament",
    shape = "GeneType"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
# If you want these palettes and have paletteer installed:
# + scale_color_paletteer_d("ggsci::category10_d3")
# + scale_fill_paletteer_d("ggsci::category10_d3")