
library(dplyr)
library(readxl)
library(tidyverse)
library(writexl)
library(ggplot2)
library(stringr)

#Load metadata
recovery_expdesign<-read.csv("Lipidomics/metadata/Sample_description guide_recovery.csv")

###Quick Processing Script for Finding Candidates of Lipids Analysis
##STEP 1: Take .csv ouputs from script, and change headers in each one so they are the same before running this script. 
######for doing the combining of all data folders from MRM QqQ MPF data
# Define the initial folder path
#need to define where the output of the Perl script that produces the .csv files of the MRM data files
#his script will look in the first folder for sub folders called "Absolute_intensity.csv" which is the compiled output of the Perl script
#This script will take that csv file and add each csv into one bit data table called "combined_data"
#==============
#Step 1:  Define the initial folder path
initial_folder <- "C:/Users/LabUser/Desktop/StemCells/lipidomic/Hek_5hrw24hrrecvr_Miranda18samples"

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
Recovery_lipidsall_v1 <- combined_data |> 
  janitor::clean_names() 

Recovery_lipidsall_v2 <- Recovery_lipidsall_v1 |>
  mutate(lipid_name = str_replace_all(lipid_name, ",", "_"))|> ##replace commas with underscore in Lipid name
  mutate(lipid_class1 = str_extract(lipid_name, "PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD")) |> # add all desired abbreviations
  mutate(lipid_class = str_extract(lipid_name, "LPC|LPE|PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD")) |>
  filter(!str_detect(lipid_class, "STD"))

Recovery_lipidsall_v3 <- Recovery_lipidsall_v2 |>
  filter(!str_detect(lipid_name, "_QUAL"))|>
  filter(!str_detect(lipid_class, "STD"))|>
  filter(!str_detect(lipid_class, "FA"))

Recovery_lipidsall_v3<- Recovery_lipidsall_v3|>
  unique()


Recovery_lipidsall_v1 <- Recovery_lipidsall_v1 |>
  mutate(lipid_name = str_replace_all(lipid_name, ",", "_"))|> ##replace commas with underscore in Lipid name
  mutate(lipid_class1 = str_extract(lipid_name, "PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD")) |> # add all desired abbreviations
  mutate(lipid_class = str_extract(lipid_name, "LPC|LPE|PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM|STD")) 
Recovery_lipidsall_v2_STDs<-Recovery_lipidsall_v1|>
  filter(lipid_class1=="STD")



##Mutate table further with TGs
Recovery_lipidsall_v4 <-Recovery_lipidsall_v3 |>
  mutate(DG_TG_acyl_chains = str_extract(lipid_name, "C\\d{1,2}:\\d")) ##adds column for acyl chain number for DGs and TGs

Recovery_lipidsall_v4 <- Recovery_lipidsall_v4 |>
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
Recovery_lipidsall_v5 <- Recovery_lipidsall_v4 |>
  rowwise() |>
  mutate(maxvalue = max(c(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18)))

Recovery_lipidsall_v5<-Recovery_lipidsall_v5|>
  mutate(
    max_divid_blank = maxvalue/blank
  )


##Plot Distribution of Lipids
ggplot(Recovery_lipidsall_v5, aes(x= max_divid_blank, fill = lipid_class1))+
  geom_histogram(bins = 40)+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Recovery: Distribution of MRMs with Signal-to-Blank",
    subtitle = "Dashed line: 1.3x "
  )+
  scale_y_continuous(trans = "log10")
  
##Plot Distribution of Lipids
ggplot(Recovery_lipidsall_v5, aes(x= max_divid_blank, y= lipid_class1, fill = lipid_class1))+
  geom_boxplot()+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Recovery: Distribution of MRMs with Signal-to-Blank",
    subtitle = "Dashed line: 1.3x "
  )
  


#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_lipidsall_v6 <- Recovery_lipidsall_v5 |>
  filter(max_divid_blank >= 1.3)

Recovery_lipidsall_v6<- Recovery_lipidsall_v6|>
  unique()

#============================================================
###STEP 3: Clean up MRM name
#this will add a column called mrm1 which will only keep 760->184 values and removes all decimal places 
Recovery_lipidsall_v6 <- Recovery_lipidsall_v6  |>
  mutate(mrm1 = str_replace(mrm, "(\\d+)\\.\\d+ -> (\\d+)\\.\\d+", "\\1 -> \\2"))

#remove duplicate values due to repeat technical inject but with higher max/blank
Recovery_lipidsall_v7 <- Recovery_lipidsall_v6 |>
  arrange(desc(max_divid_blank))|>
  distinct(mrm1, .keep_all = TRUE)

Recovery_lipidsall_v7 <- Recovery_lipidsall_v7 |>
  filter(!NL_chain %in% c("15:0"))

##Plot Distribution of Lipids
ggplot(Recovery_lipidsall_v7, aes(x= max_divid_blank, y= lipid_class1, fill = lipid_class1))+
  geom_boxplot()+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Recovery: Distribution of MRMs with Signal-to-Blank After filtering",
    subtitle = "Dashed line: 1.3x "
  )


# Find rows where 'mrm1' is duplicated
#should be all PC and PE because ran PEOx and PCOx methods
dups <- Recovery_lipidsall_v6 %>%
  group_by(mrm1) %>%
  filter(n() > 1) %>%        # keep only duplicated groups
  arrange(mrm1)

#========================================================================
###STEP 4: FIND # of IDs per lipid class or other filter
# Group the data by the desired variables and count the number of distinct MRMs

Recovery_lipidsall_v7 <- Recovery_lipidsall_v7 %>%
  mutate(
    precursor = as.numeric(str_split_fixed(mrm1, "\\s*->\\s*", 2)[,1]),
    product   = as.numeric(str_split_fixed(mrm1, "\\s*->\\s*", 2)[,2])
  )


#============
#add more filtering to remove odd chain, and other non-sense identifications






#===============
##Look at Summarization across data
Recovery_lipidsall_v7_summaryDGTG <- Recovery_lipidsall_v7  %>%
  group_by( NL_chain,lipid_class) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

Recovery_lipidsall_v7_summary <- Recovery_lipidsall_v7  %>%
  group_by( lipid_class) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

Recovery_lipidsall_v7_summary2 <- Recovery_lipidsall_v7  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

#=====================
##Step 4B:
#pivot
Recovery_lipidsall_v8<-Recovery_lipidsall_v7|>
  select(c(lipid_name, mrm,mrm1, precursor, product,lipid_class1, lipid_class, NL_chain,s10:s9))|>
  pivot_longer(
    cols = c(s10:s9),
    values_to = "absint",
    names_to = "sample"
  )


Recovery_lipidsall_v9 <- Recovery_lipidsall_v8 %>%
  left_join(recovery_expdesign, by = "sample")

Recovery_lipidsall_v9<-Recovery_lipidsall_v9|>
  mutate(
    samplename = paste(Genetype, Treatment, Rep, sep = "_")
  )

##Plot Distribution of Lipids
ggplot(Recovery_lipidsall_v9, aes(x= log2(absint), y= lipid_class1, fill = lipid_class1))+
  geom_boxplot()+
  theme_bw(base_size=12)+
  labs(
    title = "Recovery: Distribution of Intensity across Lipid Classes",
  )

ggplot(Recovery_lipidsall_v9, aes(x= log2(absint), y= samplename, fill = samplename))+
  geom_boxplot()+
  theme_bw(base_size=12)+
  labs(
    title = "Recovery: Distribution of Intensity across Samples",
  )


ggplot(Recovery_lipidsall_v9, aes(x= log2(absint), y= samplename, fill = samplename))+
  geom_boxplot()+
  theme_bw(base_size=12)+
  labs(
    title = "Recovery: Distribution of Intensity across Samples",
  )+
  facet_wrap(~lipid_class1, nrow= 5)


colnames(Recovery_lipidsall_v9)

#====================================
#PCA
# optional, only if you use the paletteer scales below
# library(paletteer)

# 1) Wide matrix
dat_wide <- Recovery_lipidsall_v9 %>%
  mutate(feature = paste(lipid_name, mrm, sep = " | ")) %>%
  select(sample, Genetype, Treatment, Rep, feature, absint) %>%
  summarise(absint = sum(absint, na.rm = TRUE),
            .by = c(sample, Genetype, Treatment, Rep, feature)) %>%
  pivot_wider(names_from = feature, values_from = absint)

meta <- dat_wide %>% select(sample, Genetype, Treatment, Rep)
X <- dat_wide %>% select(-sample, -Genetype, -Treatment, -Rep)

X_log <- log1p(as.matrix(X))
X_log[is.na(X_log)] <- 0

# 2) PCA
pca <- prcomp(X_log, center = TRUE, scale. = TRUE)

scores <- as.data.frame(pca$x) %>%
  bind_cols(meta) %>%
  mutate(group = interaction(Treatment, Genetype, drop = TRUE))

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
  group_by(group, Treatment, Genetype) %>%
  group_modify(~ {
    e <- make_ellipse(.x, level = 0.68, n = 150, ridge = 1e-6)
    if (is.null(e)) return(tibble())
    e
  }) %>%
  ungroup()

# 4) Plot: filled ellipses + points
ggplot(scores, aes(PC1, PC2, color = Treatment, shape = Genetype)) +
  geom_polygon(
    data = ellipse_df,
    aes(group = group, fill = Treatment),
    alpha = 0.18,
    color = NA
  ) +
  geom_path(
    data = ellipse_df,
    aes(group = group, color = Treatment),
    linewidth = 0.6
  ) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA",
    x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)"),
    color = "Treatment",
    fill  = "Treatment",
    shape = "Genetype"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
# If you want these palettes and have paletteer installed:
# + scale_color_paletteer_d("ggsci::category10_d3")
# + scale_fill_paletteer_d("ggsci::category10_d3")