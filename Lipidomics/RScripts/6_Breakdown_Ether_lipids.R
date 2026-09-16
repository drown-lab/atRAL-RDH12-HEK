library(dplyr)
library(tidyr)
library(stringr)

df<-read.csv("Lipidomics/output_txts/Recovery_filteredlipidtable.csv")
# Saturation class, used by the unsat summaries below (previously added by 5_Recovery_Breakdown_Lipid_Unsatlevels.R)
df <- df %>%
  mutate(
    unsat_class = case_when(
      DB_total == 0 ~ "SFA",
      DB_total == 1 ~ "MUFA",
      DB_total >= 2 ~ "PUFA",
      TRUE ~ NA_character_
    )
  )
df <- df %>%
  mutate(
    ether = str_extract(
      picked_candidate, "O-",
      
    )
  )
df <- df %>%
  mutate(lipid_class_ether = str_extract(picked_candidate,"PC O-|PE O-|PI O-|DG O-|PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM" )
  )
df<-df|>
  mutate(
  lipid_class_ether2 = paste(ether, lipid_class_ether, sep= "_"))

write.csv(df,"Lipidomics/output_txts/Recovery_filteredlipidtable2.csv" )


#Setup RDH12 Recovery tables
Recovery_RDH12_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s10, s11, s12,s13, s14, s15, s16,s17, s18)))

Recovery_RDH12_v2<-Recovery_RDH12_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_RDH12_v3<-Recovery_RDH12_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_RDH12_v4<-Recovery_RDH12_v3|>
  unique()

Recovery_RDH12__filtered_summary2 <- Recovery_RDH12_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_unsatfiltered_summary2 <- Recovery_RDH12_v4  %>%
  group_by( lipid_class_ether2, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))
#++====================
#Summarize GFP dataset
Recovery_GFP_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s1, s2, s3, s4, s5, s6,s7, s8, s9)))

Recovery_GFP_v2<-Recovery_GFP_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_GFP_v3<-Recovery_GFP_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_GFP_v4<-Recovery_GFP_v3|>
  unique()

Recovery_GFP__filtered_summary2 <- Recovery_GFP_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_GFP_unsatfiltered_summary2 <- Recovery_GFP_v4  %>%
  group_by( lipid_class_ether2, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))
#==============
#Setup RDH12 Recovery tables
Recovery_RDH12_100_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s13, s14, s15)))

Recovery_RDH12_100_v2<-Recovery_RDH12_100_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_RDH12_100_v3<-Recovery_RDH12_100_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_RDH12_100_v4<-Recovery_RDH12_100_v3|>
  unique()

Recovery_RDH12_200_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s16, s17, s18)))

Recovery_RDH12_200_v2<-Recovery_RDH12_200_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_RDH12_200_v3<-Recovery_RDH12_200_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_RDH12_200_v4<-Recovery_RDH12_200_v3|>
  unique()

Recovery_RDH12_veh_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s10, s11, s12)))

Recovery_RDH12_veh_v2<-Recovery_RDH12_veh_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_RDH12_veh_v3<-Recovery_RDH12_veh_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_RDH12_veh_v4<-Recovery_RDH12_veh_v3|>
  unique()

#######################################
#Setup up GFP recovery tables
Recovery_Ctrl_100_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s4, s5, s6)))

Recovery_Ctrl_100_v2<-Recovery_Ctrl_100_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_Ctrl_100_v3<-Recovery_Ctrl_100_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_Ctrl_100_v4<-Recovery_Ctrl_100_v3|>
  unique()

Recovery_Ctrl_200_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s7, s8, s9)))

Recovery_Ctrl_200_v2<-Recovery_Ctrl_200_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_Ctrl_200_v3<-Recovery_Ctrl_200_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_Ctrl_200_v4<-Recovery_Ctrl_200_v3|>
  unique()

Recovery_Ctrl_veh_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s1, s2, s3)))

Recovery_Ctrl_veh_v2<-Recovery_Ctrl_veh_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Recovery_Ctrl_veh_v3<-Recovery_Ctrl_veh_v2|>
  filter(max_divid_blank >= 1.3)

Recovery_Ctrl_veh_v4<-Recovery_Ctrl_veh_v3|>
  unique()
#=======================================
Recovery_RDH12_Veh_filtered_summary2 <- Recovery_RDH12_veh_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_Veh_filtered_summary2 <- Recovery_RDH12_Veh_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether2 = factor(lipid_class_ether2, levels = lipid_class_ether2))  


Recovery_RDH12_100_filtered_summary2 <- Recovery_RDH12_100_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_100_filtered_summary2 <- Recovery_RDH12_100_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether2 = factor(lipid_class_ether2, levels = lipid_class_ether2))  


Recovery_RDH12_200_filtered_summary2 <- Recovery_RDH12_200_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_200_filtered_summary2 <- Recovery_RDH12_200_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether2 = factor(lipid_class_ether2, levels = lipid_class_ether2))  


Recovery_Ctrl_veh_filtered_summary2 <- Recovery_Ctrl_veh_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_veh_filtered_summary2 <- Recovery_Ctrl_veh_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether2 = factor(lipid_class_ether2, levels = lipid_class_ether2)) 

Recovery_Ctrl_100_filtered_summary2 <- Recovery_Ctrl_100_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_100_filtered_summary2 <- Recovery_Ctrl_100_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether2 = factor(lipid_class_ether2, levels = lipid_class_ether2)) 

Recovery_Ctrl_200_filtered_summary2 <- Recovery_Ctrl_200_v4  %>%
  group_by( lipid_class_ether2) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_200_filtered_summary2 <- Recovery_Ctrl_200_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether2 = factor(lipid_class_ether2, levels = lipid_class_ether2)) 

#============================================
summarize_lipids <- function(df, label) {
  df %>%
    group_by(lipid_class_ether2) %>%
    summarise(distinct_precursor = n_distinct(precursor), .groups = "drop") %>%
    mutate(group = label)
}
combined_summary <- bind_rows(
  summarize_lipids(Recovery_RDH12_veh_v4,  "RDH12_Veh"),
  summarize_lipids(Recovery_RDH12_100_v4,  "RDH12_100"),
  summarize_lipids(Recovery_RDH12_200_v4,  "RDH12_200"),
  summarize_lipids(Recovery_Ctrl_veh_v4,   "Ctrl_Veh"),
  summarize_lipids(Recovery_Ctrl_100_v4,   "Ctrl_100"),
  summarize_lipids(Recovery_Ctrl_200_v4,   "Ctrl_200")
)

colnames(combined_summary)
combined_summary2<-combined_summary|>
  pivot_wider(
    names_from = "group",
    values_from = "distinct_precursor"
  )
write.csv(combined_summary2, "Lipidomics/output_txts/Recovery_lipidclassEther_eachTreatment_summary.csv")