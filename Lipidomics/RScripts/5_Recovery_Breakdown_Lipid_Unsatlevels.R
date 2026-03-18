##Breakdown Acute data into Satand UnsatLevel by Treatment
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)
library(ggsci)

df<-read.csv("Lipidomics/output_txts/Recovery_filteredlipidtable.csv")
df <- df %>%
  mutate(
    unsat_class = case_when(
      DB_total == 0 ~ "SFA",
      DB_total == 1 ~ "MUFA",
      DB_total >= 2 ~ "PUFA",
      TRUE ~ NA_character_
    )
  )
#======================
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
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_unsatfiltered_summary2 <- Recovery_RDH12_v4  %>%
  group_by( lipid_class, unsat_class) %>%
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
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_GFP_unsatfiltered_summary2 <- Recovery_GFP_v4  %>%
  group_by( lipid_class, unsat_class) %>%
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
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_Veh_filtered_summary2 <- Recovery_RDH12_Veh_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class = factor(lipid_class, levels = lipid_class))  


Recovery_RDH12_100_filtered_summary2 <- Recovery_RDH12_100_v4  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_100_filtered_summary2 <- Recovery_RDH12_100_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class = factor(lipid_class, levels = lipid_class))  


Recovery_RDH12_200_filtered_summary2 <- Recovery_RDH12_200_v4  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_200_filtered_summary2 <- Recovery_RDH12_200_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class = factor(lipid_class, levels = lipid_class))  


Recovery_Ctrl_veh_filtered_summary2 <- Recovery_Ctrl_veh_v4  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_veh_filtered_summary2 <- Recovery_Ctrl_veh_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class = factor(lipid_class, levels = lipid_class)) 

Recovery_Ctrl_100_filtered_summary2 <- Recovery_Ctrl_100_v4  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_100_filtered_summary2 <- Recovery_Ctrl_100_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class = factor(lipid_class, levels = lipid_class)) 

Recovery_Ctrl_200_filtered_summary2 <- Recovery_Ctrl_200_v4  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_200_filtered_summary2 <- Recovery_Ctrl_200_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class = factor(lipid_class, levels = lipid_class)) 

#============================================
summarize_lipids <- function(df, label) {
  df %>%
    group_by(lipid_class) %>%
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
write.csv(combined_summary2, "Lipidomics/output_txts/Recovery_lipidclass_eachTreatment_summary.csv")
#=====================================================

#Complex Upset Plot
## 1) Build one combined table from the underlying precursor-level data
all_data <- bind_rows(
  Recovery_RDH12_veh_v4 %>% mutate(group = "RDH12_Veh"),
  Recovery_RDH12_100_v4 %>% mutate(group = "RDH12_100"),
  Recovery_RDH12_200_v4 %>% mutate(group = "RDH12_200"),
  Recovery_Ctrl_veh_v4  %>% mutate(group = "Ctrl_Veh"),
  Recovery_Ctrl_100_v4  %>% mutate(group = "Ctrl_100"),
  Recovery_Ctrl_200_v4  %>% mutate(group = "Ctrl_200")
) %>%
  select(precursor, lipid_class, group) %>%
  distinct()

## 2) Optional check: each precursor should ideally map to one lipid class
all_data %>%
  distinct(precursor, lipid_class) %>%
  dplyr::count(precursor) %>%
  filter(n > 1)

## 3) Order lipid classes globally by abundance for nicer stacked bars
lipid_order <- all_data %>%
  dplyr::count(lipid_class, sort = TRUE) %>%
  pull(lipid_class)

## 4) Convert to presence/absence matrix for UpSet
presence <- all_data %>%
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from = group,
    values_from = value,
    values_fill = FALSE
  ) %>%
  mutate(
    lipid_class = factor(lipid_class, levels = lipid_order)
  )

## 5) ComplexUpset plot
upset(
  presence,
  intersect = c("RDH12_Veh", "RDH12_100", "RDH12_200",
                "Ctrl_Veh", "Ctrl_100", "Ctrl_200"),
  name = "Group",
  width_ratio = 0.2,
  sort_intersections_by = "cardinality",
  min_size = 1,
  base_annotations = list(
    "Intersection size" =
      intersection_size(text = list(size = 3.8)) +
      ylab("Number of precursors") +
      theme(
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  ),
  annotations = list(
    "Lipid class composition" = (
      ggplot(mapping = aes(fill = lipid_class)) +
        geom_bar(stat = "count", position = "stack", width = 0.8) +
        ylab("Lipid class composition") +
        theme_bw(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right"
        ) +
        scale_fill_d3(palette = "category20")
    )
  )
)
#=======================================
Recovery_RDH12_Veh_filtered_summary2 <- Recovery_RDH12_veh_v4  %>%
  group_by( lipid_class, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_100_filtered_summary2 <- Recovery_RDH12_100_v4  %>%
  group_by( lipid_class, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_RDH12_200_filtered_summary2 <- Recovery_RDH12_200_v4  %>%
  group_by( lipid_class, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_veh_filtered_summary2 <- Recovery_Ctrl_veh_v4  %>%
  group_by( lipid_class, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_100_filtered_summary2 <- Recovery_Ctrl_100_v4  %>%
  group_by( lipid_class, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Recovery_Ctrl_200_filtered_summary2 <- Recovery_Ctrl_200_v4  %>%
  group_by( lipid_class, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))
#============================================
summarize_lipids <- function(df, label) {
  df %>%
    group_by(lipid_class, unsat_class) %>%
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
write.csv(combined_summary2, "Lipidomics/output_txts/Recovery_lipidclass_eachTreatment_UnSatLevel_summary.csv")

#++++++++++++++++++++++++++++++
#Create Plots

combined_FA <- combined_summary %>%
  group_by(group, lipid_class) %>%
  mutate(
    proportion = distinct_precursor / sum(distinct_precursor)
  ) %>%
  ungroup()


ggplot(combined_FA,
       aes(x = group, y = proportion, fill = unsat_class)) +
  geom_col() +
  facet_wrap(~lipid_class) +
  theme_bw(base_size = 12) +
  labs(
    y = "Proportion",
    fill = "FA class"
  )
ggplot(combined_FA,
       aes(x = proportion, y = group, fill = unsat_class)) +
  geom_col() +
  facet_wrap(~lipid_class) +
  theme_bw(base_size = 12) +
  labs(
    y = "Proportion",
    fill = "FA class"
  )

write.csv(combined_FA, "Lipidomics/output_txts/Recovery_lipidclass_eachTreatment_UnSatLevelProportion_summary.csv")


combined_total <- combined_FA %>%
  group_by(group, unsat_class) %>%
  summarise(total = sum(distinct_precursor), .groups = "drop") %>%
  group_by(group) %>%
  mutate(prop = total / sum(total))

ggplot(combined_total,
       aes(x = group, y = prop, fill = unsat_class)) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(y = "Overall FA composition")
