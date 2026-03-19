
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)
library(ggsci)

##Breakdown Acute data into Satand UnsatLevel by Treatment

df<-read.csv("Lipidomics/output_txts/Acute_filteredlipidtable.csv")

df <- df %>%
  mutate(
    unsat_class = case_when(
      DB_total == 0 ~ "SFA",
      DB_total == 1 ~ "MUFA",
      DB_total >= 2 ~ "PUFA",
      TRUE ~ NA_character_
    )
  )
Acute100_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s4, s5, s6)))

Acute100_v2<-Acute100_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
  #Standard cut-off 
  Acute100_v3<-Acute100_v2|>
  filter(max_divid_blank >= 1.3)

Acute100_v4<-Acute100_v3|>
  unique()
  
  ##Plot Distribution of Lipids
  ggplot(Acute100_v2, aes(x= max_divid_blank, fill = lipid_class1))+
  geom_histogram(bins = 40)+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Acute: Distribution of MRMs with Signal-to-Blank",
    subtitle = "Dashed line: 1.3x "
  )+
  scale_y_continuous(trans = "log10")
  
  Acute100_filtered_summary <- Acute100_v4  %>%
    group_by( lipid_class) %>%
    summarise(distinct_mrms = n_distinct(mrm1))
  
  Acute100_filtered_summary2 <- Acute100_v4  %>%
    group_by( lipid_class) %>%
    summarise(distinct_precursor = n_distinct(precursor))
  colnames(Acute_filtered_summary2)
  
  Acute100_filtered_summary2 <- Acute100_filtered_summary2 %>%
    arrange(desc(distinct_precursor)) %>%
    mutate(lipid_class = factor(lipid_class, levels = lipid_class))
  
  Acute200_v1<-df|>
    rowwise() |>
    mutate(maxvalue = max(c(s7, s8, s9)))
  
  Acute200_v2<-Acute200_v1|>
    mutate(
      max_divid_blank = maxvalue/blank
    )
  #Filter by 30% higher than the blank
  #Standard cut-off 
  Acute200_v3<-Acute200_v2|>
    filter(max_divid_blank >= 1.3)
  
  Acute200_v4<-Acute200_v3|>
    unique()
  
  Acute200_filtered_summary <- Acute200_v4  %>%
    group_by( lipid_class) %>%
    summarise(distinct_mrms = n_distinct(mrm1))
  
  Acute200_filtered_summary2 <- Acute200_v4  %>%
    group_by( lipid_class) %>%
    summarise(distinct_precursor = n_distinct(precursor))

  Acute200_filtered_summary2 <- Acute200_filtered_summary2 %>%
    arrange(desc(distinct_precursor)) %>%
    mutate(lipid_class = factor(lipid_class, levels = lipid_class))
  ##########+===================
  ##########+

  AcuteVeh_v1<-df|>
    rowwise() |>
    mutate(maxvalue = max(c(s1, s2, s3)))
  
  AcuteVeh_v2<-AcuteVeh_v1|>
    mutate(
      max_divid_blank = maxvalue/blank
    )
  #Filter by 30% higher than the blank
  #Standard cut-off 
  AcuteVeh_v3<-AcuteVeh_v2|>
    filter(max_divid_blank >= 1.3)
  
  AcuteVeh_v4<-AcuteVeh_v3|>
    unique()
  
  AcuteVeh_filtered_summary <- AcuteVeh_v4  %>%
    group_by( lipid_class) %>%
    summarise(distinct_mrms = n_distinct(mrm1))
  
  AcuteVeh_filtered_summary2 <- AcuteVeh_v4  %>%
    group_by( lipid_class) %>%
    summarise(distinct_precursor = n_distinct(precursor))
  
  AcuteVeh_filtered_summary2 <- AcuteVeh_filtered_summary2 %>%
    arrange(desc(distinct_precursor)) %>%
    mutate(lipid_class = factor(lipid_class, levels = lipid_class))  
  
  write.csv(AcuteVeh_filtered_summary2, "Lipidomics/output_txts/AcuteVeh_LipidClass_summary.csv")
  write.csv(Acute100_filtered_summary2,"Lipidomics/output_txts/Acute100_LipidClass_summary.csv")
 write.csv(Acute200_filtered_summary2, "Lipidomics/output_txts/Acute200_LipidClass_summary.csv") 
  
  Acute100_v4
  Acute200_v4
  colnames(AcuteVeh_v4)
  

 
##Upset Plot Across Lipid IDS for treatments  
  all_data <- bind_rows(
    AcuteVeh_v4  %>% mutate(group = "Veh"),
    Acute100_v4 %>% mutate(group = "100"),
    Acute200_v4 %>% mutate(group = "200")
  )
  
  presence <- all_data %>%
    distinct(precursor, group) %>%
    mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = group, values_from = value, values_fill = 0)
  

  upset(
    presence,
    intersect = c("Veh", "100", "200")
  )
  

  all_data <- bind_rows(
    AcuteVeh_v4  %>% mutate(group = "Veh"),
    Acute100_v4  %>% mutate(group = "100"),
    Acute200_v4  %>% mutate(group = "200")
  ) %>%
    select(precursor, lipid_class, group) %>%
    distinct()
  
  presence <- all_data %>%
    mutate(value = TRUE) %>%
    pivot_wider(
      names_from = group,
      values_from = value,
      values_fill = FALSE
    )
  
  upset(
    presence,
    intersect = c("Veh", "100", "200"),
    base_annotations = list(
      "Intersection size" = intersection_size()
    ),
    annotations = list(
      "Lipid class" = (
        ggplot(mapping = aes(fill = lipid_class)) +
          geom_bar(stat = "count", position = "stack")
      )
    )
  )
  
  ##Upset Plot Across Lipid IDS for treatments  
  all_data <- bind_rows(
    AcuteVeh_v4  %>% mutate(group = "Veh"),
    Acute100_v4 %>% mutate(group = "100"),
    Acute200_v4 %>% mutate(group = "200")
  )
  
  presence <- all_data %>%
    distinct(precursor, group) %>%
    mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = group, values_from = value, values_fill = 0)
  

  upset(
    presence,
    intersect = c("Veh", "100", "200")
  )
  

  all_data <- bind_rows(
    AcuteVeh_v4  %>% mutate(group = "Veh"),
    Acute100_v4  %>% mutate(group = "100"),
    Acute200_v4  %>% mutate(group = "200")
  ) %>%
    select(precursor, lipid_class, group) %>%
    distinct()
  
  presence <- all_data %>%
    mutate(value = TRUE) %>%
    pivot_wider(
      names_from = group,
      values_from = value,
      values_fill = FALSE
    )
  
  upset(
    presence,
    intersect = c("Veh", "100", "200"),
    base_annotations = list(
      "Intersection size" = intersection_size()
    ),
    annotations = list(
      "Lipid class" = (
        ggplot(mapping = aes(fill = lipid_class)) +
          geom_bar(stat = "count", position = "stack")
      )
    )
  )
  
  
  
  all_data %>%
    distinct(precursor, lipid_class) %>%
    dplyr::count(precursor) %>%
    filter(n > 1)  
  

  ## UpSet plot across lipid IDs for treatments, annotated by lipid class
  all_data <- bind_rows(
    AcuteVeh_v4  %>% mutate(group = "Veh"),
    Acute100_v4  %>% mutate(group = "100"),
    Acute200_v4  %>% mutate(group = "200")
  ) %>%
    select(precursor, lipid_class, group) %>%
    distinct()
  
  presence <- all_data %>%
    mutate(value = TRUE) %>%
    pivot_wider(
      names_from = group,
      values_from = value,
      values_fill = FALSE
    ) %>%
    mutate(
      lipid_class = factor(lipid_class)
    )
  
  upset(
    presence,
    intersect = c("Veh", "100", "200"),
    name = "Treatment",
    width_ratio = 0.18,
    sort_intersections_by = "cardinality",
    min_size = 1,
    base_annotations = list(
      "Intersection size" = intersection_size(
        text = list(size = 4)
      ) +
        ylab("Number of precursors") +
        theme(
          axis.title.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    ),
    annotations = list(
      "Lipid class" = (
        ggplot(mapping = aes(fill = lipid_class)) +
          geom_bar(stat = "count", position = "stack", width = 0.8) +
          ylab("Lipid class composition") +
          theme_bw(base_size = 12) +
          theme(
            panel.grid.minor = element_blank(),
            legend.position = "right",
            axis.title.x = element_blank()
          ) +
          scale_fill_d3(palette = "category20")
      )
    )
  )
#############################
  #Save FAclass type
  
  Acute200_filtered_summary_FAclass <- Acute200_v4  %>%
    group_by( lipid_class, unsat_class) %>%
    summarise(distinct_precursor = n_distinct(precursor))
  
  
  Acute100_filtered_summary_FAclass <- Acute100_v4  %>%
    group_by( lipid_class, unsat_class) %>%
    summarise(distinct_precursor = n_distinct(precursor))
  
  AcuteVeh_filtered_summary_FAclass <- AcuteVeh_v4  %>%
    group_by( lipid_class, unsat_class) %>%
    summarise(distinct_precursor = n_distinct(precursor))
  
  library(dplyr)
  
  AcuteVeh_filtered_summary_FAclass  <- AcuteVeh_filtered_summary_FAclass  %>% mutate(Treatment = "Veh")
  Acute100_filtered_summary_FAclass <- Acute100_filtered_summary_FAclass %>% mutate(Treatment = "100")
  Acute200_filtered_summary_FAclass <- Acute200_filtered_summary_FAclass %>% mutate(Treatment = "200")
  
  combined_FA <- bind_rows(
    AcuteVeh_filtered_summary_FAclass,
    Acute100_filtered_summary_FAclass,
    Acute200_filtered_summary_FAclass
  )
  
  combined_FA <- combined_FA %>%
    group_by(Treatment, lipid_class) %>%
    mutate(
      proportion = distinct_precursor / sum(distinct_precursor)
    ) %>%
    ungroup()
  

  ggplot(combined_FA,
         aes(x = Treatment, y = proportion, fill = unsat_class)) +
    geom_col() +
    facet_wrap(~lipid_class) +
    theme_bw(base_size = 12) +
    labs(
      y = "Proportion",
      fill = "FA class"
    )
  
  combined_total <- combined_FA %>%
    group_by(Treatment, unsat_class) %>%
    summarise(total = sum(distinct_precursor), .groups = "drop") %>%
    group_by(Treatment) %>%
    mutate(prop = total / sum(total))
  
  ggplot(combined_total,
         aes(x = Treatment, y = prop, fill = unsat_class)) +
    geom_col() +
    theme_bw(base_size = 12) +
    labs(y = "Overall FA composition")

  library(dplyr)
  library(ggplot2)
  
  combined_FA <- bind_rows(
    AcuteVeh_filtered_summary_FAclass  %>% mutate(Treatment = "Veh"),
    Acute100_filtered_summary_FAclass  %>% mutate(Treatment = "100"),
    Acute200_filtered_summary_FAclass  %>% mutate(Treatment = "200")
  ) %>%
    group_by(Treatment, lipid_class) %>%
    mutate(prop = distinct_precursor / sum(distinct_precursor)) %>%
    ungroup() %>%
    mutate(
      unsat_class = factor(unsat_class, levels = c("PUFA", "MUFA", "SFA")),
      Treatment = factor(Treatment, levels = c("Veh", "100", "200"))
    )
  
  ggplot(combined_FA,
         aes(x = prop, y = Treatment, fill = unsat_class)) +
    geom_col() +
    facet_wrap(~lipid_class , nrow = 2) +
    theme_bw(base_size = 12) +
    labs(
      x = "Proportion",
      y = "Treatment",
      fill = "Unsaturation class"
    )

 write.csv(combined_FA, "Lipidomics/output_txts/Acute_byTreatment_FAClass_summarize.csv") 

  
  combined_FA <- bind_rows(
    AcuteVeh_filtered_summary_FAclass  %>% mutate(Treatment = "Veh"),
    Acute100_filtered_summary_FAclass  %>% mutate(Treatment = "100"),
    Acute200_filtered_summary_FAclass  %>% mutate(Treatment = "200")
  ) %>%
    group_by(Treatment, lipid_class) %>%
    mutate(prop = distinct_precursor / sum(distinct_precursor)) %>%
    ungroup() %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Veh", "100", "200")),
      unsat_class = factor(unsat_class, levels = c("SFA", "MUFA", "PUFA"))
    ) %>%
    arrange(lipid_class, Treatment) %>%
    mutate(
      y_axis = paste(lipid_class, Treatment, sep = " - ")
    )
  
  ggplot(combined_FA,
         aes(x = prop, y = y_axis, fill = unsat_class)) +
    geom_col() +
    theme_bw(base_size = 12) +
    labs(
      x = "Proportion within lipid class",
      y = NULL,
      fill = "Unsaturation class"
    )
  