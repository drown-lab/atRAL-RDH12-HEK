

colnames(df_with_unsat_updated)
[1] "lipid_name"        "mrm"               "blank"             "s1"                "s2"                "s3"               
[7] "s4"                "s5"                "s6"                "s7"                "s8"                "s9"               
[13] "lipid_class1"      "lipid_class"       "DG_TG_acyl_chains" "NL_chain"          "maxvalue"          "max_divid_blank"  
[19] "mrm1"              "precursor"         "product"           "picked_candidate"  "C_total"           "DB_total"         
[25] "NL_C"              "NL_DB"             "cand_score"

df<-read.csv("Lipidomics/output_txts/Acute_filteredlipidtable.csv")

colnames(df)

Acute_filtered_summary <- df  %>%
  group_by( lipid_class) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

Acute_filtered_summary2 <- df  %>%
  group_by( lipid_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))
colnames(Acute_filtered_summary2)

Acute_filtered_summary2 <- Acute_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class = factor(lipid_class, levels = lipid_class))


Acute_filtered_summary2 <- Acute_filtered_summary2 %>%
  mutate(
    percent = distinct_precursor / sum(distinct_precursor),
    label = paste0(lipid_class, " (", scales::percent(percent, accuracy = 1), ")")
  )

ggplot(Acute_filtered_summary2, 
       aes(x = "", y = distinct_precursor, fill = lipid_class)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = scales::percent(percent)),
            position = position_stack(vjust = 0.5),
            size = 3) +
  labs(fill = "Lipid Class")
##########################################
#Obtain SFA, MUFA, PUFA

df <- df %>%
  mutate(
    unsat_class = case_when(
      DB_total == 0 ~ "SFA",
      DB_total == 1 ~ "MUFA",
      DB_total >= 2 ~ "PUFA",
      TRUE ~ NA_character_
    )
  )

Acute_filtered_summary_FAclass <- df  %>%
  group_by( lipid_class, unsat_class) %>%
  summarise(distinct_precursor = n_distinct(precursor))
colnames(Acute_filtered_summary_FAclass)


Acute_filtered_summary_FAclass <- Acute_filtered_summary_FAclass %>%
  mutate(
    unsat_class = factor(unsat_class, levels = c("SFA", "MUFA", "PUFA"))
  )
write.csv(Acute_filtered_summary_FAclass, "Lipidomics/output_txts/Acute_LipidUnsatClass_summary.csv")

ggplot(Acute_filtered_summary_FAclass,
       aes(x = lipid_class, y = distinct_precursor, fill = unsat_class)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Lipid class",
    y = "Number of distinct precursors",
    fill = "FA class"
  )

ggplot(Acute_filtered_summary_FAclass,
       aes(x = lipid_class, y = distinct_precursor, fill = unsat_class)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12)


ggplot(df, aes(x=C_total, y= DB_total, color =lipid_class1))+
  geom_point(alpha=0.8, size = 2)+
  theme_bw(base_size = 12)+
  facet_wrap(~lipid_class1)

ggplot(df, aes(x=C_total, y= DB_total, color =max_divid_blank))+
  geom_point(alpha=0.8, size = 2)+
  theme_bw(base_size = 12)+
  facet_wrap(~lipid_class1)

ggplot(df, aes(x = C_total, y = DB_total, color = maxvalue)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_viridis_c(trans = "log10", option = "plasma") +
  theme_bw(base_size = 12) +
  facet_wrap(~lipid_class1)

ggplot(df, aes(x = C_total, y = DB_total, color = max_divid_blank)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_viridis_c(trans = "log10", option = "plasma") +
  theme_bw(base_size = 12) +
  facet_wrap(~lipid_class1)

plot_df<-df|>
  filter(lipid_class1 %in% c("PC", "SM", "PE", "PS","PI", "TG", "DG"))


ggplot(plot_df, aes(x = C_total, y = DB_total, color = maxvalue)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_viridis_c(trans = "log10", option = "plasma") +
  theme_bw(base_size = 12) +
  facet_wrap(~lipid_class1, nrow=2)+
  theme(
    panel.grid.minor = element_blank())+
   scale_y_continuous(breaks = seq(0, max(plot_df$DB_total, na.rm = TRUE), by = 2))


plot_count <- plot_df %>%
  group_by(lipid_class1, C_total, DB_total) %>%
  summarise(
    n_lipids = n(),
    maxvalue = max(maxvalue, na.rm = TRUE),  # or mean() if preferred
    .groups = "drop"
  )

ggplot(plot_count, aes(x = C_total, y = DB_total)) +
  geom_point(aes(size = n_lipids, color = maxvalue), alpha = 0.7) +
  scale_color_viridis_c(trans = "log10", option = "plasma") +
  scale_size_continuous(range = c(2, 8)) +  # adjust point sizes
  theme_bw(base_size = 12) +
  facet_wrap(~lipid_class1, nrow = 2) +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(
    breaks = seq(0, max(plot_count$DB_total, na.rm = TRUE), by = 2)
  ) +
  labs(
    size = "Number of lipids",
    color = "Max intensity"
  )