library(dplyr)
library(tidyr)

# ----------------------------
# 1) Per protein × condition : 0/3, 1/3, 2/3, 3/3
# ----------------------------
protein_class_by_cond <- status_all %>%
  mutate(detected_int = as.integer(detected)) %>%   # TRUE/FALSE -> 1/0
  group_by( Protein.Group, condition) %>%
  summarise(
    n_reps   = dplyr::n(),                 # should be 3
    n_detect = sum(detected_int, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  # keep only conditions with exactly 3 reps (otherwise you get Present_2of2 etc.)
 # filter(n_reps == 3) %>%
  mutate(
    missing_class_3rep = case_when(
      n_detect == 0 ~ "MNAR_0of3",
      n_detect == 1 ~ "MAR_1of3",
      n_detect == 2 ~ "MAR_2of3",
      n_detect == 3 ~ "Present_3of3",
      TRUE ~ NA_character_
    ),
    missing_class_3rep = factor(
      missing_class_3rep,
      levels = c("MNAR_0of3", "MAR_1of3", "MAR_2of3", "Present_3of3")
    )
  )

# peek
head(protein_class_by_cond)

# ----------------------------
# 2) Counts per condition 
# ----------------------------
counts_by_condition <- protein_class_by_cond %>%
  dplyr::count( condition, missing_class_3rep, name = "n_proteins") %>%
  tidyr::pivot_wider(
    names_from  = missing_class_3rep,
    values_from = n_proteins,
    values_fill = 0
  ) %>%
  arrange( condition)

counts_by_condition

# ----------------------------
# 3) Overall totals 
# ----------------------------
counts_overall <- protein_class_by_cond %>%
  dplyr::count( missing_class_3rep, name = "n_protein_condition_blocks") %>%
  arrange( missing_class_3rep)

counts_overall

# ----------------------------
# 4) OPTIONAL: wide table, one row per protein, one column per condition
#    (useful to merge into DEP outputs by Protein.Group )
# ----------------------------
protein_class_wide <- protein_class_by_cond %>%
  select( Protein.Group, condition, missing_class_3rep) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = missing_class_3rep,
    names_prefix = "missclass_"
  )

protein_class_wide %>% head()

# ----------------------------
# 5) OPTIONAL: export
# ----------------------------
# write.csv(counts_by_condition, "Proteomic_output_txts/missingclass_counts_by_condition.csv", row.names = FALSE)
# write.csv(counts_overall,      "Proteomic_output_txts/missingclass_counts_overall.csv", row.names = FALSE)
# write.csv(protein_class_by_cond,"Proteomic_output_txts/protein_missingclass_by_condition.csv", row.names = FALSE)
# write.csv(protein_class_wide,  "Proteomic_output_txts/protein_missingclass_wide.csv", row.names = FALSE)
