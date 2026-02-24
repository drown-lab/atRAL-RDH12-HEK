

###Join data_results and protein missing class
#protein_class_wide

library(dplyr)

# Make sure the join keys are the same type
protein_class_wide2 <- protein_class_wide %>%
  mutate(
    Protein.Group = as.character(Protein.Group),
    
  )

data_results2 <- data_results %>%
  mutate(
    Gene = as.character(name),
    
  )

# Left join: DEP results keep all rows, missingness columns added
data_results_w_missingclass <- data_results2 %>%
  left_join(
    protein_class_wide2,
    by = c("Gene" = "Protein.Group")
  )

# Check
dplyr::glimpse(data_results_w_missingclass)
head(data_results_w_missingclass)



# Optional: how many DEP rows did NOT find a missingness entry?
sum(is.na(data_results_w_missingclass$missclass_RDH12_control_atRAL5hr ))
sum(is.na(data_results_w_missingclass$missclass_RDH12_control_atRAL5hr.24h_recvr))
sum(is.na(data_results_w_missingclass$missclass_RDH12_200_atRAL5hr.24h_recvr))
sum(is.na(data_results_w_missingclass$missclass_RDH12_100_atRAL5hr))

#write.csv(data_results_w_missingclass, file = "Proteomic_output_txts/DEPresults_v1.csv")
write.csv(data_results_w_missingclass, path_DEP_results)
