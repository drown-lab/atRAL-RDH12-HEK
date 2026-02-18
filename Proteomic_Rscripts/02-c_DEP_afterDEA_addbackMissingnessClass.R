library(dplyr)

miss_lookup <- det_by_cond %>%
  transmute(
    protein,
    condition,
    detect_class = as.character(detect_class)
  ) %>%
  distinct()

conds <- c(
  "RDH12_100_atRAL5hr","RDH12_200_atRAL5hr","RDH12_control_atRAL5hr",
  "GFP_100_atRAL5hr.24h_recvr","GFP_200_atRAL5hr.24h_recvr","GFP_control_atRAL5hr.24h_recvr",
  "RDH12_100_atRAL5hr.24h_recvr","RDH12_200_atRAL5hr.24h_recvr","RDH12_control_atRAL5hr.24h_recvr"
)

data_results_miss <- data_results

for (cnd in conds) {
  
  tmp <- miss_lookup %>%
    filter(condition == cnd) %>%
    select(protein, detect_class)
  
  colname_new <- paste0(cnd, "_missing_class")
  names(tmp)[names(tmp) == "detect_class"] <- colname_new
  
  data_results_miss <- data_results_miss %>%
    left_join(tmp, by = c("name" = "protein"))
}

# sanity check: should now be filled
data_results_miss %>%
  select(name, ID, ends_with("_missing_class")) %>%
  head()



saveRDS(data_results_miss, file = "RData/DEP_results_wMissingnessclass.rds")

write.csv(data_results_miss,"Proteomic_output_txts/DEP_results_wMissingnessclass.csv")
