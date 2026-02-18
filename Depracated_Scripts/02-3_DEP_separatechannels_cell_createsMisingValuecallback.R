library(dplyr)
library(stringr)

###After DEP analysis want to infer the missing values back

# det_by_cond has: channel, Protein.Group, condition, detect_class
#this table is made in Missingness for after DEP
miss_lookup <- det_by_cond %>%
  transmute(
    Protein.Group,
    channel,
    condition = paste0("X", condition),   # remove paste0 if already has X
    detect_class = as.character(detect_class)
  ) %>%
  distinct()

conds <- c("X24h_UT_mglps", "X6h_UT_mglps",
           "X24h_MG132", "X6h_MG132")

data_results_miss <- data_results

for (cnd in conds) {
  
  tmp <- miss_lookup %>%
    filter(condition == cnd) %>%
    select(Protein.Group, channel, detect_class)
  
  colname_new <- paste0(cnd, "_missing_class")
  
  names(tmp)[names(tmp) == "detect_class"] <- colname_new
  
  data_results_miss <- data_results_miss %>%
    left_join(
      tmp,
      by = c("ID" = "Protein.Group", "channel" = "channel")
    )
}



data_results_miss <- data_results_miss %>%
  mutate(across(
    .cols = ends_with("_missing_class"),
    .fns  = ~ ifelse(is.na(.x), "MNAR_0", .x)
  ))
