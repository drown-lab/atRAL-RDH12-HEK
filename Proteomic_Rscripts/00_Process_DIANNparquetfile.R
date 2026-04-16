
##read parquet
library(arrow)

HekatRAL_combndDIANN_dinj_v1<- read_parquet("/report.parquet")
HekatRAL_combndDIANN_dinj_v2 <-HekatRAL_combndDIANN_dinj_v1


HekatRAL_combndDIANN_dinj_v2 <-HekatRAL_combndDIANN_dinj_v2 |>
  filter(Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & Channel.Q.Value <= 0.05)|>
  filter(Proteotypic == 1)|>
  filter(!grepl("cRAP", Protein.Ids, ignore.case = TRUE))

#filter unique, and keep MaxLFQ values for quant
HekatRAL_combndDIANN_dinj_v3 <- HekatRAL_combndDIANN_dinj_v2 |>
  select(Run, Protein.Group, Channel, Genes, PG.MaxLFQ) |>
  unique() |>
  pivot_wider(values_from = PG.MaxLFQ, names_from = Run)

HekatRAL_combndDIANN_dinj_v3 <- HekatRAL_combndDIANN_dinj_v2 |>
  select(Run, Protein.Group, Genes, PG.MaxLFQ) |>
  unique()

HekatRAL_combndDIANN_dinj_v4 <- HekatRAL_combndDIANN_dinj_v3 %>%
  mutate(
    Treatment = case_when(
      str_detect(Run, "100u") ~ "100",
      str_detect(Run, "200u") ~ "200",
      str_detect(Run, "EtOH") ~ "control",
      str_detect(Run, "Control") ~ "control",
      
      TRUE ~ NA_character_
    ),
    Rep = case_when(
      str_detect(Run, "Rep1") ~ "a",
      str_detect(Run, "Rep2") ~ "b",
      str_detect(Run, "Rep3") ~ "c",
      str_detect(Run, "RepA") ~ "a",
      str_detect(Run, "RepB") ~ "b",
      str_detect(Run, "RepC") ~ "c",
      TRUE ~ NA_character_
    ),
    Genetype = case_when(
      str_detect(Run, "RDH12")~"RDH12",
      str_detect(Run, "GFP")~"GFP",
      str_detect(Run, "0703")~"RDH12",
    ),
    ExpType = 
      case_when(
        str_detect(Run, "RDH12")~"atRAL5hr+24h_recvr",
        str_detect(Run, "GFP")~"atRAL5hr+24h_recvr",
        str_detect(Run, "0703")~"atRAL5hr",
      )
  )

HekatRAL_combndDIANN_dinj_v4 <- HekatRAL_combndDIANN_dinj_v4 %>%
  mutate(
    Rep = case_when(
      str_detect(Run, "RepA") ~ "a",
      str_detect(Run, "RepB") ~ "b",
      str_detect(Run, "RepC") ~ "c",
      TRUE ~ NA_character_
    )
  )


write.csv(HekatRAL_combndDIANN_dinj_v4, "HekRDH12andGFP_combinedDatasets_atRAL5hr_with24hrRecvry_rawdata.csv")