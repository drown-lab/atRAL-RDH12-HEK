Hek293_PPF_v1

Hek293_PPF_v2<-Hek293_PPF_v1|>
  pivot_longer(
    cols = c(GFP_1:RDH12_3),
    names_to = "SampleName",
    values_to = "PG.MaxLFQ"
  )

Hek293_PPF_v3 <- Hek293_PPF_v2 %>%
  mutate(
    CellType = case_when(
      str_detect(SampleName, "RDH12") ~ "RDH12",
      str_detect(SampleName, "GFP") ~ "GFP",
   
      TRUE ~ NA_character_
    ),
    Fraction = case_when(
      str_detect(SampleName, "_1") ~ "1",
      str_detect(SampleName, "_2") ~ "2",
      str_detect(SampleName, "_3") ~ "3",
      str_detect(SampleName, "-2") ~ "2",
      
      TRUE ~ NA_character_
    ))


Organoid_v5_fixed <- Orgnoid_v5 |>
  mutate(
    Protein.Ids = NA_character_,
    Protein.Names = NA_character_,
    First.Protein.Description = NA_character_
  ) |>
  select(colnames(Hek293_PPF_v3))

combined_proteome <- bind_rows(
  Hek293_PPF_v3,
  Organoid_v5_fixed
)

colnames(combined_proteome)
table(combined_proteome$CellType)
table(combined_proteome$SampleName)


