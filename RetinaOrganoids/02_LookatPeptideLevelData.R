
####
#organoid Peptide level data

Orgnoid_peptide_v1<-Orgnoid_v3
Orgnoid_peptide_v2 <- Orgnoid_peptide_v1 %>%
  mutate(
    CellType = case_when(
      str_detect(Run, "^1_Slot1-19") ~ "EP",
      str_detect(Run, "^2_Slot1-20") ~ "EP",
      str_detect(Run, "^3_Slot1-21") ~ "EP",
      str_detect(Run, "^4_Slot1-22") ~ "EP",
      str_detect(Run, "^5_Slot1-23") ~ "WO",
      str_detect(Run, "^6_Slot1-24") ~ "WO",
      str_detect(Run, "^7_Slot1-25") ~ "WO",
      str_detect(Run, "^8_Slot1-26") ~ "WO",
      TRUE ~ NA_character_
    ),
    Fraction = case_when(
      str_detect(Run, "^1_") ~ "1",
      str_detect(Run, "^2_") ~ "2",
      str_detect(Run, "^3_") ~ "3",
      str_detect(Run, "^4_") ~ "4",
      str_detect(Run, "^5_") ~ "1",
      str_detect(Run, "^6_") ~ "2",
      str_detect(Run, "^7_") ~ "3",
      str_detect(Run, "^8_") ~ "4",
      TRUE ~ NA_character_
    ),
    SampleName = case_when(
      str_detect(Run, "^1_") ~ "EP_F1",
      str_detect(Run, "^2_") ~ "EP_F2",
      str_detect(Run, "^3_") ~ "EP_F3",
      str_detect(Run, "^4_") ~ "EP_F4",
      str_detect(Run, "^5_") ~ "WO_F1",
      str_detect(Run, "^6_") ~ "WO_F2",
      str_detect(Run, "^7_") ~ "WO_F3",
      str_detect(Run, "^8_") ~ "WO_F4",
      TRUE ~ NA_character_
    )
  )

count_pg <- Orgnoid_peptide_v2|>
  filter(!is.na(Precursor.Normalised))|>
  filter(!Precursor.Normalised==0)|>
  group_by(SampleName)|>
  summarise(count = n_distinct(Precursor.Id))

Orgnoid_peptide_v3 <- Orgnoid_peptide_v2|>
  filter(!is.na(Precursor.Normalised))|>
  filter(!Precursor.Normalised==0)

Orgnoid_peptide_v3<-Orgnoid_peptide_v3|>
  select(Precursor.Id, Precursor.Charge, Stripped.Sequence, Precursor.Mz, Precursor.Quantity, Precursor.Normalised, RT, IM, CellType, Fraction, SampleName, Run)|>
  unique()


list_upset <- list(
  EP_F1 = filter(Orgnoid_peptide_v3, SampleName == "EP_F1" & !is.na(Precursor.Normalised))$Precursor.Id,
  EP_F2 = filter(Orgnoid_peptide_v3, SampleName == "EP_F2" & !is.na(Precursor.Normalised))$Precursor.Id,    
  EP_F3 = filter(Orgnoid_peptide_v3, SampleName == "EP_F3" & !is.na(Precursor.Normalised))$Precursor.Id,    
  EP_F4 = filter(Orgnoid_peptide_v3, SampleName == "EP_F4" & !is.na(Precursor.Normalised))$Precursor.Id,    
  WO_F1 = filter(Orgnoid_peptide_v3, SampleName == "WO_F1" & !is.na(Precursor.Normalised))$Precursor.Id,    
  WO_F2 = filter(Orgnoid_peptide_v3, SampleName == "WO_F2" & !is.na(Precursor.Normalised))$Precursor.Id,    
  WO_F3 = filter(Orgnoid_peptide_v3, SampleName == "WO_F3" & !is.na(Precursor.Normalised))$Precursor.Id,    
  WO_F4 = filter(Orgnoid_peptide_v3, SampleName == "WO_F4" & !is.na(Precursor.Normalised))$Precursor.Id    
)

# Run UpSetR
# 2. Convert to input for UpSetR
upset_data <- fromList(list_upset)

# 3. Plot
upset(upset_data, order.by = "freq", nsets = 8, nintersects =20)
#upset(upset_data, order.by = "freq", nsets = )


# build sets
ep <- Orgnoid_peptide_v3 |>
  filter(str_detect(SampleName, "^EP_"), !is.na(Precursor.Normalised)) |>
  pull(Precursor.Id) |>
  unique()

wo <- Orgnoid_peptide_v3 |>
  filter(str_detect(SampleName, "^WO_"), !is.na(Precursor.Normalised)) |>
  pull(Precursor.Id) |>
  unique()

venn_plot <- venn.diagram(
  x = list(EP = ep, WO = wo),
  filename = NULL,
  fill = c("lightskyblue", "#FF34B3"),
  alpha = 0.5,
  cex = 1.3,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  margin = 0.08
)

grid.newpage()
grid.draw(venn_plot)

#=====================
#Creation of Venn Diagrams


list_ep <- list(
  EP_F1 = Orgnoid_peptide_v3 |> filter(SampleName == "EP_F1", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique(),
  EP_F2 = Orgnoid_peptide_v3 |> filter(SampleName == "EP_F2", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique(),
  EP_F3 = Orgnoid_peptide_v3 |> filter(SampleName == "EP_F3", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique(),
  EP_F4 = Orgnoid_peptide_v3 |> filter(SampleName == "EP_F4", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique()
)

venn_plot <- venn.diagram(
  x = list_ep,
  filename = NULL,
  fill = c("#FF6DB6", "#009292", "#920000", "#24FF24"),
  alpha = 0.45,
  cex = 0.9,
  cat.cex = 0.9,
  margin = 0.08
)

grid.newpage()
grid.draw(venn_plot)


list_ep <- list(
  WO_F1 = Orgnoid_peptide_v3 |> filter(SampleName == "WO_F1", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique(),
  WO_F2 = Orgnoid_peptide_v3 |> filter(SampleName == "WO_F2", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique(),
  WO_F3 = Orgnoid_peptide_v3 |> filter(SampleName == "WO_F3", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique(),
  WO_F4 = Orgnoid_peptide_v3 |> filter(SampleName == "WO_F4", !is.na(Precursor.Normalised)) |> pull(Precursor.Id) |> unique()
  
)

venn_plot <- venn.diagram(
  x = list_ep,
  filename = NULL,
  fill = c("#FF6DB6", "#009292", "#920000", "#24FF24"),
  alpha = 0.45,
  cex = 0.9,
  cat.cex = 0.9,
  margin = 0.08
)

grid.newpage()
grid.draw(venn_plot)
#=======================================================================
# this gets the proteins unique to EP and to WO
# EP-only proteins
ep_only <- setdiff(ep, wo)

# WO-only proteins
wo_only <- setdiff(wo, ep)

# Combine into a tidy table
unique_tbl_pep <- tibble::tibble(
  Precursor.Id = c(ep_only, wo_only),
  Unique_to = c(
    rep("EP", length(ep_only)),
    rep("WO", length(wo_only))
  )
)

unique_tbl_pep



# 2) Add back abundance rows where observed (many-to-many)
unique_long_pep <- unique_tbl_pep %>%
  left_join(
    Orgnoid_peptide_v3 %>%
     
      mutate(Group = case_when(
        str_detect(SampleName, "^EP_") ~ "EP",
        str_detect(SampleName, "^WO_") ~ "WO",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(Group)) %>%
      select(Precursor.Id, SampleName, Group, CellType, Fraction, Precursor.Normalised),
    by = "Precursor.Id",
    relationship = "many-to-many"
  ) %>%
  # keep only the rows that match the uniqueness label (optional but usually desired)
  filter(Group == Unique_to)

unique_long_pep



#$+++++++++++++++++++++

###Plots of Overal Abundance across all samples
ggplot(unique_long_pep, aes(x= Fraction, y=log10(Precursor.Normalised), fill= Fraction))+
  geom_boxplot(position = position_dodge())+
  facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)

ggplot(unique_long_pep, aes(x= Fraction, y=log10(Precursor.Normalised), fill= Fraction))+
  geom_violin(position = position_dodge())+
  facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)

ggplot(unique_long_pep, aes(x= CellType, y=log10(Precursor.Normalised), fill= CellType))+
  geom_boxplot(position = position_dodge())+
  # facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)


ggplot(unique_long_pep, aes(x= CellType, y=log10(Precursor.Normalised), fill= CellType))+
  geom_violin(position = position_dodge())+
  # facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)
