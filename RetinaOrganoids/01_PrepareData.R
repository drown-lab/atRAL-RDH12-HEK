##read parquet
##load libraries
library(arrow)
library(reshape2)
library(readxl)
library(tidyverse)
library(readr)
library(dplyr)
#######################

Orgnoid_v1 <- read_parquet("C:/Users/LabUser/Desktop/StemCells/RamsCollab_RetinaOrganoids/DIANN2.3.2_output/report.parquet")
exists("Orgnoid_v1")
colnames(Orgnoid_v1)

#filter PG based on certain values
##filter for only proteotypic and get rid of "crap" or decoy protein IDs

Orgnoid_v2 <- Orgnoid_v1 |>
  filter(Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & Channel.Q.Value <= 0.05)|>
  filter(Proteotypic == 1)|>
  filter(!grepl("cRAP", Protein.Ids, ignore.case = TRUE))

#known human contaminants:
contaminant_proteins <- c(
  "Q6E0U4", "P20930", "Q5D862", "Q86YZ3", "P04264", "P13645", "Q99456",
  "P13646", "P02533", "P19012", "P08779", "Q04695", "P05783", "P08727",
  "P35908", "P35900", "Q8N1A0", "Q9C075", "Q2M2I5", "Q7Z3Z0", "Q7Z3Y9",
  "Q7Z3Y8", "Q7Z3Y7", "P12035", "Q15323", "Q14532", "O76009", "Q14525",
  "O76011", "Q92764", "O76013", "O76014", "O76015", "Q6A163", "P19013",
  "Q6A162", "P13647", "P02538", "P04259", "P48668", "Q3KNV1", "P08729",
  "Q3SY84", "Q14CN4", "Q86Y46", "Q7RTS7", "O95678", "Q01546", "Q7Z794",
  "Q8N1N4", "Q5XKE5", "P05787", "Q6KB66", "Q14533", "Q9NSB4", "P78385",
  "Q9NSB2", "P78386", "O43790", "A6NCN2", "P35527", "P60331", "P60014",
  "P60412", "P60413", "P60368", "P60369", "P60372", "P60370", "P60371",
  "P60409", "P60410", "P60411", "Q07627", "Q8IUC1", "P59990", "P59991",
  "P60328", "P60329", "Q8IUG1", "Q8IUC0", "Q52LG2", "Q3SY46", "Q3LI77",
  "P0C5Y4", "Q9BYS1", "Q3LI76", "A8MUX0", "Q9BYP8", "Q8IUB9", "Q3LHN2",
  "Q7Z4W3", "Q3LI73", "Q3LI72", "Q3LI70", "Q3SYF9", "Q3LI54", "Q3LI63",
  "Q3LI61", "Q9BYU5", "Q3LI58", "Q3LI59", "Q3LHN1", "Q9BYT5", "Q3MIV0",
  "P0C7H8", "A1A580", "Q9BYR9", "Q3LI83", "Q3LHN0", "Q6PEX3", "Q3LI81",
  "A8MX34", "Q9BYR8", "Q9BYR7", "Q9BYR6", "Q9BYQ7", "Q9BYQ6", "Q9BQ66",
  "Q9BYR5", "Q9BYR4", "Q9BYR3", "Q9BYR2", "Q9BYQ5", "Q9BYR0", "Q9BYQ9",
  "Q9BYQ8", "Q6L8H4", "Q6L8G5", "Q6L8G4", "Q701N4", "Q6L8H2", "Q6L8H1",
  "Q701N2", "Q6L8G9", "Q6L8G8", "O75690", "P26371", "Q3LI64", "Q3LI66",
  "Q3LI67", "Q8IUC3", "Q8IUC2", "A8MXZ3", "Q9BYQ4", "Q9BYQ3", "Q9BYQ2",
  "A8MVA2", "A8MTY7", "Q9BYQ0", "Q9BYP9"
)

Orgnoid_v3 <- Orgnoid_v2 |>
  filter(!Protein.Group %in% contaminant_proteins)

Orgnoid_v4 <- Orgnoid_v3 |>
  select(Run, Protein.Group, Genes, PG.MaxLFQ) |>
  unique()

Orgnoid_v5 <- Orgnoid_v4 %>%
  select(-Run)
  

Orgnoid_v4 <- Orgnoid_v4 %>%
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

##obtain count of values in the table
count_pg <- Orgnoid_v4|>
  filter(!is.na(PG.MaxLFQ))|>
  filter(!PG.MaxLFQ==0)|>
  group_by(SampleName)|>
  summarise(count = n_distinct(Protein.Group))

Orgnoid_v5<-Orgnoid_v4|>
  filter(!is.na(PG.MaxLFQ))|>
  filter(!PG.MaxLFQ==0)


##Make list for upset plot

list_upset <- list(
  EP_F1 = filter(Orgnoid_v5, SampleName == "EP_F1" & !is.na(PG.MaxLFQ))$Protein.Group,
  EP_F2 = filter(Orgnoid_v5, SampleName == "EP_F2" & !is.na(PG.MaxLFQ))$Protein.Group,    
  EP_F3 = filter(Orgnoid_v5, SampleName == "EP_F3" & !is.na(PG.MaxLFQ))$Protein.Group,    
  EP_F4 = filter(Orgnoid_v5, SampleName == "EP_F4" & !is.na(PG.MaxLFQ))$Protein.Group,    
  WO_F1 = filter(Orgnoid_v5, SampleName == "WO_F1" & !is.na(PG.MaxLFQ))$Protein.Group,    
  WO_F2 = filter(Orgnoid_v5, SampleName == "WO_F2" & !is.na(PG.MaxLFQ))$Protein.Group,    
  WO_F3 = filter(Orgnoid_v5, SampleName == "WO_F3" & !is.na(PG.MaxLFQ))$Protein.Group,    
  WO_F4 = filter(Orgnoid_v5, SampleName == "WO_F4" & !is.na(PG.MaxLFQ))$Protein.Group    
)

# Run UpSetR
# 2. Convert to input for UpSetR
upset_data <- fromList(list_upset)

# 3. Plot
upset(upset_data, order.by = "freq", nsets = 8, nintersects =20)
upset(upset_data, order.by = "freq", nsets = )

library(dplyr)
library(stringr)
library(VennDiagram)
library(grid)

# build sets
ep <- Orgnoid_v5 |>
  filter(str_detect(SampleName, "^EP_"), !is.na(PG.MaxLFQ)) |>
  pull(Protein.Group) |>
  unique()

wo <- Orgnoid_v5 |>
  filter(str_detect(SampleName, "^WO_"), !is.na(PG.MaxLFQ)) |>
  pull(Protein.Group) |>
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
  EP_F1 = Orgnoid_v5 |> filter(SampleName == "EP_F1", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  EP_F2 = Orgnoid_v5 |> filter(SampleName == "EP_F2", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  EP_F3 = Orgnoid_v5 |> filter(SampleName == "EP_F3", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  EP_F4 = Orgnoid_v5 |> filter(SampleName == "EP_F4", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique()
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
  WO_F1 = Orgnoid_v5 |> filter(SampleName == "WO_F1", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  WO_F2 = Orgnoid_v5 |> filter(SampleName == "WO_F2", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  WO_F3 = Orgnoid_v5 |> filter(SampleName == "WO_F3", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  WO_F4 = Orgnoid_v5 |> filter(SampleName == "WO_F4", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique()

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

#++++++++
#+Doesn't work
list_ep <- list(
  EP_F1 = Orgnoid_v4 |> filter(SampleName == "EP_F1", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  EP_F2 = Orgnoid_v4 |> filter(SampleName == "EP_F2", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  EP_F3 = Orgnoid_v4 |> filter(SampleName == "EP_F3", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  EP_F4 = Orgnoid_v4 |> filter(SampleName == "EP_F4", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  WO_F1 = Orgnoid_v4 |> filter(SampleName == "WO_F1", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  WO_F2 = Orgnoid_v4 |> filter(SampleName == "WO_F2", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  WO_F3 = Orgnoid_v4 |> filter(SampleName == "WO_F3", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique(),
  WO_F4 = Orgnoid_v4 |> filter(SampleName == "WO_F4", !is.na(PG.MaxLFQ)) |> pull(Protein.Group) |> unique()
  
)

venn_plot <- venn.diagram(
  x = list_ep,
  filename = NULL,
  fill = c("#FF6DB6", "#009292", "#920000", "#24FF24","#FFFF6D","#B6DBFF","#FFB6DB","#006DDB"),
  alpha = 0.45,
  cex = 0.9,
  cat.cex = 0.9,
  margin = 0.08
)

grid.newpage()
grid.draw(venn_plot)

#===============================
#Creates upset plots for each EP or WO
list_upset <- list(
  EP_F1 = filter(Orgnoid_v5, SampleName == "EP_F1" & !is.na(PG.MaxLFQ))$Protein.Group,
  EP_F2 = filter(Orgnoid_v5, SampleName == "EP_F2" & !is.na(PG.MaxLFQ))$Protein.Group,    
  EP_F3 = filter(Orgnoid_v5, SampleName == "EP_F3" & !is.na(PG.MaxLFQ))$Protein.Group,    
  EP_F4 = filter(Orgnoid_v5, SampleName == "EP_F4" & !is.na(PG.MaxLFQ))$Protein.Group    
 
)

# Run UpSetR
# 2. Convert to input for UpSetR
upset_data <- fromList(list_upset)

# 3. Plot
upset(upset_data, order.by = "freq", nsets = 4, nintersects =20)
upset(upset_data, order.by = "freq", nsets = )


list_upset <- list(
    
  WO_F1 = filter(Orgnoid_v5, SampleName == "WO_F1" & !is.na(PG.MaxLFQ))$Protein.Group,    
  WO_F2 = filter(Orgnoid_v5, SampleName == "WO_F2" & !is.na(PG.MaxLFQ))$Protein.Group,    
  WO_F3 = filter(Orgnoid_v5, SampleName == "WO_F3" & !is.na(PG.MaxLFQ))$Protein.Group,    
  WO_F4 = filter(Orgnoid_v5, SampleName == "WO_F4" & !is.na(PG.MaxLFQ))$Protein.Group    
)

# Run UpSetR
# 2. Convert to input for UpSetR
upset_data <- fromList(list_upset)

# 3. Plot
upset(upset_data, order.by = "freq", nsets = 4, nintersects =20)

#==========================================================================
# this gets the proteins unique to EP and to WO
# EP-only proteins
ep_only <- setdiff(ep, wo)

# WO-only proteins
wo_only <- setdiff(wo, ep)

# Combine into a tidy table
unique_tbl <- tibble::tibble(
  Protein.Group = c(ep_only, wo_only),
  Unique_to = c(
    rep("EP", length(ep_only)),
    rep("WO", length(wo_only))
  )
)

unique_tbl


# 1) Define unique sets
ep_only <- setdiff(ep, wo)
wo_only <- setdiff(wo, ep)

unique_tbl <- tibble::tibble(
  Protein.Group = c(ep_only, wo_only),
  Unique_to     = c(rep("EP", length(ep_only)), rep("WO", length(wo_only)))
)

# 2) Add back abundance rows where observed (many-to-many)
unique_long <- unique_tbl %>%
  left_join(
    Orgnoid_v5 %>%
      filter(!is.na(PG.MaxLFQ)) %>%  # only observed
      filter(!PG.MaxLFQ==0)|>
      mutate(Group = case_when(
        str_detect(SampleName, "^EP_") ~ "EP",
        str_detect(SampleName, "^WO_") ~ "WO",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(Group)) %>%
      select(Protein.Group,Genes, SampleName, Group, CellType, Fraction, PG.MaxLFQ),
    by = "Protein.Group",
    relationship = "many-to-many"
  ) %>%
  # keep only the rows that match the uniqueness label (optional but usually desired)
  filter(Group == Unique_to)

unique_long

unique_long2<-unique_long|>
  filter(!PG.MaxLFQ==0)
  
#$+++++++++++++++++++++

###Plots of Overal Abundance across all samples
ggplot(unique_long, aes(x= SampleName, y=log10(PG.MaxLFQ), fill= Fraction))+
  geom_boxplot(position = position_dodge())+
  facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)

ggplot(unique_long, aes(x= SampleName, y=log10(PG.MaxLFQ), fill= Fraction))+
  geom_violin(position = position_dodge())+
  facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)

ggplot(unique_long, aes(x= CellType, y=log10(PG.MaxLFQ), fill= CellType))+
  geom_boxplot(position = position_dodge())+
 # facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)


ggplot(unique_long, aes(x= CellType, y=log10(PG.MaxLFQ), fill= CellType))+
  geom_violin(position = position_dodge())+
  # facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)

ggplot(Organoids_v5_ppm, aes(x= CellType, y=log10(abundance), fill= CellType))+
  geom_boxplot(position = position_dodge())+
  # facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)

ggplot(Organoids_v2, aes(x= SampleName, y=log10(PG.MaxLFQ), fill= Fraction))+
  geom_boxplot(position = position_dodge())+
   facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)

ggplot(Organoids_v2, aes(x= CellType, y=log10(PG.MaxLFQ)))+
  geom_boxplot(position = position_dodge())+
 # facet_wrap(~CellType, ncol=2, scales = "free_x")+
  theme_bw(base_size = 12)