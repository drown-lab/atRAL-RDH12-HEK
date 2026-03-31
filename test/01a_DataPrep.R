


#Read In raw data files
SILACHL_v1 <- read_parquet("Cell_HalfLife/quant_DIANN/report.parquet")
exp_designOG<-read.csv("Cell_HalfLife/metadata/HLdata2.csv")

#filter out low confidence proteins
SILACHL_v2<-SILACHL_v1|>
  filter(
  Q.Value <= 0.01,
  PG.Q.Value <= 0.05,
  Lib.Q.Value <= 0.01,
  Lib.PG.Q.Value <= 0.01,
  Channel.Q.Value <= 0.05,
  Proteotypic == 1,
  !grepl("cRAP", Protein.Ids, ignore.case = TRUE)
)

SILACHL_v3 <- SILACHL_v2 |>
  dplyr:: select(Run, Protein.Group, Channel, Genes, PG.MaxLFQ) |>
  filter(Channel !='') |>
  unique()



#remove extra contaminatns keratin protiens as found from doi.org/10.1021/acs.jproteome.2c00145
contam_protein_list <- c(
  #From Keratin
  "Q6E0U4", "P20930", "Q5D862", "Q86YZ3", "P04264", "P13645", "Q99456", "P13646",
  "P02533", "P19012", "P08779", "Q04695", "P05783", "P08727", "P35908", "P35900",
  "Q8N1A0", "Q9C075", "Q2M2I5", "Q7Z3Z0", "Q7Z3Y9", "Q7Z3Y8", "Q7Z3Y7", "P12035",
  "Q15323", "Q14532", "O76009", "Q14525", "O76011", "Q92764", "O76013", "O76014",
  "O76015", "Q6A163", "P19013", "Q6A162", "P13647", "P02538", "P04259", "P48668",
  "Q3KNV1", "P08729", "Q3SY84", "Q14CN4", "Q86Y46", "Q7RTS7", "O95678", "Q01546",
  "Q7Z794", "Q8N1N4", "Q5XKE5", "P05787", "Q6KB66", "Q14533", "Q9NSB4", "P78385",
  "Q9NSB2", "P78386", "O43790", "A6NCN2", "P35527", "P60331", "P60014", "P60412",
  "P60413", "P60368", "P60369", "P60372", "P60370", "P60371", "P60409", "P60410",
  "P60411", "Q07627", "Q8IUC1", "P59990", "P59991", "P60328", "P60329", "Q8IUG1",
  "Q8IUC0", "Q52LG2", "Q3SY46", "Q3LI77", "P0C5Y4", "Q9BYS1", "Q3LI76", "A8MUX0",
  "Q9BYP8", "Q8IUB9", "Q3LHN2", "Q7Z4W3", "Q3LI73", "Q3LI72", "Q3LI70", "Q3SYF9",
  "Q3LI54", "Q3LI63", "Q3LI61", "Q9BYU5", "Q3LI58", "Q3LI59", "Q3LHN1", "Q9BYT5",
  "Q3MIV0", "P0C7H8", "A1A580", "Q9BYR9", "Q3LI83", "Q3LHN0", "Q6PEX3", "Q3LI81",
  "A8MX34", "Q9BYR8", "Q9BYR7", "Q9BYR6", "Q9BYQ7", "Q9BYQ6", "Q9BQ66", "Q9BYR5",
  "Q9BYR4", "Q9BYR3", "Q9BYR2", "Q9BYQ5", "Q9BYR0", "Q9BYQ9", "Q9BYQ8", "Q6L8H4",
  "Q6L8G5", "Q6L8G4", "Q701N4", "Q6L8H2", "Q6L8H1", "Q701N2", "Q6L8G9", "Q6L8G8",
  "O75690", "P26371", "Q3LI64", "Q3LI66", "Q3LI67", "Q8IUC3", "Q8IUC2", "A8MXZ3",
  "Q9BYQ4", "Q9BYQ3", "Q9BYQ2", "A8MVA2", "A8MTY7", "Q9BYQ0", "Q9BYP9",
  #skin
  "P81605","Q08554","Q15517","P31151","Q9NZT1","Q96DR8",
  #From FBS
  "P01023","P02774","P00734","P01024","P02771","P04114","Q06033","P25311", "P00747"
  
)

SILACHL_v4 <- SILACHL_v3 %>%
  filter(!Protein.Group %in% contam_protein_list)
print(unique(SILACHL_v4$Run))

SILACHL_v5 <- SILACHL_v4 %>%
  left_join(
    exp_designOG,
    by = c("Run" = "Run"))
    
  
colnames(SILACHL_v5)


SILACHL_v6 <- SILACHL_v5 |>
  unique() |>
  pivot_wider(values_from = PG.MaxLFQ, names_from = Channel, values_fill = 0) |>
  mutate(RHIA = H/(L+H), logsumAbund = log(L+H, 2)) 

SILACHL_v6<-SILACHL_v6|>
  mutate(Time = as.numeric(Time)) 


ggplot(SILACHL_v6, aes(x = Time, y = RHIA, fill = Protein.Group)) +
  geom_line(alpha = 0.01, show.legend = FALSE, base_size = 14, base_family= "Arial") +
  theme_bw(base_size = 14)+
  labs(title= "All Proteins RHIA vs Time")
ggsave(filename = file.path(path_figures,"ProteinsRHIAoverTime.svg"), width = 7, height = 5)

SILACHL_v7<-SILACHL_v6|>
  filter(logsumAbund >0)

ggplot(SILACHL_v7, aes(x = Time, y = RHIA, fill = Protein.Group)) +
  geom_line(alpha = 0.01, show.legend = FALSE, base_size = 14) +
  theme_bw(base_size = 12)+
  labs(title= "All Proteins RHIA vs Time")
ggsave(filename = file.path(path_figures,"ProteinsRHIAoverTime.svg"), width = 7, height = 5)

#call to df to run source scripts
df<-SILACHL_v7
colnames(df)

#Run to produce Presence Absence Figures
source("Cell_HalfLife/Script_HalfLife/01b_Missingness_PresenceAbsence.R")

#Run to produce Missingess Summarization
source("Cell_HalfLife/Script_HalfLife/01c_Missingness_Summarization_Xof3.R")

#Run to plot Missingness Abundance Plot
source("Cell_HalfLife/Script_HalfLife/01d_Missingness_SumAbundancePlot.R")



instability <- SILACHL_v7 |>
  filter(!is.na(RHIA)) |>
  group_by(Protein.Group) |>
  arrange(Time) |>
  summarize(quant.days = n(),
            pathlength = pathlength(Time, RHIA),
            swings = has.swings(RHIA),
            deltaRIA = total.RHIA.changes(RHIA),
            extreme_values = length(RHIA[RHIA == 0]) + length(RHIA[RHIA == 1]),#how many 0 and 1s
            extreme_fraction = extreme_values/quant.days
            ) 


instability <- SILACHL_v7 |>
  filter(!is.na(RHIA)) |>
  group_by(Protein.Group) |>
  arrange(Time, .by_group = TRUE) |>
  summarise(
    quant.days = n(),
    deltaRIA = total.RHIA.changes(RHIA),
    extreme_values = sum(RHIA == 0 | RHIA == 1, na.rm = TRUE),
    
    swing_stats(RHIA),
    extreme_fraction = extreme_values/quant.days,
    
    .groups = "drop"
  )

#Number of time points
ggplot(instability, aes(x = quant.days)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  labs(x = "# quantified timepoints", y = "Proteins")

#QC plots looking over potential unstable protein quantification
#Checks Cgane in RHIA over time points sampled
ggplot(instability, aes(x = deltaRIA)) +
  geom_histogram(bins = 50) +
  theme_bw()

ggplot(instability, aes(x = deltaRIA, y= n_neg_swings ,color = n_neg_swings > 0)) +
  geom_point(alpha = 0.3) +
  theme_bw()

ggplot(instability, aes(x = deltaRIA, y= n_neg_swings ,color = n_neg_swings ==1)) +
  geom_point(alpha = 0.3) +
  theme_bw()
ggplot(instability, aes(x = n_swings, y= n_neg_swings ,color = n_neg_swings ==1)) +
  geom_point(alpha = 0.3) +
  theme_bw()

ggplot(instability, aes(x = deltaRIA, y= n_swings ,color = n_swings >=2)) +
  geom_point(alpha = 0.3) +
  theme_bw()

ggplot(instability, aes(x = deltaRIA, y= n_swings ,color = n_neg_swings ==1)) +
  geom_point(alpha = 0.3) +
  theme_bw()

ggplot(instability, aes(x = quant.days, color = n_swings)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  theme_bw()

ggplot(instability, aes(x = quant.days, y = deltaRIA)) +
  geom_point(alpha = 0.3) +
  theme_bw()

ggplot(instability, aes(x = extreme_fraction)) +
  geom_histogram(bins = 50) +
  theme_bw()

ggplot(instability, aes(x = n_swings)) +
  geom_bar() +
  theme_bw()
ggplot(instability, aes(x = n_neg_swings)) +
  geom_bar() +
  theme_bw()

ggplot(instability, aes(x = deltaRIA, fill = n_swings)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  theme_bw()

ggplot(instability, aes(x = extreme_values)) +
  geom_histogram(binwidth = 1) +
  theme_bw()
library(GGally)

GGally::ggpairs(
  instability %>%
    select(quant.days, deltaRIA, extreme_fraction, n_swings, n_neg_swings)
)
#Get Suspect POIs
#gets proteins with <3 identifications across samples, and looks for proteins with 3 swings and 1 neg swing
SuspectProts <- instability |>
  filter(
    quant.days < 3 |
      (n_swings == 3 & n_neg_swings == 1)
  )

#filter out suspect proteins
SILACHL_v8 <- SILACHL_v7 |>
  filter(!Protein.Group %in% SuspectProts$Protein.Group)

ggplot(SILACHL_v8, aes(x = Time, y = RHIA, fill = Protein.Group)) +
  geom_line(alpha = 0.01, show.legend = FALSE, base_size = 14) +
  theme_bw(base_size = 12)+
  labs(title= "All Proteins RHIA vs Time")

write.csv(SILACHL_v8, "Cell_HalfLife/output_txts/SILAC_filtered_rawdata.csv")

#################
#Plot Bip
bip <- SILACHL_v8 |>
  filter(Protein.Group == 'P11021')

#Drome Bip P29844
#human Bip P11021

ggplot(bip, aes(x = Time, y = RHIA)) +
  geom_point() +
  theme_bw()+
  labs(
    title="Bip" 
  ) +
  lims(x = c(0,50), y = c(0,1))

bip.fixed <- minpack.lm::nlsLM(RHIA ~ aa_two_compartment(Time, k0),
                               data = bip,
                               start = list(k0 = 0.05))
summary(bip.fixed)

x <- data.frame(Time = seq(0,50,0.2))
x$RHIA <- predict(bip.fixed, x)


ggplot(bip, aes(x = Time, y = (RHIA))) +
  geom_point() +
  geom_line(data = x) +
  labs(x = "Time (hr)", y = "Heavy Isotope Incorporation") +
  theme_bw(base_size = 15) +
  lims(x = c(0,50), y = c(0,1))
