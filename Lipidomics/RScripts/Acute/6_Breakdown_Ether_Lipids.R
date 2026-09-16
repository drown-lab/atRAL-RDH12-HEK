
##Get Ether Lipids
library(dplyr)
library(stringr)
library(ggplot2)

df<-read.csv("Lipidomics/output_txts/Acute_filteredlipidtable.csv")

df <- df %>%
  mutate(
    ether = str_extract(
      picked_candidate, "O-",
     
    )
  )
df <- df %>%
  mutate(lipid_class_ether = str_extract(picked_candidate,"PC O-|PE O-|PI O-|DG O-|PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM" )
         )
write.csv(df,"Lipidomics/output_txts/Acute_filteredlipidtable.csv" )

Acute100_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s4, s5, s6)))

Acute100_v2<-Acute100_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Acute100_v3<-Acute100_v2|>
  filter(max_divid_blank >= 1.3)

Acute100_v4<-Acute100_v3|>
  unique()

##Plot Distribution of Lipids
ggplot(Acute100_v2, aes(x= max_divid_blank, fill = lipid_class1))+
  geom_histogram(bins = 40)+
  geom_vline(xintercept = 1.3, linetype = "dashed")+
  theme_bw(base_size=12)+
  scale_x_continuous(trans = "log2")+
  labs(
    title = "Acute: Distribution of MRMs with Signal-to-Blank",
    subtitle = "Dashed line: 1.3x "
  )+
  scale_y_continuous(trans = "log10")

Acute100_filtered_summary <- Acute100_v4  %>%
  group_by( lipid_class_ether) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

Acute100_filtered_summary2 <- Acute100_v4  %>%
  group_by( lipid_class_ether) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Acute100_filtered_summary2 <- Acute100_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether = factor(lipid_class_ether, levels = lipid_class_ether))

Acute200_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s7, s8, s9)))

Acute200_v2<-Acute200_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
Acute200_v3<-Acute200_v2|>
  filter(max_divid_blank >= 1.3)

Acute200_v4<-Acute200_v3|>
  unique()

Acute200_filtered_summary <- Acute200_v4  %>%
  group_by( lipid_class_ether) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

Acute200_filtered_summary2 <- Acute200_v4  %>%
  group_by( lipid_class_ether) %>%
  summarise(distinct_precursor = n_distinct(precursor))

Acute200_filtered_summary2 <- Acute200_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether = factor(lipid_class_ether, levels = lipid_class_ether))
##########+===================
##########+

AcuteVeh_v1<-df|>
  rowwise() |>
  mutate(maxvalue = max(c(s1, s2, s3)))

AcuteVeh_v2<-AcuteVeh_v1|>
  mutate(
    max_divid_blank = maxvalue/blank
  )
#Filter by 30% higher than the blank
#Standard cut-off 
AcuteVeh_v3<-AcuteVeh_v2|>
  filter(max_divid_blank >= 1.3)

AcuteVeh_v4<-AcuteVeh_v3|>
  unique()

AcuteVeh_filtered_summary <- AcuteVeh_v4  %>%
  group_by( lipid_class_ether) %>%
  summarise(distinct_mrms = n_distinct(mrm1))

AcuteVeh_filtered_summary2 <- AcuteVeh_v4  %>%
  group_by( lipid_class_ether) %>%
  summarise(distinct_precursor = n_distinct(precursor))

AcuteVeh_filtered_summary2 <- AcuteVeh_filtered_summary2 %>%
  arrange(desc(distinct_precursor)) %>%
  mutate(lipid_class_ether = factor(lipid_class_ether, levels = lipid_class_ether))  

write.csv(AcuteVeh_filtered_summary2, "Lipidomics/output_txts/AcuteVeh_LipidClassEther_summary.csv")
write.csv(Acute100_filtered_summary2,"Lipidomics/output_txts/Acute100_LipidClassEther_summary.csv")
write.csv(Acute200_filtered_summary2, "Lipidomics/output_txts/Acute200_LipidClassEther_summary.csv") 

Acute100_v4
Acute200_v4
colnames(AcuteVeh_v4)


