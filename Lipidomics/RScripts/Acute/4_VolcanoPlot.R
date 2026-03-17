
library(ggrepel)  # for nice text labels
library(dplyr)
library(tidyr)
library(stringr)


colnames(data_results_v1)

data_results_v1 <- data_results |>
  mutate(lipid_class1 = str_extract(name, "PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM")) |> # add all desired abbreviations
  mutate(lipid_class = str_extract(name, "LPC|LPE|PC|TG|PE|PS|PI|PG|CAR|CE|Cer|DG|FA|SM"))


#  Pivot p-values + ratios into long format
dea_long <- data_results_v1 %>%
  pivot_longer(
    cols = matches("(_p\\.val$|_p\\.adj$|_ratio$)"),
    names_to = c("contrast", "stat"),
    names_pattern = "(.*)_(p\\.val|p\\.adj|ratio)$",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )

dea_long_v1<-dea_long|>
  select(name, ID, lipid_class1, lipid_class, contrast, p.val, p.adj, ratio)

write.csv(dea_long_v1, "Lipidomics/output_txts/Volcanotable_Acute_v1.csv")

print(unique(dea_long_v1$contrast))
#This code below is for strucutring your volcano plot
#User Input needed:
# 1: define your contrast and optional favorite proteins
contrast <- "RDH12_200_vs_RDH12_Veh" #need to define contrast because normally have multiple contrast tested in DEP code

#2: if you have some favorite name names you want to visualize call them here
favorite_IDs <- c("DG 34:1 NL 18:1", "PC 32:0", "DG 32:0 NL 16:0", 
                   "PC 28:0" 
)  # example

#This will format the table based on  you want and contrast you want
res_longH <- dea_long_v1 %>%
  filter(contrast =="RDH12_200_vs_RDH12_Veh")|> #need to call contrast of interest
  transmute(
    name,
    ratio = ratio,
    padj  = p.val,
    sig   = padj < 0.1, #set signifcance thresholds
    direction = case_when(
      ratio >  log2(1.3) & sig ~ "Up", #set signifcance thresholds
      ratio < -log2(1.3) & sig ~ "Down", #set signifcance thresholds
      TRUE ~ "NS"
    )
  )

# pick top 10 most significant to be labeled other than favorited Ids; can adjust value
top10 <- res_longH %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  pull(name)

# merge top10 + favorites for labeling
res_longH <- res_longH %>%
  mutate(label = if_else(name %in% c(top10, favorite_IDs), name, NA_character_))

##Fancy volcano
ggplot(res_longH, aes(x = ratio, y = -log10(padj))) +
  # base layer: all points
  geom_point(aes(color = direction), alpha = 0.6, size = 2) +
  
  # outline for labeled points
  geom_point(
    data = subset(res_longH, !is.na(label)),
    shape = 21, size = 2.8, stroke = 0.5,
    color = "black", fill = NA
  ) +
  
  # significance cutoffs
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed", color = "grey40") +
  
  # labels with connecting lines
  geom_text_repel(
    aes(label = label),
    max.overlaps = 25,
    size = 3.5,
    color = "black",
    segment.color = "grey20",
    segment.size = 0.3, #line thickness
    box.padding = 0.4,#space around label
    point.padding = 0.3, #label and point space
    force = 2, #fix overlap
    min.segment.length = 0 
  ) +
  
  theme_bw() +
  labs(
    title = paste("Volcano:", contrast),
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  scale_color_manual(
    values = c("Up" = "purple", "Down" = "#00B8B8", "NS" = "grey70"),
    name = "Regulation"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text( size = 12),
    axis.title = element_text(size = 12)
  )