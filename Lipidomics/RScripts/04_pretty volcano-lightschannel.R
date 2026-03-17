library(ggrepel)  # for nice text labels
library(dplyr)
library(tidyr)
library(stringr)

#This workflow requires some knowledge about how to work the script
#Please review code carefully

#=======================================================================
#this section is if you are using the DEP output procuded as described in next two lines
#First read in DEP output file
#this is the file produced in `04_Data_results_with_missing_class_output`
#this file has many columns name, ID, p.val, p.adj, ratio, significant, channel, missclass
Paper1_DEPresults_v1<-read.csv("Cell_Perturbed/Paper1/Paper1_DEPresults_v1.csv")
colnames(Paper1_DEPresults_v1)


#  Pivot p-values + ratios into long format
dep_long <- Paper1_DEPresults_v1 %>%
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

# Clean missing class (convert to long as well)
miss_long <- Paper1_DEPresults_v1 %>%
  select(ID, channel, starts_with("missclass_")) %>%
  pivot_longer(
    cols = starts_with("missclass_"),
    names_to = "condition",
    names_prefix = "missclass_",
    values_to = "missing_class"
  )

# Join missing class onto tidy DEP table that is formated correctly for volcano plot making
dep_long_final <- dep_long %>%
  left_join(miss_long, by = c("ID", "channel"))

dep_long_final

colnames(dep_long_final)
dep_long_final2<-dep_long_final|>
  select(channel, Gene, contrast, p.val, p.adj, condition, ratio, missing_class)|>
  unique()
dep_long_final3<-dep_long_final2|>
  select(channel, Gene, contrast, p.val, p.adj, ratio)|>
  unique()

#write a csv file of the volcano friednly table
write.csv(dep_long_final3, "Cell_Perturbed/Paper1/output_txts/Volcanotable_v1.csv")

#=================================================================================
#This code below is for strucutring your volcano plot

#User Input needed:
# 1: define your contrast and optional favorite proteins
contrast <- "X24h_MG132_vs_X6h_MG132" #need to define contrast because normally have multiple contrast tested in DEP code

#2: if you have some favorite Gene names you want to visualize call them here
favorite_IDs <- c("HTT", "ATF4", "HSPA5", "SNCA", 
                  "HERPUD1" 
                  )  # example

#This will format the table based on channel you want and contrast you want
res_longH <- dep_long_final3 %>%
  filter(channel == "L") %>% #need to call channel of interest
  filter(contrast =="X24h_MG132_vs_X6h_MG132")|> #need to call contrast of interest
  transmute(
    Gene,
    channel,
    ratio = ratio,
    padj  = p.adj,
    sig   = padj < 0.054, #set signifcance thresholds
    direction = case_when(
      ratio >  log2(2) & sig ~ "Up", #set signifcance thresholds
      ratio < -log2(2) & sig ~ "Down", #set signifcance thresholds
      TRUE ~ "NS"
    )
  )

# pick top 10 most significant to be labeled other than favorited Ids; can adjust value
top10 <- res_longH %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  pull(Gene)

# merge top10 + favorites for labeling
res_longH <- res_longH %>%
  mutate(label = if_else(Gene %in% c(top10, favorite_IDs), Gene, NA_character_))

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
  geom_vline(xintercept = c(-log2(2), log2(2)),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.054),
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