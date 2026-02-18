

###Devender Analysis for RDH12 5hr 200vsEtOH
DE_RDH12_5h_200_vs_Control<-read_csv("Z:/data/Projects/Rams_Collab_RDH12experiments/atRAL_experiments/proteomics/DEPandEnrichment_DevenderAnalysis/rdh12_5h/DE_RDH12_5h_200_vs_Control.csv")
colnames(DE_RDH12_5h_200_vs_Control)

gene_mapping<-read_xlsx("Z:/data/Projects/Rams_Collab_RDH12experiments/atRAL_experiments/proteomics/DEPandEnrichment_DevenderAnalysis/rdh12_5h/geneidmapping_2025_12_02.xlsx")

colnames(gene_mapping)


gene_mapping1<-gene_mapping|>
  mutate(
    ProteinID=Entry
  )

##GET GeneIDs since these files only have uniprot id
DE_RDH12_5h_200_vs_Control_c2 <- DE_RDH12_5h_200_vs_Control %>%
  left_join(gene_mapping1, by = "ProteinID")
colnames(DE_RDH12_5h_200_vs_Control_c2)

###Get/Clean GeneID
DE_RDH12_5h_200_vs_Control_c2<-DE_RDH12_5h_200_vs_Control_c2|>
  mutate(
    GeneID = str_remove(`Entry Name`, "_HUMAN$")
  )


# define your contrast and optional favorite proteins
contrast <- "RDH12 5hr 200 vs EtOH"
favorite_IDs <- c(
  "RDH10", "RDH11", "RDH12", "RDH13", "RDH14"
)  # example
##inpute table of interest
tableofinterest<-DE_RDH12_5h_200_vs_Control_c2   

colnames(tableofinterest)

res_longH <- tableofinterest    %>%
  transmute(
    GeneID,
    ratio = logFC,
   # padj=  P.Value,
    padj  = adj.P.Val,
    sig   = padj <= 0.055,
    direction = case_when(
      ratio >=  log2(1.5) & sig ~ "Up",
      ratio <= -log2(1.5) & sig ~ "Down",
      TRUE ~ "NS"
    )
  )

# pick top 10 most significant
top10 <- res_longH %>%
  arrange(padj) %>%
  slice_head(n = 12) %>%
  pull(GeneID)

# merge top10 + favorites for labeling
res_longH <- res_longH %>%
  mutate(label = if_else(GeneID %in% c(top10, favorite_IDs), GeneID, NA_character_))

##FAncy volcano
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
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.055),
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
    y = "-log10(adj. p-value)"
  ) +
  scale_color_manual(
    values = c("Up" = "purple", "Down" = "#00B8B8", "NS" = "grey70"),
    name = "Regulation"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text( size = 12),
    axis.title = element_text(size = 12)
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#######################
###################3
##################

# ---- pick top 5 Up and top 5 Down by padj ----
top_up <- res_longH %>%
  filter(direction == "Up") %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  pull(GeneID)

top_down <- res_longH %>%
  filter(direction == "Down") %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  pull(GeneID)

# ---- label top 5 up + top 5 down + any favorites ----
res_plot <- res_longH %>%
  mutate(
    label = if_else(
      GeneID %in% c(top_up, top_down, favorite_IDs),
      GeneID,
      NA_character_
    )
  )

# ---- Volcano plot ----
ggplot(res_plot, aes(x = ratio, y = -log10(padj))) +
  geom_point(aes(color = direction), alpha = 0.6, size = 2) +
  
  geom_point(
    data = subset(res_plot, !is.na(label)),
    shape = 21, size = 2.8, stroke = 0.5,
    color = "black", fill = NA
  ) +
  
  geom_vline(
    xintercept = c(-log2(1.5), log2(1.5)),
    linetype = "dashed", color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.055),
    linetype = "dashed", color = "grey40"
  ) +
  
  geom_text_repel(
    aes(label = label),
    na.rm = TRUE,
    seed = 1,
    max.overlaps = 25,
    max.time = 2,
    box.padding = 0.6,
    point.padding = 0.4,
    force = 8,
    force_pull = 0.2,
    min.segment.length = 0,
    segment.color = "grey20",
    segment.size = 0.3,
    size = 3.5,
    color = "black"
  ) +
  
  theme_bw() +
  labs(
    title = paste("Volcano:", contrast),
    x = "log2 Fold Change",
    y = "-log10(adj. p-value)"
  ) +
  scale_color_manual(
    values = c("Up" = "purple", "Down" = "#00B8B8", "NS" = "grey70"),
    name = "Regulation"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
