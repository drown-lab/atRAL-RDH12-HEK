library(ggrepel)  # for nice text labels
#"RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
#"RDH12_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
#"RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr"))
#"RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
#"RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr",
#"RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
#"RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
#"RDH12_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
#"RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr",
#"RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr",
#"RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
#"RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
#"GFP_200_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" ,
#"GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr"



# define your contrast and optional favorite proteins
contrast <- "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr"
favorite_IDs <- c(
  "RDH10", "RDH11", "RDH12", "RDH13", "RDH14"
)  # example
##inpute table of interest
tableofinterest<-data_results   

colnames(tableofinterest)
res_longH <- tableofinterest    %>%
  transmute(
    name,
    ratio = !!sym(paste0(contrast, "_ratio")),
    padj  = !!sym(paste0(contrast, "_p.val")),
    #padj  = !!sym(paste0(contrast, "_p.adj")),
    
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
  pull(name)

# merge top10 + favorites for labeling
res_longH <- res_longH %>%
  mutate(label = if_else(name %in% c(top10, favorite_IDs), name, NA_character_))

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
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


############################################################
# ---- pick top 5 Up and top 5 Down by padj ----
top_up <- res_longH %>%
  filter(direction == "Up") %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  pull(name)

top_down <- res_longH %>%
  filter(direction == "Down") %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  pull(name)

# ---- label top 5 up + top 5 down + any favorites ----
res_plot <- res_longH %>%
  mutate(
    label = if_else(
      name %in% c(top_up, top_down, favorite_IDs),
      name,
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
    y = "-log10(p-value)"
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


