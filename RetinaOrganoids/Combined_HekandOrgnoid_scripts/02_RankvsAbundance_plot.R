
library(dplyr)
library(ggplot2)
library(ggrepel)

ggplot(combined_proteome3, aes(x=log10_ppm, y= log2(abundance)))+
  geom_point()+
  facet_wrap(~CellType)+
  theme_bw(base_size = 12)



df <- combined_proteome3

# ----------------------------
# Choose which IDs to highlight
# ----------------------------
target_genes <- c("RDH12")    # <-- edit
target_genes <- c("RDH12", "RDH11", "RDH10", "RDH5", "RDH14")
target_genes <- c("RDH12", "RCVRN", "CRX", "NRL", "GNAT1", "GNAT2", "ARR3",
                "SAG")
# target_proteins <- c("P12345", "Q9XYZ1")  # optional if you want Protein.Group matching too

# ----------------------------
# Prep: clean + rank
# ----------------------------
plot_df <- df %>%
  mutate(
    Genes = as.character(Genes),
    Protein.Group = as.character(Protein.Group),
    abundance = as.numeric(abundance),
    log2_abundance = log2(abundance),
    is_target = Genes %in% target_genes # | Protein.Group %in% target_proteins
  ) %>%
  filter(is.finite(log2_abundance)) %>%
  arrange(desc(abundance)) %>%
  mutate(rank = row_number())

plot_df <- df %>%
  mutate(
    Genes = as.character(Genes),
    Protein.Group = as.character(Protein.Group),
    abundance = as.numeric(abundance),
    log2_abundance = log2(abundance),
    is_target = Genes %in% target_genes # | Protein.Group %in% target_proteins
  ) %>%
  filter(is.finite(log2_abundance)) %>%
  arrange(desc(abundance)) %>%
  group_by(CellType)|>
  mutate(rank = row_number())

# ----------------------------
# Scatter plot: rank vs log2 abundance
# ----------------------------
p <- ggplot(plot_df, aes(x = rank, y = log2_abundance)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_point(
    data = plot_df %>% filter(is_target),
    aes(x = rank, y = log2_abundance, color = Genes),
    size = 3
  ) +
  #ggrepel::geom_text_repel(
   # data = plot_df %>% filter(is_target),
    #aes(label = Genes),
   # size = 3,
    #max.overlaps = Inf
 # ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  labs(
    title = "Protein abundance rank plot",
    x = "Rank (highest abundance → lowest)",
    y = "log2(abundance)"
  )+
  facet_wrap(~CellType, ncol = 1)

print(p)

plot_df<-plot_df|>
  filter(!CellType=="GFP")

p <- ggplot(plot_df, aes(x = rank, y = log2_abundance)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_point(
    data = plot_df %>% filter(is_target),
    aes(x = rank, y = log2_abundance, color = Genes),
    size = 3
  ) +
  #ggrepel::geom_text_repel(
  # data = plot_df %>% filter(is_target),
  #aes(label = Genes),
  # size = 3,
  #max.overlaps = Inf
  # ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  labs(
    title = "Protein abundance rank plot",
    x = "Rank (highest abundance → lowest)",
    y = "log2(abundance)"
  )+
  facet_wrap(~CellType, ncol = 1)

print(p)