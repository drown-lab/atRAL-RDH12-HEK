#venn diagram

library(dplyr)
library(tidyr)
library(ggVennDiagram)
library(ggplot2)

df <- combined_proteome3

# Choose the ID you want to define "presence"
# Use Protein.Group (unique protein IDs) OR Genes
id_col <- "Protein.Group"   # change to "Genes" if you prefer gene-level overlap

# Keep proteins that are truly present (optional rule)
df2 <- df %>%
  filter(!is.na(.data[[id_col]]), .data[[id_col]] != "")

# Build list of sets: one vector per CellType
prot_sets <- df2 %>%
  distinct(CellType, id = .data[[id_col]]) %>%
  group_by(CellType) %>%
  summarise(ids = list(unique(id)), .groups = "drop") %>%
  { setNames(.$ids, .$CellType) }

# ---- If you have >4 CellTypes, pick which ones to plot
# Example: keep the 4 CellTypes with most proteins
if (length(prot_sets) > 4) {
  message("More than 4 CellTypes found; keeping the top 4 by set size for Venn.")
  sizes <- sapply(prot_sets, length)
  keep_names <- names(sort(sizes, decreasing = TRUE))[1:4]
  prot_sets <- prot_sets[keep_names]
}


library(ggvenn)
#FFFF00
#984EA3

cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")

p <- ggvenn(
  prot_sets,
  fill_color = cols,
  fill_alpha = 0.5,
  stroke_size = 0.8,
  set_name_size = 4
) +
  theme_void()

print(p)

