library(dplyr)
library(tidyr)
library(UpSetR)

# ---------------------------------------------------------
# Choose what "present in a condition" means:
#   present_cutoff = 2  -> present if 2/3 or 3/3 (recommended)
#   present_cutoff = 3  -> present only if 3/3
#   present_cutoff = 1  -> present if 1/3 or more
# ---------------------------------------------------------
present_cutoff <- 1

# ---------------------------------------------------------
# 1) Make a protein x condition binary matrix (1 = present, 0 = absent)
#    Using det_by_cond that you already created
# ---------------------------------------------------------
upset_wide <- det_by_cond |>
  mutate(present = as.integer(n_detect >= present_cutoff)) |>
  select(protein, condition, present) |>
  tidyr::pivot_wider(
    names_from  = condition,
    values_from = present,
    values_fill = 0
  )

# UpSetR wants a plain data.frame; rownames optional but nice
upset_df <- as.data.frame(upset_wide)
rownames(upset_df) <- upset_df$protein
upset_df$protein <- NULL

# ---------------------------------------------------------
# 2) Draw the UpSet plot
# ---------------------------------------------------------
UpSetR::upset(
  upset_df,
  nsets = ncol(upset_df),         # show all conditions
  nintersects = 30,               # top 30 intersections (adjust)
  order.by = "freq",
  decreasing = TRUE,
  keep.order = TRUE,              # keep condition column order as in upset_df
  sets.bar.color = "grey30",
  main.bar.color = "grey10",
  matrix.color = "grey30",
  text.scale = c(1.4, 1.4, 1.1, 1.1, 1.3, 1.2)
)
