data_diff <- test_diff(
  data_imp, 
  type = "manual",
  test = c(
    # treatment effect within RDH12
    "RDH12_100_vs_RDH12_Veh",
    "RDH12_200_vs_RDH12_Veh",
    "RDH12_200_vs_RDH12_100"
  
  )
)

dea <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1.5))

plot_pca(dea, x=1, y=2, n=248, point_size = 4)

##Plot initial volcano plots
plot_volcano(dea, contrast = "RDH12_100_vs_RDH12_Veh", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12_100_vs_RDH12_Veh")

plot_volcano(dea, contrast = "RDH12_200_vs_RDH12_Veh", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12_200_vs_RDH12_Veh")

plot_volcano(dea, contrast = "RDH12_200_vs_RDH12_100", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12_200_vs_RDH12_100")



plot_single(dea, proteins = "LPC 18:1", type = "centered") +
  labs(title = "LPC 18:1", x = "Condition")

plot_single(dea, proteins = "DG 34:1 NL 20:0", type = "centered") +
  labs(title = "DG 34:1 NL 20:0", x = "Condition")


# Generate a results table
data_results <- get_results(dea)
write.csv(data_results,"Lipidomics/output_txts/DEA_output_v1.csv")

data_results_long<- get_df_long(dea)

data_results_wide<- get_df_wide(dea)

# Plot the Pearson correlation matrix
plot_cor(dea, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


###Plot centered  log2intensity heatmaps

png(  filename = file.path(path_figures, "Heatmap_log2CentInt.png"),
      width = 8, height = 8, units = "in", res = 1200)

DEP::plot_heatmap(
  dea,
  type = "centered",
  kmeans = TRUE,
  k = 2,
  col_limit = 6,
  show_row_names = TRUE,
  row_font_size = 1,
  indicate = c("condition", "replicate")
)

dev.off()


library(pheatmap)
library(RColorBrewer)

mat_centered <- data_results %>%
  select(name, RDH12_Veh_centered, RDH12_100_centered, RDH12_200_centered) %>%
  distinct(name, .keep_all = TRUE) %>%
  tibble::column_to_rownames("name") %>%
  as.matrix()

pheatmap(
  mat_centered,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  breaks = seq(-max(abs(mat_centered), na.rm = TRUE),
               max(abs(mat_centered), na.rm = TRUE),
               length.out = 101),
  show_rownames = TRUE,
  fontsize_row = 5,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA
)

colnames(data_results)
##log2FC
mat_centered <- data_results %>%
  select(name, RDH12_100_vs_RDH12_Veh_ratio,RDH12_200_vs_RDH12_Veh_ratio, RDH12_200_vs_RDH12_100_ratio) %>%
  distinct(name, .keep_all = TRUE) %>%
  tibble::column_to_rownames("name") %>%
  as.matrix()

pheatmap(
  mat_centered,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  breaks = seq(-max(abs(mat_centered), na.rm = TRUE),
               max(abs(mat_centered), na.rm = TRUE),
               length.out = 101),
  show_rownames = TRUE,
  fontsize_row = 5,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA
)
