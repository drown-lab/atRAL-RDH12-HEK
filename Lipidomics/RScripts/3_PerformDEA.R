data_diff <- test_diff(
  data_normvsn,
  type = "manual",
  test = c(
    # treatment effect within Control
    "Control_100_vs_Control_Veh",
    "Control_200_vs_Control_Veh",
    "Control_200_vs_Control_100",
    
    # treatment effect within RDH12
    "RDH12._100_vs_RDH12._Veh",
    "RDH12._200_vs_RDH12._Veh",
    "RDH12._200_vs_RDH12._100",
    
    # genotype effect at each treatment
    "RDH12._Veh_vs_Control_Veh",
    "RDH12._100_vs_Control_100",
    "RDH12._200_vs_Control_200"
  )
)

dea <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1.5))

plot_pca(dea, x=1, y=2, n=363, point_size = 4)

##Plot initial volcano plots
plot_volcano(dea, contrast = "RDH12._200_vs_RDH12._Veh", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12._200_vs_RDH12._Veh")

plot_volcano(dea, contrast = "Control_200_vs_Control_Veh", label_size = 3, add_names = TRUE) +
  labs(title = " Control_200_vs_Control_Veh")

plot_volcano(dea, contrast = "RDH12._Veh_vs_Control_Veh", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12._Veh_vs_Control_Veh")


plot_volcano(dea, contrast = "RDH12._200_vs_Control_200", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12._200_vs_Control_200")


plot_single(dea, proteins = "PC 34:5", type = "centered") +
labs(title = "PC 34:5", x = "Condition")

# Generate a results table
data_results <- get_results(dea)
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
  k = 4,
  col_limit = 6,
  show_row_names = TRUE,
  row_font_size = 1,
  indicate = c("condition", "replicate")
)

dev.off()


