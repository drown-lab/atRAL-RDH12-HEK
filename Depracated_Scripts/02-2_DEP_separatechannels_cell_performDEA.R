

###Perform Limma DEP from Bioconductor workflow and visualize initial results
#https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html

#treatment effect 24h
X24h_MG132_vs_X24h_UT_mglps
#treatment effect 6h
X6h_MG132_vs_X6h_UT_mglps
#time effect within UT
X24h_UT_mglps_vs_X6h_UT_mglps
#time effect with MG132
X24h_MG132_vs_X6h_MG132


# Test all possible comparisons of samples
data_H_diff <- test_diff(data_H_imp_mixed, type = "manual",
                         test = c("X24h_MG132_vs_X24h_UT_mglps", "X6h_MG132_vs_X6h_UT_mglps","X24h_UT_mglps_vs_X6h_UT_mglps", "X24h_MG132_vs_X6h_MG132" ))
dep_H <- add_rejections(data_H_diff, alpha = 0.055, lfc = log2(1.5))

data_L_diff <- test_diff(data_L_imp_mixed, type = "manual",
                         test = c("X24h_MG132_vs_X24h_UT_mglps", "X6h_MG132_vs_X6h_UT_mglps","X24h_UT_mglps_vs_X6h_UT_mglps", "X24h_MG132_vs_X6h_MG132" ))
dep_L <- add_rejections(data_L_diff, alpha = 0.055, lfc = log2(1.5))

##Plot initial volcano plots
plot_volcano(dep_H, contrast = "X24h_MG132_vs_X24h_UT_mglps", label_size = 3, add_names = TRUE) +
  labs(title = "Heavy channel - X24h_MG132_vs_X24h_UT_mglps")
#ggsave("figures_separate_DEP/Volcano_Heavy_WvsB_D7.png", width = 5, height = 5)
plot_volcano(dep_L, contrast = "X24h_MG132_vs_X24h_UT_mglps", label_size = 3, add_names = TRUE) +
  labs(title = "Light channel - X24h_MG132_vs_X24h_UT_mglps")

plot_volcano(dep_H, contrast = "X24h_UT_mglps_vs_X6h_UT_mglps", label_size = 3, add_names = TRUE) +
  labs(title = "Heavy channel - X24h_UT_mglps_vs_X6h_UT_mglps")
plot_volcano(dep_L, contrast = "X24h_UT_mglps_vs_X6h_UT_mglps", label_size = 3, add_names = TRUE) +
  labs(title = "Light channel - X24h_UT_mglps_vs_X6h_UT_mglps")
#ggsave("figures_separate_DEP/Volcano_Light_WvsB_D7.png", width = 5, height = 5)

#plot centered protein abudance for individual proteins
plot_single(dep_H, proteins = "HSPA6", type = "centered") +
  labs(title = "Heavy channel", x = "Condition")
plot_single(dep_L, proteins = "HSPA6", type = "centered") +
  labs(title = "Light channel", x = "Condition")

plot_single(dep_H, proteins = "HSPA5", type = "centered") +
  labs(title = "Heavy channel", x = "Condition")
plot_single(dep_L, proteins = "HSPA5", type = "centered") +
  labs(title = "Light channel", x = "Condition")


#plot actual fold change for indivudal proteins
plot_single(dep_H, proteins = c("HSPA5", "HSPA6"))+
  labs(title = "Heavy Channel")
plot_single(dep_L, proteins = c("HSPA5", "HSPA6"))+
  labs(title = "Light Channel")



# Generate a results table
data_H_results <- get_results(dep_H)
data_H_results$channel <- "H"
data_L_results <- get_results(dep_L)
data_L_results$channel <- "L"
data_results <- rbind(data_L_results, data_H_results)

# Plot the first and second principal components
plot_pca(dep_H, x = 1, y = 2, n = 500, point_size = 4)
plot_pca(dep_L,x = 1, y = 2, n = 500, point_size = 4)

###Plot centered  log2intensity heatmaps
png(filename = "Cell_Perturbed/Figures/Heatmap_Heavy.png",
    width = 6, height = 6, units = "in", res = 1200)

DEP::plot_heatmap(
  dep_H,
  type = "centered",
  kmeans = TRUE,
  k = 6,
  col_limit = 6,
  show_row_names = TRUE,
  row_font_size = 1,
  indicate = c("condition", "replicate")
)

dev.off()


png(filename = "Cell_Perturbed/Figures/Heatmap_Light.png",
    width = 6, height = 6, units = "in", res = 1200)

DEP::plot_heatmap(
  dep_L,
  type = "centered",
  kmeans = TRUE,
  k = 6,
  col_limit = 6,
  show_row_names = TRUE,
  row_font_size = 1,
  indicate = c("condition", "replicate")
)

dev.off()
