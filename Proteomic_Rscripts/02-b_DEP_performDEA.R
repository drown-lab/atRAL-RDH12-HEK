

###Perform Limma DEP from Bioconductor workflow and visualize initial results
#https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html

as.data.frame(colData(data_imp_mixed)) |> 
  dplyr::pull(condition) |> 
  unique()

#treatment effect acute
RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr
RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr
RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr

#treatment effect control recovery
GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr
GFP_200_atRAL5hr.24h_recv_vs_GFP_control_atRAL5hr.24h_recvr
GFP_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr

#treatment effect RDH12 recovery
RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr
RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr
RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr.24h_recvr


#recovery effect RDH12
RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr
RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_200_atRAL5hr
RDH12_control_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr

#genotype effect
RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr
RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr
RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr


data_imp_mixed




# Test all possible comparisons of samples
data_diff <- test_diff(data_imp_mixed, type = "manual",
                         test = c(#treatment effect acute
                           "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
                           "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr",
                           "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
                           
                           #treatment effect control recovery
                           "GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr",
                           "GFP_200_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr",
                           "GFP_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                           
                           #treatment effect RDH12 recovery
                           "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                           "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                           "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr.24h_recvr",
                           
                           
                           #recovery effect RDH12
                           "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr",
                           "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_200_atRAL5hr",
                           "RDH12_control_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr",
                           
                           #genotype effect
                           "RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                           "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr",
                           "RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" ))

dep <- add_rejections(data_diff, alpha = 0.055, lfc = log2(1.5))

##Plot initial volcano plots
plot_volcano(dep, contrast = "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr")
#ggsave("figures_separate_DEP/Volcano_Heavy_WvsB_D7.png", width = 5, height = 5)
plot_volcano(dep, contrast = "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 3, add_names = TRUE) +
  labs(title = "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr")

#ggsave("figures_separate_DEP/Volcano_Light_WvsB_D7.png", width = 5, height = 5)

#plot centered protein abudance for individual proteins
plot_single(dep, proteins = "HSPA6", type = "centered") +
  labs(title = "HSPA6", x = "Condition")


plot_single(dep, proteins = "HSPA5", type = "centered") +
  labs(title = "HSPA5", x = "Condition")



#plot actual fold change for indivudal proteins
plot_single(dep, proteins = c("HSPA5", "HSPA6"))+
  labs(title = "Channel")
plot_single(dep, proteins = c("RDH12", "RDH11"))+
  labs(title = "Channel")


# Generate a results table
data_results <- get_results(dep)


# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

###Plot centered  log2intensity heatmaps
png(filename = "Proteomic_Figs/Heatmap.png",
    width = 6, height = 8, units = "in", res = 1200)

DEP::plot_heatmap(
  dep,
  type = "centered",
  kmeans = TRUE,
  k = 4,
  col_limit = 6,
  show_row_names = TRUE,
  row_font_size = 1,
  indicate = c("condition", "replicate")
)

dev.off()


