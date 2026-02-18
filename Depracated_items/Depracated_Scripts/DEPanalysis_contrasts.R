
###Volcano plot from DEP code
print(as.data.frame(colData(data_imp)))


#Tests for RDH12 vs GFP recovery experiment
data_diff <- test_diff(data_imp, type = "manual", 
                       test = c("RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr",
                                "RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr",
                                "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                                "GFP_200_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" ,
                                "GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr"))
dep <- add_rejections(data_diff, alpha = 0.054, lfc = log2(2))

plot_volcano(dep, contrast = "RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr", label_size = 2, add_names = TRUE)

plot_volcano(dep, contrast = "RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr", label_size = 2, add_names = TRUE)

data_results <- get_results(dep)

RDHvsGFP_DEP_MM_v2<-data_results

##################################
##TEST for RDH12-expressing vs Control 

data_diff_acute <- test_diff(data_imp_acute, type = "manual", 
                       test = c("RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
                                "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr",
                                "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr"
                                ))
dep <- add_rejections(data_diff_acute, alpha = 1, lfc = log2(1.5))

plot_volcano(dep, contrast = "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr", label_size = 2, add_names = TRUE)


data_results <- get_results(dep)

RDH12vsEtOH_DEPimputecorrect_MM_v1<-data_results


######if all contrast together
data_diff <- test_diff(data_imp, type = "manual", 
                       test = c("RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr",
                                "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr",
                                "RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr",
                                "RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr",
                                "RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr",
                                "RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                                "RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr",
                                "GFP_200_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr" ,
                                "GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr"
                       ))

dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1.5))
plot_volcano(dep, contrast = "RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 2, add_names = TRUE)

data_results <- get_results(dep)

DEP_mixedimpute_corrected_MM_v1 <-data_results
