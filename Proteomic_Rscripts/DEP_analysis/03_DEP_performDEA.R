

###Perform Limma DEP from Bioconductor workflow and visualize initial results
#https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html

as.data.frame(colData(data_imp_mixed)) |> 
  dplyr::pull(condition) |> 
  unique()

#treatment effect acute
#RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr
#RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr
#RDH12_200_atRAL5hr_vs_RDH12_100_atRAL5hr

#treatment effect control recovery
#GFP_100_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr
#GFP_200_atRAL5hr.24h_recv_vs_GFP_control_atRAL5hr.24h_recvr
#GFP_200_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr

#treatment effect RDH12 recovery
#RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr
#RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr.24h_recvr
#RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr.24h_recvr


#recovery effect RDH12
#RDH12_100_atRAL5hr.24h_recvr_vs_RDH12_100_atRAL5hr
#RDH12_200_atRAL5hr.24h_recvr_vs_RDH12_200_atRAL5hr
#RDH12_control_atRAL5hr.24h_recvr_vs_RDH12_control_atRAL5hr

#genotype effect
#RDH12_100_atRAL5hr.24h_recvr_vs_GFP_100_atRAL5hr.24h_recvr
#RDH12_200_atRAL5hr.24h_recvr_vs_GFP_200_atRAL5hr.24h_recvr
#RDH12_control_atRAL5hr.24h_recvr_vs_GFP_control_atRAL5hr.24h_recvr


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

dep <- add_rejections(data_diff, alpha = 0.055, lfc = log2(2))

##Plot initial volcano plots
plot_volcano(dep, contrast = "RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 3, add_names = TRUE) +
  labs(title = " RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr")
plot_volcano(dep, contrast = "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr", label_size = 3, add_names = TRUE) +
  labs(title = "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr")


#plot centered protein abudance for individual proteins
#plot_single(dep, proteins = "HSPA6", type = "centered") +
 # labs(title = "HSPA6", x = "Condition")


#plot_single(dep, proteins = "HSPA5", type = "centered") +
  #labs(title = "HSPA5", x = "Condition")



#plot actual fold change for indivudal proteins
#plot_single(dep, proteins = c("HSPA5", "HSPA6"))+
 # labs(title = "Channel")
#plot_single(dep, proteins = c("RDH12", "RDH11"))+
 # labs(title = "Channel")


# Generate a results table
data_results <- get_results(dep)

# Protein-level response contrasts:
#   [(RDH12_dose - RDH12_vehicle) - (GFP_dose - GFP_vehicle)]
# GFP samples are present for the 24 h recovery condition, so these contrasts use
# the matched RDH12 recovery conditions.
add_response_contrasts <- function(data_results, se, alpha = 0.055, lfc = log2(2)) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("The limma package is required for response contrasts.")
  }

  contrast_definitions <- c(
    RDH12_100_atRAL5hr.24h_recvr_minus_vehicle_vs_GFP_100_atRAL5hr.24h_recvr_minus_vehicle =
      "(RDH12_100_atRAL5hr.24h_recvr - RDH12_control_atRAL5hr.24h_recvr) - (GFP_100_atRAL5hr.24h_recvr - GFP_control_atRAL5hr.24h_recvr)",
    RDH12_200_atRAL5hr.24h_recvr_minus_vehicle_vs_GFP_200_atRAL5hr.24h_recvr_minus_vehicle =
      "(RDH12_200_atRAL5hr.24h_recvr - RDH12_control_atRAL5hr.24h_recvr) - (GFP_200_atRAL5hr.24h_recvr - GFP_control_atRAL5hr.24h_recvr)"
  )

  expr_mat <- SummarizedExperiment::assay(se)
  sample_info <- as.data.frame(SummarizedExperiment::colData(se))

  if (!"condition" %in% colnames(sample_info)) {
    stop("Missing required colData column: condition")
  }

  required_conditions <- c(
    "RDH12_100_atRAL5hr.24h_recvr",
    "RDH12_200_atRAL5hr.24h_recvr",
    "RDH12_control_atRAL5hr.24h_recvr",
    "GFP_100_atRAL5hr.24h_recvr",
    "GFP_200_atRAL5hr.24h_recvr",
    "GFP_control_atRAL5hr.24h_recvr"
  )
  missing_conditions <- setdiff(required_conditions, unique(as.character(sample_info$condition)))
  if (length(missing_conditions) > 0) {
    stop(
      "Cannot calculate response contrasts because these conditions are missing:\n  ",
      paste(missing_conditions, collapse = "\n  ")
    )
  }

  design <- model.matrix(~ 0 + condition, data = sample_info)
  colnames(design) <- sub("^condition", "", colnames(design))

  fit <- limma::lmFit(expr_mat, design)
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_definitions,
    levels = design
  )
  colnames(contrast_matrix) <- names(contrast_definitions)
  fit_contrasts <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))

  feature_info <- as.data.frame(SummarizedExperiment::rowData(se)) %>%
    dplyr::transmute(
      name = as.character(.data$name),
      ID = as.character(.data$ID)
    )

  for (contrast_name in colnames(contrast_matrix)) {
    contrast_tbl <- limma::topTable(
      fit_contrasts,
      coef = contrast_name,
      number = Inf,
      sort.by = "none"
    ) %>%
      dplyr::transmute(
        name = feature_info$name,
        ID = feature_info$ID,
        !!paste0(contrast_name, "_p.val") := .data$P.Value,
        !!paste0(contrast_name, "_p.adj") := .data$adj.P.Val,
        !!paste0(contrast_name, "_significant") :=
          !is.na(.data$adj.P.Val) & !is.na(.data$logFC) &
          .data$adj.P.Val <= alpha & abs(.data$logFC) >= lfc,
        !!paste0(contrast_name, "_ratio") := .data$logFC
      )

    data_results <- data_results %>%
      dplyr::left_join(contrast_tbl, by = c("name", "ID"))
  }

  data_results
}

data_results <- add_response_contrasts(data_results, data_imp_mixed)


# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

###Plot centered  log2intensity heatmaps

png(  filename = file.path(path_figures, "Heatmap_log2CentInt.png"),
      width = 8, height = 8, units = "in", res = 1200)

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


