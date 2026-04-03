# =========================================================
# RUN ALL – Paper2 (no setwd needed)
# =========================================================

#rm(list = ls())
#cat("\014")

#Step 0) Need to load necessary libraries before running
#Step 1) Make sure you setup 00_Setup_Input_Files_and_Output_paths.R before running
#Step 2) Perform this script

# Where the scripts actually live relative to repo root
base_dir <- file.path(getwd(), "Proteomic_Rscripts","DEP_analysis")
stopifnot(dir.exists(base_dir))

scripts <- c(
  "00_Setup_Input_Files_and_Output_paths.R",
  "01_Datasetup_for_DEPanalysis.R",
  #"01a_QC_HeatMap_Missingness.R",
  #"01b_QC_PCA_RDH12_acuteAtRAL.R",
  #"01b_QC_PCA_Recovery.R",
  "02-a_DEP_imputationMixedmodel.R",
  "02-b_Missingnesssummarization.R",
  "02-c_MissingnessAbundance.R",
  "03_DEP_performDEA.R",
  "04_Data_results_with_missing_class_output.R",
  "04_WriteCsv_DEP_missingclass_BHadjusment.R"
)

start_time <- Sys.time()

for (s in scripts) {
  
  full_path <- file.path(base_dir, s)
  
  message("\n==============================")
  message("Running: ", full_path)
  message("==============================\n")
  
  if (!file.exists(full_path)) {
    stop(
      "Script not found:\n  ", full_path,
      "\n\nFiles in base_dir are:\n  ",
      paste(list.files(base_dir), collapse = "\n  ")
    )
  }
  
  tryCatch(
    source(full_path, local = FALSE),
    error = function(e) stop("Error in script: ", s, "\n", e$message)
  )
}

end_time <- Sys.time()

message("\nALL SCRIPTS COMPLETED")
message("Total runtime: ",
        round(difftime(end_time, start_time, units = "mins"), 2),
        " minutes")
