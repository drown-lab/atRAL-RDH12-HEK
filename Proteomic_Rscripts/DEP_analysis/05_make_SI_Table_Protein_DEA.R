# Derive SI Table 2 (protein differential abundance) from the DEP pipeline output.
# Reproduces the column layout of the previously hand-made SI_Table_Protein_DEA.csv
# (Gene name, Uniprot ID, <contrast>_p.val / _p.adj / _ratio, <condition>_centered, missclass_<condition>)
# with the condition labels used in the manuscript (GFP -> WT, control -> Veh), and now also carries
# the <contrast>_basis columns written by 04_WriteCsv_DEP_missingclass_BHadjusment.R.
#
# Run after 99_RunAll_Scripts_working.R; then run Proteomic_Rscripts/XICs/build_SI_support_table.R.

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })

in_file  <- "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv"
out_files <- c("Proteomic_output_txts/SI_Table_Protein_DEA.csv",
               "manuscript/SI_materials/SI_Tables/SI_Table_Protein_DEA.csv")

res <- read_csv(in_file, show_col_types = FALSE)

# dplyr:: prefixes are deliberate: after the DEP pipeline has loaded SummarizedExperiment, S4Vectors masks
# dplyr::rename/first and BiocGenerics masks several base functions, which breaks unqualified calls.
si <- res |>
  dplyr::select(-dplyr::any_of(c("significant", "Gene")), -dplyr::ends_with("_significant")) |>
  dplyr::rename(`Gene name` = name, `Uniprot ID` = ID)

# manuscript condition labels
names(si) <- names(si) |> str_replace_all("GFP", "WT") |> str_replace_all("control", "Veh")

# column order: ids, p.val, p.adj, ratio, basis, centered, missclass
ord <- c("Gene name", "Uniprot ID",
         str_subset(names(si), "_p\\.val$"), str_subset(names(si), "_p\\.adj$"),
         str_subset(names(si), "_ratio$"),   str_subset(names(si), "_basis$"),
         str_subset(names(si), "_centered$"), str_subset(names(si), "^missclass_"))
stopifnot(setequal(ord, names(si)))
si <- si[, ord]

# write_csv() does not create directories, and the manuscript/ target is gitignored, so on any
# machine without it the run used to abort here with the first file already written.
for (f in out_files) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  write_csv(si, f, na = "NA")
}
message("Wrote SI Table 2: ", nrow(si), " protein groups x ", ncol(si), " columns")

# summary of what the both-imputed rule changed, at the thresholds used in the Results text
for (ct in str_remove(str_subset(names(si), "_p\\.adj$"), "_p\\.adj$")) {
  sig <- !is.na(si[[paste0(ct, "_p.adj")]]) & si[[paste0(ct, "_p.adj")]] < 0.01 & abs(si[[paste0(ct, "_ratio")]]) >= 1
  cat(sprintf("%-70s up %4d  down %4d  (both-imputed removed: %d)\n", ct,
              sum(sig & si[[paste0(ct, "_ratio")]] > 0), sum(sig & si[[paste0(ct, "_ratio")]] < 0),
              sum(si[[paste0(ct, "_basis")]] == "both_imputed", na.rm = TRUE)))
}
