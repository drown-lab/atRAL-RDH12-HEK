# Derive SI Table 2 (protein differential abundance) from the DEP pipeline output.
# Reproduces the column layout of the previously hand-made SI_Table_Protein_DEA.csv
# (Gene name, Uniprot ID, <contrast>_p.val / _p.adj / _ratio, <condition>_centered, missclass_<condition>)
# with the condition labels used in the manuscript (GFP -> WT, control -> Veh), and now also carries
# the <contrast>_basis columns written by 04_WriteCsv_DEP_missingclass_BHadjusment.R.
#
# Runs as the last step of 99_RunAll_Scripts_working.R (it is in that script's `scripts` vector);
# afterwards run Proteomic_Rscripts/XICs/build_SI_support_table.R.

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

# Thresholds for the counts quoted in the Results text. These are deliberately stricter than the
# alpha_cutoff = 0.05 / lfc_cutoff = log2(1.5) that 04_WriteCsv_DEP_missingclass_BHadjusment.R uses
# for its `<contrast>_significant` columns. Those columns are read by no downstream script and are
# dropped from the SI table above, so the Results counts are defined here and only here -- but the
# stale flags do still ship inside data_results_w_missingclass_BH_readjusted.csv, where they will
# not reproduce these numbers. Keep the two in view of each other when either is changed.
RESULTS_ALPHA          <- 0.01
RESULTS_MIN_ABS_LOG2FC <- 1        # 2-fold

# summary of what the both-imputed rule changed, at the thresholds used in the Results text
cat(sprintf("\nSignificance for the counts below: p.adj < %g and |log2 FC| >= %g (%.3g-fold)\n",
            RESULTS_ALPHA, RESULTS_MIN_ABS_LOG2FC, 2^RESULTS_MIN_ABS_LOG2FC))
for (ct in str_remove(str_subset(names(si), "_p\\.adj$"), "_p\\.adj$")) {
  sig <- !is.na(si[[paste0(ct, "_p.adj")]]) & si[[paste0(ct, "_p.adj")]] < RESULTS_ALPHA &
    abs(si[[paste0(ct, "_ratio")]]) >= RESULTS_MIN_ABS_LOG2FC
  cat(sprintf("%-70s up %4d  down %4d  (both-imputed removed: %d)\n", ct,
              sum(sig & si[[paste0(ct, "_ratio")]] > 0), sum(sig & si[[paste0(ct, "_ratio")]] < 0),
              sum(si[[paste0(ct, "_basis")]] == "both_imputed", na.rm = TRUE)))
}
