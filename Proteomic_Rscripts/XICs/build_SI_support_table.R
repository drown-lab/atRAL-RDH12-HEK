# Build per-protein-group evidence/support columns from the deposited DIA-NN report and
# append them to the protein DEA supplementary table.
#
# Outputs (in si_dir):
#   SI_Table_Protein_Support.csv           - one row per protein group (8,116 after the V12 filter chain)
#   SI_Table_Protein_DEA_with_support.csv  - SI Table 2 + support columns + per-contrast quantification basis
#
# Support columns
#   Distinct_proteotypic_peptides  distinct stripped sequences (proteotypic precursors only, q-value filters as in V12)
#   Distinct_precursors            distinct precursor ids (sequence + modifications + charge)
#   Library_proteotypic_peptides   distinct proteotypic peptides for the group in the refined (MBR) library
#   Runs_detected_of_27            runs in which >= 1 proteotypic precursor passed the filters
#   Detected_<condition>_of_3      per-condition replicate detection counts (9 conditions)
#   Median_DIANN_Evidence          median of DIA-NN "Evidence" (fragment-ion coelution score) over identified precursor-runs
#   Median_MS1_profile_corr        median DIA-NN Ms1.Profile.Corr over identified precursor-runs
#   Single_peptide_protein_group   TRUE if Distinct_proteotypic_peptides == 1
#   Support_tier                   A: >= 2 proteotypic peptides
#                                  B: single peptide, median Evidence >= 4  (chromatograms show coeluting fragments)
#                                  C: single peptide, median Evidence <  4  (weak / noise-level; not used as individual evidence)
#   <contrast>_basis               "measured" (both conditions have >= 1 observed value), "num_imputed" / "den_imputed"
#                                  (one side is MNAR_0of3, i.e. entirely QRILC-imputed), "both_imputed"
#
# The Evidence threshold of 4 was chosen after inspecting XICs for 17 named proteins: every chromatogram with
# convincing fragment coelution scored >= 4.3 and every noise-level one <= 3.1 (see Revision_plan_single_peptide_XICs.docx).

suppressPackageStartupMessages({ library(arrow); library(dplyr); library(tidyr); library(stringr); library(readr) })
# Bioconductor packages loaded by the DEP pipeline (matrixStats, S4Vectors, IRanges, BiocGenerics) mask several
# dplyr verbs (count, first, rename, desc, slice, ...). The masked ones used here are called with dplyr:: explicitly,
# so the script behaves the same in a fresh session and after 99_RunAll_Scripts_working.R.

report_orig <- "F:/MS_Temp/atRAL_manuscript/proteomics/report.parquet"
report_lib  <- "F:/MS_Temp/atRAL_manuscript/proteomics/report_lib.parquet"
si_dir      <- "C:/Users/bryon/atRAL-Stress-Rams-Collab/manuscript/SI_materials/SI_Tables"
dea_file    <- file.path(si_dir, "SI_Table_Protein_DEA.csv")
evidence_threshold <- 4

# ---- 1. precursor-level rows, filtered exactly as in 00_Process_DIANNparquetfile.R ----------------------
cols <- c("Run", "Precursor.Id", "Stripped.Sequence", "Proteotypic", "Protein.Group", "Protein.Ids", "Genes",
          "Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "Channel.Q.Value", "Evidence", "Ms1.Profile.Corr")
rep <- read_parquet(report_orig, col_select = all_of(cols)) |>
  filter(Q.Value <= 0.01, PG.Q.Value <= 0.05, Lib.Q.Value <= 0.01, Lib.PG.Q.Value <= 0.01, Channel.Q.Value <= 0.05,
         Proteotypic == 1, !grepl("cRAP", Protein.Ids, ignore.case = TRUE)) |>
  mutate(Run = basename(Run) |> str_remove("\\.raw$"))

# ---- 2. sample annotation from run names --------------------------------------------------------------
sample_from_run <- function(run) {
  acute <- str_match(run, "Hek_(Control|100uMatRAL|200uMatRAL)_Rep([ABC])")
  recov <- str_match(run, "Hek(GFP|RDH12)_(EtOH|100uM|200uM)_Rep(\\d)")
  genotype  <- ifelse(!is.na(acute[,1]), "RDH12", ifelse(recov[,2] == "GFP", "WT", "RDH12"))
  treatment <- ifelse(!is.na(acute[,1]),
                      c(Control = "Veh", `100uMatRAL` = "100", `200uMatRAL` = "200")[acute[,2]],
                      c(EtOH = "Veh", `100uM` = "100", `200uM` = "200")[recov[,3]])
  phase <- ifelse(!is.na(acute[,1]), "atRAL5hr", "atRAL5hr.24h_recvr")
  paste(genotype, treatment, phase, sep = "_")            # matches the DEA / missclass condition names
}
rep <- rep |> mutate(condition = sample_from_run(Run))
stopifnot(n_distinct(rep$Run) == 27, !any(is.na(rep$condition)))
conditions <- c("RDH12_Veh_atRAL5hr", "RDH12_100_atRAL5hr", "RDH12_200_atRAL5hr",
                "RDH12_Veh_atRAL5hr.24h_recvr", "RDH12_100_atRAL5hr.24h_recvr", "RDH12_200_atRAL5hr.24h_recvr",
                "WT_Veh_atRAL5hr.24h_recvr", "WT_100_atRAL5hr.24h_recvr", "WT_200_atRAL5hr.24h_recvr")
stopifnot(setequal(unique(rep$condition), conditions))

# ---- 3. per-protein-group summaries --------------------------------------------------------------------
pg_summary <- rep |>
  group_by(Protein.Group) |>
  summarise(Gene_names = dplyr::first(Genes),
            Distinct_proteotypic_peptides = n_distinct(Stripped.Sequence),
            Distinct_precursors = n_distinct(Precursor.Id),
            Runs_detected_of_27 = n_distinct(Run),
            Median_DIANN_Evidence = median(Evidence, na.rm = TRUE),
            Median_MS1_profile_corr = median(Ms1.Profile.Corr, na.rm = TRUE),
            .groups = "drop")

per_cond <- rep |>
  dplyr::distinct(Protein.Group, condition, Run) |>
  dplyr::count(Protein.Group, condition, name = "n") |>
  mutate(condition = factor(condition, levels = conditions)) |>
  complete(Protein.Group, condition, fill = list(n = 0L)) |>
  pivot_wider(names_from = condition, values_from = n, names_glue = "Detected_{condition}_of_3")

lib <- read_parquet(report_lib, col_select = c("Precursor.Id", "Stripped.Sequence", "Proteotypic", "Decoy", "Protein.Group")) |>
  filter(Decoy == 0, Proteotypic == 1) |>
  dplyr::distinct(Protein.Group, Stripped.Sequence) |>
  dplyr::count(Protein.Group, name = "Library_proteotypic_peptides")

support <- pg_summary |>
  left_join(lib, by = "Protein.Group") |>
  left_join(per_cond, by = "Protein.Group") |>
  mutate(Single_peptide_protein_group = Distinct_proteotypic_peptides == 1,
         Support_tier = case_when(
           Distinct_proteotypic_peptides >= 2 ~ "A",
           Median_DIANN_Evidence >= evidence_threshold ~ "B",
           TRUE ~ "C")) |>
  relocate(Library_proteotypic_peptides, .after = Distinct_precursors) |>
  arrange(Protein.Group)

write_csv(support, file.path(si_dir, "SI_Table_Protein_Support.csv"))
message(sprintf("Support table: %d protein groups; single-peptide %d (tier B %d, tier C %d)",
                nrow(support), sum(support$Single_peptide_protein_group),
                sum(support$Support_tier == "B"), sum(support$Support_tier == "C")))

# ---- 4. append to the DEA table with a per-contrast quantification-basis column ------------------------
dea <- read_csv(dea_file, show_col_types = FALSE)
stopifnot("Uniprot ID" %in% names(dea))

contrasts <- names(dea) |> str_subset("_p\\.adj$") |> str_remove("_p\\.adj$")
basis_of <- function(num_class, den_class) {
  ni <- num_class == "MNAR_0of3"; di <- den_class == "MNAR_0of3"
  case_when(is.na(num_class) | is.na(den_class) ~ NA_character_,
            ni & di ~ "both_imputed", ni ~ "num_imputed", di ~ "den_imputed", TRUE ~ "measured")
}
for (ct in contrasts) {
  if (paste0(ct, "_basis") %in% names(dea)) next          # already provided by 04_WriteCsv_DEP_missingclass_BHadjusment.R
  parts <- str_split_fixed(ct, "_vs_", 2)
  num_col <- paste0("missclass_", parts[1]); den_col <- paste0("missclass_", parts[2])
  if (!all(c(num_col, den_col) %in% names(dea))) { warning("missing missclass columns for ", ct); next }
  dea[[paste0(ct, "_basis")]] <- basis_of(dea[[num_col]], dea[[den_col]])
}

dea_out <- dea |>
  left_join(support |> dplyr::select(-Gene_names), by = c("Uniprot ID" = "Protein.Group")) |>
  dplyr::relocate(Distinct_proteotypic_peptides, Distinct_precursors, Library_proteotypic_peptides, Runs_detected_of_27,
           Median_DIANN_Evidence, Median_MS1_profile_corr, Single_peptide_protein_group, Support_tier,
           .after = `Uniprot ID`)
stopifnot(!any(is.na(dea_out$Distinct_proteotypic_peptides)))   # every DEA row must be found in the report

write_csv(dea_out, file.path(si_dir, "SI_Table_Protein_DEA_with_support.csv"))
message("Wrote ", file.path(si_dir, "SI_Table_Protein_DEA_with_support.csv"))

# quick sanity report for the proteins discussed in the revision plan
named <- c("MGST2","SLC11A2","TUBB2A","TUBB4B","ALG3","AGPAT3","SLC7A11","CYBA","CHKA",
           "NFE2L2","SLC39A7","BCL2L12","PIDD1","GSTA4","SLC39A1","SLC30A5","SLC31A1","HMOX1","SRXN1","RDH12")
print(support |> filter(Gene_names %in% named) |>
        select(Gene_names, Distinct_proteotypic_peptides, Runs_detected_of_27, Median_DIANN_Evidence, Support_tier) |>
        arrange(Support_tier, desc(Median_DIANN_Evidence)), n = 30)
