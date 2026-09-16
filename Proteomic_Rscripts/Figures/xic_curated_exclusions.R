# =========================================================
# XIC-curated protein exclusions and single-peptide labelling for the
# main-text heatmaps (Fig 3H, 5E, 6A) and the supporting heatmap scripts.
#
# Source this file from a figure script (working directory = repo root):
#   source("Proteomic_Rscripts/Figures/xic_curated_exclusions.R")
#
# The --xic 90 re-search (plot_xics_single_peptide.R;
# manuscript/SI_materials/XICs/XIC_curated_named_proteins.pdf) showed several named proteins
# resting on noise-level chromatograms or entirely imputed values. Verdicts are recorded in
# manuscript/SI_materials/Revision_plan_single_peptide_XICs.docx; the REMOVE set is listed
# here once so every heatmap applies the same curation without editing its pathway gene lists.
# An explicit list, not a Support_tier filter: GSTA4, SLC31A1 and CHKA are tier A on peptide
# count but noise-level in the XICs.
#
# Provides
#   xic_excluded_genes     named character vector: gene symbol -> reason
#   apply_xic_curation()   drops excluded genes from a data frame and reports what was removed
#   single_peptide_ids()   Uniprot IDs of single-peptide protein groups, read from
#                          Proteomic_output_txts/SI_Table_Protein_Support.csv (tracked), else the
#                          manuscript/ copy (built by XICs/build_SI_support_table.R)
#   label_single_peptide() appends "*" to single-peptide protein group labels
#                          (figure legend: "* single-peptide protein group; see SI XIC data")
# =========================================================

xic_excluded_genes <- c(
  NFE2L2  = "single peptide, Evidence 2.3, no coeluting fragments in any run; recovery increase is an imputed-floor artefact",
  SLC39A7 = "single 26-residue 3+ peptide, Evidence 3.1, no coelution",
  SLC39A1 = "single peptide, 5 of 27 runs, Evidence 2.9, noise",
  SLC30A5 = "single peptide, Evidence 3.1, MS1 profile correlation 0 (no MS1 support)",
  SLC31A1 = "two peptides but <= 5 runs each, Evidence 4.0 / 3.0, sparse and weak",
  GSTA4   = "two peptides, both noise-level (Evidence 2.9 / 2.5); recovery 'surge' is imputed floor",
  CHKA    = "two peptides, Evidence 2.0 / 2.7, single dominant fragment only; ETNK1 (3 peptides) carries the Kennedy-pathway claim",
  BCL2L12 = "single peptide, Evidence 2.6, noise",
  PIDD1   = "single peptide, Evidence 2.0, noise",
  CYBA    = "MNAR_0of3 in all three acute conditions; acute depletion is a QRILC imputed-floor draw, which also inflates the recovery cells by ~1.2 log2 through the per-protein centering"
)

# Tracked copy first, manuscript working copy second. manuscript/ is gitignored, so only the
# tracked copy makes the "*" markers reproducible on a fresh clone.
support_table_candidates <- c(
  "Proteomic_output_txts/SI_Table_Protein_Support.csv",
  "manuscript/SI_materials/SI_Tables/SI_Table_Protein_Support.csv"
)
support_table_file <- {
  found <- support_table_candidates[file.exists(support_table_candidates)]
  if (length(found)) found[1] else support_table_candidates[1]
}

# Remove excluded genes from a data frame; `gene_col` is the column holding gene symbols.
apply_xic_curation <- function(df, gene_col = "Gene", label = "") {
  hit <- as.character(df[[gene_col]]) %in% names(xic_excluded_genes)
  if (any(hit)) {
    message(sprintf("XIC curation%s: removed %d protein(s): %s",
                    if (nzchar(label)) paste0(" [", label, "]") else "",
                    sum(hit), paste(sort(unique(as.character(df[[gene_col]][hit]))), collapse = ", ")))
  } else {
    message(sprintf("XIC curation%s: nothing to remove", if (nzchar(label)) paste0(" [", label, "]") else ""))
  }
  df[!hit, , drop = FALSE]
}

# Uniprot IDs (Protein.Group strings, as in the DEP table `ID` column) of single-peptide groups.
single_peptide_ids <- function(path = support_table_file) {
  if (!file.exists(path)) {
    stop("Support table not found. Looked for:\n  ",
         paste(support_table_candidates, collapse = "\n  "),
         "\nRun Proteomic_Rscripts/XICs/build_SI_support_table.R first. Without it the ",
         "figures would ship with no single-peptide '*' markers but a legend claiming them.")
  }
  sup <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(c("Protein.Group", "Single_peptide_protein_group") %in% names(sup)))
  sup$Protein.Group[sup$Single_peptide_protein_group %in% c(TRUE, "TRUE")]
}

# Display label with "*" appended for single-peptide protein groups.
label_single_peptide <- function(gene_label, uniprot_id, sp_ids = single_peptide_ids(), mark = "*") {
  is_sp <- !is.na(uniprot_id) & uniprot_id %in% sp_ids
  ifelse(is_sp, paste0(gene_label, mark), gene_label)
}
