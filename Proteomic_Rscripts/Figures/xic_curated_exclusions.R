# =========================================================
# XIC-curated protein exclusions and single-peptide labelling for the
# main-text heatmaps (Fig 3H, 5E, 6A) and the supporting heatmap scripts.
#
# Source this file from a figure script (working directory = repo root):
#   source("Proteomic_Rscripts/Figures/xic_curated_exclusions.R")
#
# Why this exists
#   The DIA-NN re-search with --xic 90 (Proteomic_Rscripts/XICs/plot_xics_single_peptide.R;
#   manuscript/SI_materials/XICs/XIC_curated_named_proteins.pdf) showed that several proteins
#   named in the manuscript are supported only by noise-level chromatograms or by entirely
#   imputed values. The verdicts are recorded in
#   manuscript/SI_materials/Revision_plan_single_peptide_XICs.docx. Proteins with a REMOVE
#   verdict (or a SOFTEN verdict that drops them from a figure) are listed here ONCE, with
#   the reason, so every heatmap script applies the same curation and the pathway gene
#   lists in those scripts remain the original, uncurated membership lists.
#
#   Note that this is deliberately an explicit list rather than a Support_tier filter:
#   GSTA4, SLC31A1 and CHKA have two peptides (tier A by the automatic rule) but both
#   peptides are noise-level in the XICs, so a tier filter alone would not remove them.
#
# What it provides
#   xic_excluded_genes     named character vector: gene symbol -> one-line reason
#   apply_xic_curation()   drops excluded genes from a data frame and reports what was removed
#   single_peptide_ids()   Uniprot IDs (Protein.Group) of single-peptide protein groups, read from
#                          manuscript/SI_materials/SI_Tables/SI_Table_Protein_Support.csv
#                          (built by Proteomic_Rscripts/XICs/build_SI_support_table.R)
#   label_single_peptide() appends "*" to the display label of single-peptide protein groups
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
  PIDD1   = "single peptide, Evidence 2.0, noise"
)

support_table_file <- "manuscript/SI_materials/SI_Tables/SI_Table_Protein_Support.csv"

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
    warning("Support table not found (", path, "); single-peptide labels will not be applied. ",
            "Run Proteomic_Rscripts/XICs/build_SI_support_table.R first.")
    return(character(0))
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
