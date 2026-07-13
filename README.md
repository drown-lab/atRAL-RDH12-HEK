# Characterizing Response to atRAL Toxicity

This repository contains proteomic and lipidomic analysis code for atRAL toxicity experiments from the Rams-Drown-Ferreira collaboration. The cell model is HEK293T cells that are either WT or overexpress RDH12. Cells were treated with atRAL and collected either 5 hours after treatment or after a 24 hour recovery period with a media switch.

The repo is a working analysis project, but the manuscript-facing pieces are now documented separately from exploratory, deprecated, and quality-control scripts.

## Manuscript Analysis Map

Use [MANUSCRIPT_OUTPUTS.md](MANUSCRIPT_OUTPUTS.md) as the detailed map of scripts, figure/table outputs, and manuscript status.

Core manuscript-facing analysis paths:

-   `Proteomic_Rscripts/DEP_analysis/` - canonical proteomics DEP workflow. Start with `00_Setup_Input_Files_and_Output_paths.R`, then run `99_RunAll_Scripts_working.R`.
-   `Proteomic_Rscripts/Figures/` - proteomics figure-generating scripts, including inferred main-text and supplemental figure candidates.
-   `Proteomic_Figs/` - proteomics figures and figure input tables. Manuscript-ready pathway heatmaps are grouped in `Proteomic_Figs/enriched/1_Maintext_version/`.
-   `Proteomic_output_txts/` - proteomics tables, SI tables, DEA outputs, PANGEA exports, and enrichment-ready outputs.
-   `Proteomic_Figs/Proteomic_PCAs/` - proteomics PCA figure exports.
-   `Lipidomics/RScripts/` - recovery lipidomics processing, DEA, summaries, and figure scripts.
-   `Lipidomics/RScripts/Acute/` - acute lipidomics processing, DEA, summaries, and figure scripts.
-   `Lipidomics/Figures/` - lipidomics figure outputs.
-   `Lipidomics/output_txts/` - lipidomics DEA, abundance, and summary tables.
-   `GOanalysis/scripts/` and `PathwayAnalysis/` - GO and Reactome enrichment analysis scripts and output folders. Outputs moved into `GOanalysis/output/clusterProfiler_Not for manuscript/` and `PathwayAnalysis/output/Reactome_not for manuscript/` are retained for provenance but are not manuscript-facing.
-   `PPF_datasets/` - PPF-specific proteomics comparison scripts and outputs.

## Organization

This workflow is for differential expression analysis of proteomic and lipidomic datasets.

### Folder Guide

-   `quant_DIANN_outputs/` - DIA-NN quantification output used as proteomics input.
-   `Proteomic_Rscripts/` - proteomics analysis and figure scripts.
-   `Proteomic_Rscripts/Figures/depracated/` - older proteomics figure scripts retained for provenance but no longer treated as active manuscript figure drivers.
-   `Proteomic_Rscripts/DEP_analysis/` - canonical proteomics DEP pipeline.
-   `Proteomic_Figs/` - proteomics figures, QC plots, heatmaps, volcano plots, and figure source tables.
-   `Proteomic_output_txts/` - proteomics tables and text-like outputs.
-   `Proteomic_Figs/Proteomic_PCAs/` - proteomics PCA exports.
-   `Lipidomics/` - lipidomics scripts, figures, metadata, and output tables.
-   `GOanalysis/` - GO enrichment scripts and clusterProfiler outputs.
-   `PathwayAnalysis/` - Reactome pathway outputs.
-   `PPF_datasets/` - PPF dataset- basal condition of HEK293T WT and RDH12 cells- preprocessing and comparison outputs.
-   `RData/` - serialized intermediate R objects.
-   `Environments/` - dependency/environment helper script.

## Differential Expression Analysis

-   Proteomics DEP workflow: `Proteomic_Rscripts/DEP_analysis/README.md`
-   Proteomics run-all script: `Proteomic_Rscripts/DEP_analysis/99_RunAll_Scripts_working.R`
-   Lipidomics recovery workflow: `Lipidomics/RScripts/`
-   Lipidomics acute workflow: `Lipidomics/RScripts/Acute/`

## Noncanonical And Exploratory Code

The cleanup policy for this repository is conservative: preserve old scripts and document their status instead of deleting them. See [ARCHIVE_NOTES.md](ARCHIVE_NOTES.md) for the current list of deprecated, test, exploratory, or stale files.

Known noncanonical locations include:

-   `test/`
-   `Lipidomics/test/`
-   `not sure if used/`
-   `Proteomic_Rscripts/Depracated/`
-   `Proteomic_Rscripts/Figures/depracated/`
-   `Proteomic_Rscripts/Other_QCscripts/`

## Version, Dependencies, Packages

Package and dependency notes are in `dependciesAndPackages_info.txt`. Additional R session information is stored in output folders such as `Proteomic_output_txts/sessionInfo.txt`.
