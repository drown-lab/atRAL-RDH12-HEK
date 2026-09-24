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
-   `GOanalysis/scripts/` and `PathwayAnalysis/` - GO and Reactome enrichment analysis scripts and output folders. Manuscript outputs are in `GOanalysis/output/clusterProfiler_batch/` and `PathwayAnalysis/output/Reactome_batch/`.
-   `PPF_datasets/` - PPF-specific proteomics comparison scripts and outputs.

## Organization

This workflow is for differential expression analysis of proteomic and lipidomic datasets.

### Folder Guide

-   `quant_DIANN_outputs/` - DIA-NN quantification output used as proteomics input.
-   `Proteomic_Rscripts/` - proteomics analysis and figure scripts.
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

## Raw Data Locations

Run scripts with the repo root as the working directory (open `atRALExps.Rproj`); paths inside the repo are relative to it. Raw inputs too large for git are read from `raw_data/` (gitignored), or from the folder named by an environment variable:

| Input | Environment variable | Default | Used by |
|----|----|----|----|
| DIA-NN `report.parquet`, `report_lib.parquet`, `DIANN_xic_rerun/` | `ATRAL_REPORT_DIR` | `raw_data/diann` | `Proteomic_Rscripts/00_Process_DIANNparquetfile.R`; `Proteomic_Rscripts/XICs/` |
| Acute lipid MRM exports (`Absolute_intensity.csv` folders) | `ATRAL_LIPID_ACUTE_DIR` | `raw_data/lipids/acute` | `Lipidomics/RScripts/Acute/1_Initial_Processing_MRMdata_script.R` |
| Recovery lipid MRM exports | `ATRAL_LIPID_RECOVERY_DIR` | `raw_data/lipids/recovery` | `Lipidomics/RScripts/1_Initial_Processing_MRMdata_script.R` |

You only need the variables for the inputs you rerun, and none if you copy the data into the default folders above.

### Setting the variables

1.  Open your user-level `.Renviron` from the R console. On Windows it is usually `C:/Users/<you>/Documents/.Renviron`; `path.expand("~/.Renviron")` prints the exact location.

    ``` r
    file.edit(path.expand("~/.Renviron"))
    ```

2.  Add one line per variable, pointing at the folder on this machine (the paths below are examples). Use forward slashes; spaces in paths are fine without quotes:

    ```
    ATRAL_REPORT_DIR=F:/MS_Temp/atRAL_manuscript/proteomics
    ATRAL_LIPID_ACUTE_DIR=D:/lipidomics/20250702_HekCells_atRALtreated_rams/Life Sciences Native LIpids
    ATRAL_LIPID_RECOVERY_DIR=D:/lipidomics/Hek_5hrw24hrrecvr_Miranda18samples
    ```

3.  Save, restart R (RStudio: *Session > Restart R*), and check:

    ``` r
    Sys.getenv("ATRAL_REPORT_DIR")
    ```

To set a variable for the current session only, run `Sys.setenv(ATRAL_REPORT_DIR = "F:/MS_Temp/atRAL_manuscript/proteomics")` before sourcing the script.

Put these in `~/.Renviron`, not in a `.Renviron` inside the project: R reads only one `.Renviron` at startup, so a project-level file would hide your user-level settings.

## Noncanonical And Exploratory Code

Scripts and outputs not used by the manuscript were pruned in September 2026 and remain recoverable from git history. See [ARCHIVE_NOTES.md](ARCHIVE_NOTES.md) for what was removed and for the working analyses kept for review.

Known noncanonical locations include:

-   `Proteomic_Rscripts/Other_QCscripts/`

## Version, Dependencies, Packages

Package versions are pinned in `renv.lock` (R 4.4.1, Bioconductor 3.20); `renv_setup.R` restores the environment. Additional R session information is stored in output folders such as `Proteomic_output_txts/sessionInfo.txt`.
