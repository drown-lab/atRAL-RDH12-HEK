---
editor_options:
  markdown:
    wrap: 72
---

# Archive Notes For Noncanonical Code

This file documents scripts and folders that are not currently treated
as the main manuscript analysis path. It also records user-created
organization changes that separate manuscript-ready outputs from
historical or not-for-manuscript outputs.

## Cleanup Policy

-   Preserve old code and generated outputs unless deletion is
    explicitly approved.
-   Prefer documenting status over deleting files, because old scripts
    may still explain how a figure or table was created.
-   Treat this file as a review queue: items can be promoted to
    manuscript-facing status or retired in a later cleanup pass.

## Deprecated Or Superseded

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Proteomic_Rscripts/Depracated/` | `deprecated` | Older proteomics DEP scripts appear superseded by `Proteomic_Rscripts/DEP_analysis/`. Folder name is misspelled in the repo; do not rename without a separate dependency check. |
| `Proteomic_Rscripts/Figures/HeatMap_Fig4wFerrop_unknownversion.R` | `deprecated` | Filename marks it as an unknown version. Prefer `HeatMap_Fig4wFerrop.R` unless manuscript review says otherwise. |
| `Proteomic_Rscripts/ProteinCount_wUporDown - Copy.R` | `deprecated` | Copy-named script. Prefer the non-copy version if this analysis is still needed. |
| `Lipidomics/RScripts/04_pretty volcano-lightschannel.R` | `deprecated` | Removed from the working tree. It previously contained hardcoded `Cell_Perturbed/Paper1/...` paths that did not match this repo structure. |
| `00_ReadFilesIN.R` | `deprecated` / `review` | Removed from the working tree. The canonical proteomics DEP workflow now uses `Proteomic_Rscripts/DEP_analysis/00_Setup_Input_Files_and_Output_paths.R` for input setup. |

## Test Or Exploratory Folders

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `test/` | `exploratory` | Contains scratch/test scripts, including an older PANGEA table script. Prefer `Proteomic_Rscripts/Create_PANGEA_tables.R` for manuscript documentation. |
| `Lipidomics/test/` | `exploratory` | Contains lipidomics figure experiment script. Not currently mapped to main analysis. |
| `not sure if used/` | `exploratory` | Folder name indicates uncertain provenance. Review before citing or deleting. |
| `Proteomic_Rscripts/Other_QCscripts/` | `QC` | Useful QC scripts, but outside the canonical DEP run-all workflow. |
| `GOanalysis/test3.R` | `exploratory` | Test-like filename outside the main `GOanalysis/scripts/` folder. |

## Not For Manuscript Output Folders

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `GOanalysis/output/clusterProfiler_Not for manuscript/` | `deprecated` / `exploratory` | User-organized location for enrichment outputs that should not be used as current manuscript outputs. |
| `PathwayAnalysis/output/Reactome_not for manuscript/` | `deprecated` / `exploratory` | User-organized location for Reactome outputs that should not be used as current manuscript outputs. |
|  |  |  |
| `Rplots.pdf` | `deprecated` | Generic R graphics-device output currently shown as deleted in the worktree. Do not restore unless needed. |

## Organized But Not Primary Manuscript Evidence

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Proteomic_Figs/BarPlots of invidiual proteins/` | `supplemental` / `review` | Individual protein barplots were moved out of the root figure folder. Useful for provenance and possible supplement, but not currently marked as main-text. |
| `Lipidomics/Figures/Acute/Barplots_individual-lipids/` | `supplemental` / `review` | Individual acute lipid barplots and boxplots were moved out of the acute figure root. |
| `Proteomic_Figs/Proteomic_PCAs/` | `QC` | Current location for proteomics PCA exports after moving them out of the former top-level `Proteomic_PCAs/` folder. |

## Stale Or Needs Review

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Proteomic_output_txts/DEPresults_test.csv` | `review` | Despite the name, this is the currently configured output in `Proteomic_Rscripts/DEP_analysis/00_Setup_Input_Files_and_Output_paths.R`. Rename only after confirming downstream scripts and manuscript tables. |
| `Proteomic_PCAs/` | `deprecated` | Former top-level PCA output folder; current PCA files are in `Proteomic_Figs/Proteomic_PCAs/`. |
| `PathwayAnalysis/output/Reactome/` | `deprecated` | Former non-batch Reactome output folder; current manuscript support should use `PathwayAnalysis/output/Reactome_batch/`. |

## 
