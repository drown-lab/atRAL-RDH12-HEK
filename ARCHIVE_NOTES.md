---
editor_options:
  markdown:
    wrap: 72
---

# Archive Notes For Noncanonical Code

This file documents scripts and folders that are not currently treated
as the main manuscript analysis path, and records what has been removed
from the working tree.

## Cleanup Policy

-   Delete code and outputs only after tracing the manuscript figures
    and tables back to their scripts and confirming nothing kept depends
    on them.
-   Record each pruning pass here. Deleted files stay recoverable from
    git history.
-   Treat the review tables below as a queue: items can be promoted to
    manuscript-facing status or retired in a later pass.

## Pruned 2026-09-16

Checked against manuscript draft V12 and `SI_doc_v2.docx` (branch
`prune-stale`). Recover any of these with
`git checkout <commit before the prune> -- <path>`.

| Removed | Reason |
|------------------------|------------------------|
| `test/`, `Lipidomics/test/`, `GOanalysis/test3.R` | Scratch scripts; several pointed at other projects (`Cell_HalfLife/`, `Cell_Perturbed/`). |
| `Proteomic_Rscripts/Figures/depracated/` (13 scripts) | Superseded by the current `Proteomic_Rscripts/Figures/` scripts. |
| `Proteomic_Rscripts/Acute_RDH12_Heatmap_*.R` (5), `GetTransporterProteins.R`, `Heatmap_clusterSigProts_*.R`, `not sure if used/protein heatmap.R` | Superseded by `HeatMap_Fig1_horizontal_maintext.R`, or never produced an output folder. |
| `Proteomic_Rscripts/ProteinCount_wUporDown - Copy.R`, `Volcano_wMissingness_callContrast.R`, `Remove_other_contams.R`, `Other_QCscripts/Upsetplot_checkingmissingness.R`, `Other_QCscripts/VarianceCheckofDEP.R` | Copies, console snippets, or logic now in `DEP_analysis/01_Datasetup_for_DEPanalysis.R`; their `DEPresults_v1/v2` and `Volcanotable_v1/v2` tables were read by nothing. |
| `GOanalysis/scripts/BioProcess.R`, `Pathways.R` | Single-contrast versions of the batch scripts. |
| Lipid `Rvr_Volcanos_stnd.R`, `Acute/HeatMap_KeyLipids.R`, `Rvr_countUpDown_lollipopperclass.R`, `Acute_Count_UpDown.R`, `Rvr_Count_UpDown.R`, `Acute/4_VolcanoPlot.R`, `not sure if used/lipid recovery heatmap.R` | Superseded by `Recovery_Volcanos_stnd.R`, `HeatMap_KeyLipids2.R`, `Recovery_Count_SigUpDownperclass.R`, and `Rvr_HeatMap.R`, or not in V12. |
| Lipid `5_Breakdown_*` / `5_Recovery_Breakdown_*` scripts and their summary CSVs | Replaced by the DBE composition figure scripts. |
| `Lipidomics/RScripts/CombinedDatasets/`, `Lipidomics/Figures/CombinedDatasets/`, `Sample_description guide_combined.csv` | Dropped batch-combined (ComBat) approach. |
| `Acute/POIs_barplot*.R`, `Acute/Acute_Scatter_100vs200_Veh_log2FC.R` and their figure folders | Not in V12. |
| `Acute/2a_PCA_beforeNorm.R`, `Acute/2a_PCA_normVSNImp.R`, `2a_QC_PCA_data_se.R` | Fig S6C/D come from `Acute/2a_PCA_normVSN.R` and `2a_QC_PCA_data_norm.R` (PC percentages match). |
| `Acute/1a_Filtering_Score_script.R` | Byte-identical to `Lipidomics/RScripts/1a_Filtering_Score_script.R`, which the acute pipeline now sources. |
| `GOanalysis/output/clusterProfiler_Not for manuscript/`, `PathwayAnalysis/output/Reactome_not for manuscript/` | Outputs already marked not for manuscript. |
| `Proteomic_output_txts/PANGEA_results_*`, `Proteomic_Figs/PANGEA_top_terms_plots/` | PANGEA web downloads and plots from deprecated scripts; Methods use clusterProfiler/ReactomePA. `PANGEA_exports/` is kept (batch GO/Reactome input). |
| Older heatmap versions in `Proteomic_Figs/enriched/`, stale exports (`Heatmap.png`, lipid UpSet/Alllipids PDFs, extra lipid count layouts), `RawData_long_LFQ_table.csv` | Superseded or regenerable; not used by the manuscript. |
| `dependciesAndPackages_info.txt` | Replaced by `renv.lock`. |

`not sure if used/lipid CV report.R` was kept and moved to
`Lipidomics/RScripts/QC_CV_lipids_rawdata.R` (Fig S6A/B).

## Deprecated Or Superseded

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Lipidomics/RScripts/04_pretty volcano-lightschannel.R` | `deprecated` | Removed from the working tree. It previously contained hardcoded `Cell_Perturbed/Paper1/...` paths that did not match this repo structure. |
| `00_ReadFilesIN.R` | `deprecated` / `review` | Removed from the working tree. The canonical proteomics DEP workflow now uses `Proteomic_Rscripts/DEP_analysis/00_Setup_Input_Files_and_Output_paths.R` for input setup. |

## Test Or Exploratory Folders

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Proteomic_Rscripts/Other_QCscripts/` | `QC` | Holds `QC_plots_UpsetandCV_rawdata.R` (Fig S2A-C), outside the canonical DEP run-all workflow. |

## Not For Manuscript Output Folders

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Rplots.pdf` | `deprecated` | Generic R graphics-device output currently shown as deleted in the worktree. Do not restore unless needed. |

## Organized But Not Primary Manuscript Evidence

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Proteomic_Figs/BarPlots of invidiual proteins/` | `supplemental` / `review` | Individual protein barplots were moved out of the root figure folder. Useful for provenance and possible supplement, but not currently marked as main-text. |
| `Proteomic_Figs/Proteomic_PCAs/` | `QC` | Current location for proteomics PCA exports after moving them out of the former top-level `Proteomic_PCAs/` folder. |

## Stale Or Needs Review

| Path | Status | Reason |
|------------------------|------------------------|------------------------|
| `Proteomic_output_txts/DEPresults_test.csv` | `review` | Despite the name, this is the currently configured output in `Proteomic_Rscripts/DEP_analysis/00_Setup_Input_Files_and_Output_paths.R`. Rename only after confirming downstream scripts and manuscript tables. |
| `Proteomic_PCAs/` | `deprecated` | Former top-level PCA output folder; current PCA files are in `Proteomic_Figs/Proteomic_PCAs/`. |
| `PathwayAnalysis/output/Reactome/` | `deprecated` | Former non-batch Reactome output folder; current manuscript support should use `PathwayAnalysis/output/Reactome_batch/`. |
| Not-in-V12 working analyses | `review` | Kept pending reviewer comments: missingness-aware volcanoes, `KEGG_batch/`, rank-abundance and ratio scatters, 100 vs 200 uM overlap and its GO run, `HeatMap_KEAP1-detoxROS.R`, `HeatMap_MultiPathways_paper.R`, `ProteinCount_wUporDown.R`, POI barplots, `PPF_datasets/02_missingness_heatmap.R`. |
