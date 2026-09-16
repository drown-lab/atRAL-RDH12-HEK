---
editor_options:
  markdown:
    wrap: 72
---

# Manuscript Outputs And Script Provenance

This is a documentation map for manuscript-facing analysis. It was
updated after the output folders were reorganized into manuscript-ready
and not-for-manuscript locations.

Status labels:

-   `main-text` - inferred main-text figure/table candidate.
-   `supplemental` - supplemental figure/table or manuscript support.
-   `QC` - quality-control output useful for validation but not core
    manuscript evidence.
-   `exploratory` - useful working analysis that needs review before
    citation.
-   `deprecated` - old, copied, stale, moved, or superseded code/output.

## Canonical Analysis Pipelines

| Analysis area | Status | Primary scripts | Primary outputs | Notes |
|----|----|----|----|----|
| Proteomics DEP | `main-text` | `Proteomic_Rscripts/DEP_analysis/00_Setup_Input_Files_and_Output_paths.R`; `Proteomic_Rscripts/DEP_analysis/99_RunAll_Scripts_working.R` | `Proteomic_output_txts/`; `Proteomic_Figs/`; `RData/Hekcells_Prep.rds` | Canonical proteomics pipeline. The setup script currently points final DEP output at `Proteomic_output_txts/DEPresults_test.csv`; keep that behavior until the file name is verified. |
| Proteomics figure scripts | `main-text` / `supplemental` | `Proteomic_Rscripts/Figures/` | `Proteomic_Figs/` | Figure scripts are separate from the DEP run-all workflow. Manuscript-ready pathway heatmaps are now grouped under `Proteomic_Figs/enriched/1_Maintext_version/`. |
| Lipidomics recovery DEA | `main-text` / `supplemental` | `Lipidomics/RScripts/1_Initial_Processing_MRMdata_script.R`; `Lipidomics/RScripts/2_Setup_for_DEA.R`; `Lipidomics/RScripts/3_PerformDEA.R` | `Lipidomics/output_txts/`; `Lipidomics/Figures/Recovery/` | Recovery lipidomics workflow and figure exports. |
| Lipidomics acute DEA | `main-text` / `supplemental` | `Lipidomics/RScripts/Acute/1_Initial_Processing_MRMdata_script.R`; `Lipidomics/RScripts/Acute/2_SetupForDEA.R`; `Lipidomics/RScripts/Acute/3_Perform_DEA.R` | `Lipidomics/output_txts/`; `Lipidomics/Figures/Acute/` | Acute lipidomics workflow and figure exports. |
| GO enrichment | `supplemental` | `GOanalysis/scripts/Batch_BioProcessOutput.R`; `GOanalysis/scripts/Fix_RDH12_200_recovery_down_BP_dotplot.R` | `GOanalysis/output/clusterProfiler_batch/` | Current manuscript support should use batch folders. The fix script shortens one term label and redraws the `RDH12_200_..._recvr_vs_RDH12_control_..._recvr_down` simplified dotplot used in Fig S7C; rerun it after the batch script. |
| Reactome pathway analysis | `supplemental` | `GOanalysis/scripts/Batch_PathwayOutput.R`; Reactome output folders under `PathwayAnalysis/` | `PathwayAnalysis/output/Reactome_batch/`; `PathwayAnalysis/output/pathway_batch_summary.csv` | Current Reactome manuscript support should use batch outputs. |
| PANGEA exports | `supplemental` | `Proteomic_Rscripts/Create_PANGEA_tables.R` | `Proteomic_output_txts/PANGEA_exports/` | Gene lists read by the GO and Reactome batch scripts. |
| PPF dataset comparison | `exploratory` / `supplemental` | `PPF_datasets/01_Prep_data.R` through `PPF_datasets/06_RDH_GFP_abundance_correlation_scatter.R` | `PPF_datasets/output_missingness/`; `PPF_datasets/output_upset/`; `PPF_datasets/output_venn/`; `PPF_datasets/output_ranked_abundance/`; `PPF_datasets/output_correlation/` | Basal HEK293T WT/RDH12 PPF comparison outputs. Mark as supplemental or exploratory until the manuscript need is confirmed. |

## Main-Text Figure Candidates

These are inferred from script names, output names, and the cleaned
manuscript-ready output folders. Review against the current manuscript
figure list before final submission.

| Figure/output candidate | Status | Script | Output location | Purpose |
|----|----|----|----|----|
| Proteomics Figure 1 pathway heatmap | `main-text` | `Proteomic_Rscripts/Figures/HeatMap_Fig1_horizontal_maintext.R` | `Proteomic_Figs/enriched/1_Maintext_version/Fig1_horizontal_maintext_heatmap_centered_sigOnly.pdf`; `Proteomic_Figs/enriched/1_Maintext_version/Fig1_horizontal_maintext_heatmap_centered_sigOnly.svg`; `Proteomic_Figs/enriched/1_Maintext_version/Fig1_horizontal_maintext_heatmap_input_table.csv` | Main-text pathway heatmap candidate. |
| Proteomics Figure 4A pathway heatmap | `main-text` | `Proteomic_Rscripts/Figures/HeatMap_Fig4wFerrop.R` | `Proteomic_Figs/enriched/1_Maintext_version/Fig4a_pathway_heatmap_centered_noSigAnno_WTfarRight.pdf`; `Proteomic_Figs/enriched/1_Maintext_version/Fig4a_pathway_heatmap_centered_noSigAnno_WTfarRight.svg`; `Proteomic_Figs/enriched/1_Maintext_version/Fig4a_pathway_heatmap_centered_noSigAnno_WTfarRight_input.csv` | Main-text pathway heatmap candidate, with ferroptosis-oriented variant. |
| Ferroptosis main-text heatmap | `main-text` | `Proteomic_Rscripts/Figures/HeatMap_Ferroptosis_maintext_Fig.R` | `Proteomic_Figs/ferroptosis_heatmap/1Ferroptosis_process_heatmap_RDH12_acute_recovery.pdf`; `Proteomic_Figs/ferroptosis_heatmap/1Ferroptosis_process_heatmap_RDH12_acute_recovery.svg`; `Proteomic_Figs/ferroptosis_heatmap/Ferroptosis_process_heatmap_input_table.csv` | Ferroptosis-focused main-text candidate. |
| RDH12 abundance rank plots | `main-text` / `supplemental` | `Proteomic_Rscripts/Figures/RDH12_abundance_rank_plots.R` | `Proteomic_Figs/RDH12_abundance_rank_plots/` | Rank-abundance plots and supporting tables. |
| Ferroptosis rank-abundance highlights | `main-text` / `supplemental` | `Proteomic_Rscripts/Figures/Rank_abundance_ferroptosis_highlights.R` | `Proteomic_Figs/rank_abundance_ferroptosis_highlights/` | Ferroptosis highlight rank plots and tables. |
| Treatment response ratio scatter | `main-text` / `supplemental` | `Proteomic_Rscripts/Figures/Scatter_treatment_response_ratios_by_genotype.R` | `Proteomic_Figs/Treatment_response_ratio_scatters/` | Genotype/treatment response comparison plot. |
| 100 uM versus 200 uM log2FC scatter | `main-text` / `supplemental` | `Proteomic_Rscripts/Figures/Scatter_100uM_vs_200uM_log2FC.R` | `Proteomic_Figs/RDH12_vs_GFP_response_contrasts/` | Dose-response comparison plot and table. |
| Individual protein abundance barplots | `supplemental` | `Proteomic_Rscripts/POI_Log2Int_acrossConditions_barplot.R` | `Proteomic_Figs/BarPlots of invidiual proteins/` | Individual protein barplots grouped out of the root figure folder. |
| Lipidomics acute class/count figures | `main-text` / `supplemental` | `Lipidomics/RScripts/Acute/Acute_Count_SigUpDownperclass.R`; `Lipidomics/RScripts/Acute/Acute_DistinctPrecursors_LipidClass_Figure.R`; `Lipidomics/RScripts/Acute/Acute_DBE_Composition_Figure.R` | `Lipidomics/Figures/Acute/` | Acute lipid class, precursor, and DBE composition figures. |
| Lipidomics recovery class/count figures | `main-text` / `supplemental` | `Lipidomics/RScripts/Recovery_Count_SigUpDownperclass.R`; `Lipidomics/RScripts/Recovery_DistinctPrecursors_LipidClass_Figure.R`; `Lipidomics/RScripts/Recovery_DBE_Composition_Figure.R` | `Lipidomics/Figures/Recovery/` | Recovery lipid class, precursor, and DBE composition figures. |

## Supplemental Tables And Outputs

| Output | Status | Script/source | Location | Notes |
|----|----|----|----|----|
| Protein DEA SI table | `supplemental` | `Proteomic_Rscripts/DEP_analysis/04_WriteCsv_DEP_missingclass_BHadjusment.R` | `Proteomic_output_txts/SI_Table_Protein_DEA.csv` | Supplemental protein differential-expression table. |
| Protein abundance SI table | `supplemental` | Proteomics DEP workflow | `Proteomic_output_txts/SI_Table_Protein Abundance.csv` | Supplemental protein abundance table. |
| Lipid acute DEA and abundance SI tables | `supplemental` | `Lipidomics/RScripts/Acute/3_Perform_DEA.R`; downstream lipid scripts | `Lipidomics/output_txts/SI_Lipid_Acute_DEA.csv`; `Lipidomics/output_txts/SI_Lipid_Acute_Abundance.csv` | Supplemental acute lipidomics tables. |
| Lipid recovery DEA and abundance SI tables | `supplemental` | `Lipidomics/RScripts/3_PerformDEA.R`; downstream lipid scripts | `Lipidomics/output_txts/SI_Lipid_Recovery_DEA.csv`; `Lipidomics/output_txts/SI_Lipid_Recovery_abundance.csv` | Supplemental recovery lipidomics tables. |
| Proteomics volcano plots | `supplemental` | `Proteomic_Rscripts/Create_Volcanos_allContrasts_toDirctry.R`; `Proteomic_Rscripts/Create_Volcanos_wMissingess_allContrasts_toDirctry.R` | `Proteomic_Figs/VolcanoPlots/StandardVolcano2/`; `Proteomic_Figs/VolcanoPlots/wMissingessInfo/` | Current volcano outputs are organized into standard and missingness-aware folders. |
| Lipidomics volcano plots | `supplemental` | `Lipidomics/RScripts/Acute/Acute_Volcanos_stnd.R`; `Lipidomics/RScripts/Recovery_Volcanos_stnd.R` | `Lipidomics/Figures/Acute/VolcanoPlots/`; `Lipidomics/Figures/Recovery/VolcanoPlots/` | Acute and recovery lipid volcano plots. |
| Proteomics PCA and missingness QC | `QC` | `Proteomic_Rscripts/DEP_analysis/01a_QC_HeatMap_Missingness.R`; `01b_QC_PCA_RDH12_acuteAtRAL.R`; `01b_QC_PCA_Recovery.R`; `02-c_MissingnessAbundance.R` | `Proteomic_Figs/Proteomic_PCAs/`; `Proteomic_Figs/PCA_AcuteRDH12.png`; `Proteomic_Figs/HeatMap_missingness.png`; `Proteomic_Figs/Missingness_Abund_Distribution.png` | PCA exports were moved from top-level `Proteomic_PCAs/` into `Proteomic_Figs/Proteomic_PCAs/`. |
| Lipidomics CV QC (Fig S6A/B) | `supplemental` | `Lipidomics/RScripts/QC_CV_lipids_rawdata.R` | `Lipidomics/Figures/QC_CV/` | Lipid CV% from raw intensities for both datasets. |
| Lipidomics PCA (Fig S6C/D) | `supplemental` | `Lipidomics/RScripts/Acute/2a_PCA_normVSN.R`; `Lipidomics/RScripts/2a_QC_PCA_data_norm.R` | `Lipidomics/Figures/Acute/QC_PCA_normVSN.pdf`; `Lipidomics/Figures/Recovery/QC_pca_normvsn_SE.pdf` | VSN-normalized PCAs. The PDFs were saved from the plot pane; the scripts only write `PCA_Recovery.png` into `path_figures`. |
| GO and Reactome dotplots | `supplemental` | `GOanalysis/scripts/`; `Proteomic_Rscripts/Figures/Supplemental_Acute_RDH12_GO_Reactome_Figure.R` | `GOanalysis/output/clusterProfiler_batch/`; `PathwayAnalysis/output/Reactome_batch/`; `Proteomic_Figs/Supplemental_Acute_RDH12_GO_Reactome/` | Use batch/current outputs for manuscript support. Not-for-manuscript outputs are documented in `ARCHIVE_NOTES.md`. |

## Not-For-Manuscript Or Historical Outputs

| Output location | Status | Notes |
|----|----|----|
| Pruned 2026-09-16 | `deprecated` | Superseded scripts, not-for-manuscript enrichment outputs, PANGEA downloads, and stale exports were deleted. See `ARCHIVE_NOTES.md`. |
| `Proteomic_PCAs/` | `deprecated` | Top-level PCA folder was emptied/removed; current PCA exports are in `Proteomic_Figs/Proteomic_PCAs/`. |
| `PathwayAnalysis/output/Reactome/` | `deprecated` | Non-batch Reactome output folder was removed; use `PathwayAnalysis/output/Reactome_batch/`. |
| `Lipidomics/RScripts/04_pretty volcano-lightschannel.R` | `deprecated` | Removed from the working tree. It previously contained hardcoded `Cell_Perturbed/Paper1/...` paths. |

## Review Flags

-   `Proteomic_output_txts/DEPresults_test.csv` is the current path
    configured in
    `Proteomic_Rscripts/DEP_analysis/00_Setup_Input_Files_and_Output_paths.R`.
    Treat it as a live configured output, not as junk, until the final
    file name is verified.
-   Several generated output folders are currently deleted or moved in
    the worktree. The manifest follows the new organized locations where
    present and marks old locations as historical.
