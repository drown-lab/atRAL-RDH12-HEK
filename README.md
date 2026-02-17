# Characterizing Response to atRAL Toxicity

This is a working repo for proteomic and lipidomic data that was generated in the atRAL toxicity experiments in the Rams collab. Experimental approach: Cell model expressing RDH12. Dose-response was performed by treatment with atRAL to induce toxicity. Cells were collected at 5hrs post-treatment, or 24hrs post-treatment with a media switch.

## Organization

This workflow is for differential expression analysis of proteomic and lipidomic datasets.

Folder `Proteomic_Rscripts`- contains scripts used to generate DEP analysis

Folder `quant_DIANN_outputs`- contains the filtered excel worksheet of protein IDs and intensity

#### How to utilize `Proteomic_Rscripts`in a step wise manner:

-   Step 1: `Datasetup_for_DEPanalysis`- takes an already filtered excel worksheet and puts data into summarized experiment table for compatibility with DEP bioconductor package; There are initial plots in this script that will look at QC visualizations of the data (currently these plots will not save to a path)

    -   QC Step 1a: `QC_HeatMap_Missingness`- plots presence/absence heatmap based on missingness per sample
    -   QC Step 1b: `QC_PCA_RDH12_actueAtRAL`- plots PCA before imputation for Acute atRAL treatment on RDH12 cells
    -   QC Step 1b: `QC_PCA_Recovery` - plots PCA before imputation for Recovery on RDH12 and Control cells
    -   QC Step 1c: `QC_Missingness_summarization` - creates tables that summarize the missingness in the dataset into missingness classes MNAR = 0/3, MAR= 1/3 or 2/3, or Present 3/3
    -   QC Step 1c: `QC_Missingness_AbundanceVisualization` - plots the log2 intensity distribution of the MAR and present 3/3

-   Step 2: mixed model imputation

#### How to utilize Lipidomic scripts in a step wise manner:

To be added

## Updates

-   As of 2/2/26: this repo only contains R scripts for proteomics analysis

## Version, Dependencies, Packages

Information on what is required/used to run scripts can be found in `dependciesAndPackages_info.txt`
