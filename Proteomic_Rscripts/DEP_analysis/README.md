---
editor_options: 
  markdown: 
    wrap: 72
---

# atRAL Toxicity Proteomics - Data Analysis

Data analysis workflow for differential expression analysis of proteins
(DEA or DEP)

## How to utilize the scripts to generate DEP output:

-   Step 1: `00_Setup_Input_Files_and_Output_paths` - requires the user
    to set the input file paths and output file paths/names. It is
    important to NOT change the table names that those file paths are
    called to or the scripts will not work
-   Step 2: `99_RunAll_Scripts_working` - This runs all the scripts
    necessary to perform analysis for proteomics data; This will create
    Figures and output_txts

#### Understand what the scripts do that Run All is performing:

-   Step 1: `01_Datasetup_for_DEPanalysis.R`- initial loading of
    filtered CSV file with raw data output

    -   QC Step 1a: `02-01a_QC_HeatMap_Missingness.R` - heatmap
        visualization of the presence/absence of missingness in the
        dataset

    -   QC Step 1b: `02-01b_QC_PCA_RDH12_acuteAtRAL.R` - creates PCA
        plot of acute atRAL data

    -   QC Step 1c: `02-QC Step 1b:`02-01b_QC_PCA_RDH12_acuteAtRAL.R\` -
        creates PCA plot of Recovery atRAL data

-   Step 2: `02-a_DEP_imputationMixedmodel.R` - performs mixed
    imputation on missing values. MAR- knn and MNAR- QRILC

    -   QC Step 2a: `02-b_Missingnesssummarization.R` - provides tabular
        results of the protein missingness in the dataset into
        missingness classes MNAR = 0/3, MAR= 1/3 or 2/3, or Present 3/3

    -   Step 2b: `02-c_MissingnessAbundance.R` - provides visualizations
        of the protein abundance for missingness in the dataset

-   Step 3: `03_DEP_performDEA` - perform the differential expression
    analysis of proteins (DEA or DEP)

-   Step 4:`04_Data_results_with_missing_class_output` - creates table
    proving the status of missing values for each protein in each
    condition

## Dependencies

Find dependencies in `output_txts` folder -\> `sessionInfo` and
`packageversions` Analysis workflows depend on quantification data from
DIANN These data files are handled by [Git
LFS](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage),
so be sure to have it installed to fetch these files.
