# Characterizing Response to atRAL Toxicity

This is a working repo for proteomic and lipidomic data that was generated in the atRAL toxicity experiments in the Rams collab. Experimental approach: Cell model expressing RDH12. Dose-response was performed by treatment with atRAL to induce toxicity. Cells were collected at 5hrs post-treatment, or 24hrs post-treatment with a media switch.

## Organization

This workflow is for differential expression analysis of proteomic and lipidomic datasets.

#### Description of Folders

-   Folder `Proteomic_Rscripts`- contains scripts used to generate DEP analysis

-   Folder `Proteomic_Figs` - contains figures related to proteomic data

-   Folder `Proteomic_output_txts` - contains key output txt-like files

-   Folder `PCAs` - contains visualizations of proteomic data

-   Folder `quant_DIANN_outputs`- contains the filtered excel worksheet of protein IDs and intensity

-   Folder `Lipidomics`- contains scripts, figures, output_txts from lipidomics datasets

#### Differential Expression Analysis:

-   Differential expression of proteins workflow (DEP): is found in `Proteomic_Rscripts/DEP_analysis`: comes with own Readme file for detailed explanation
-   Differential Expression of lipids: TBA

#### More Lipid scripts To be added

## Updates

-   As of 3/13/26: this repo only contains R scripts for proteomics and lipidomics analysis
    -   Currently assessing if can perform lipidomics DEA on the combined datasets-as there could be batch effect from both instrument collection date and a experiment-type

## Version, Dependencies, Packages

Information on what is required/used to run scripts can be found in `dependciesAndPackages_info.txt`
