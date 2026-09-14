
#Setup Input and Output paths/files 
#in the works

#IMPORTANT: Do NOT change the table names that are being created for each file path
#The scripts require those table names for the Run ALL to work
#After Setup this up run `00_RunAll_Scripts`
#====================================
#Setup input path/files
#=====================================

# 1) CSV file of Raw Data from DIANN; already filtered file by Q.Value <= 0.01,PG.Q.Value <= 0.05, Lib.Q.Value <= 0.01,
#Lib.PG.Q.Value <= 0.01, Channel.Q.Value <= 0.05, Proteotypic == 1,!grepl("cRAP", Protein.Ids, ignore.case = TRUE)
#Find csv file
#DO not change "Protmcs_HekatRALexps_v1" or scripts will not work
Protmcs_HekatRALexps_v1 <-  read_csv("quant_DIANN_outputs/HekRDH12andGFP_combinedDatasets_atRAL5hr_with24hrRecvry_rawdata.csv")

 
exists("Protmcs_HekatRALexps_v1")
colnames(Protmcs_HekatRALexps_v1)


#Metadata produced in script no files needed


###=======================================
# Setup Output paths and file names
#========================================
#This is for Where do you want output files to go
#IMPORTANT: Do NOT change the table names that are being created for each file path
#The scripts require those table names for the Run ALL to work

dir.create("RData", showWarnings = FALSE)

# 1) Setup output path where want Figures to go 
path_figures <- "Proteomic_Figs/"

# 2) Setup Output to save CSV file for final DEP results with missingness class identification
path_DEP_results<- "Proteomic_output_txts/DEPresults_test.csv"

