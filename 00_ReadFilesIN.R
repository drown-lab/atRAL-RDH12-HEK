Protmcs_HekatRALexps_v1<-read_excel("Z:/data/Projects/Rams_Collab_RDH12experiments/atRAL_experiments/proteomics/quant_data_filtered_DIANNandExcels/HekRDH12andGFP_combinedDatasets_atRAL5hr_with24hrRecvry_rawdata.xlsx")
colnames(Protmcs_HekatRALexps_v1)

###########
#export data

# Define the path where you want to save the Excel file
path <- "~/GitHub/atRALExps/output_txt_DEP_MM/DEPallcontrasts_mixedimputation_v1.csv"

# Define the table you want to export
table <- DEPallcontrasts_mixedimputation_v1      
#pg_volcano_avgFCs
# Export the table to Excel
#write_xlsx(table, path)

write_csv(table,path)