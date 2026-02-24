# Create folder if needed
dir.create("Proteomic_output_txts",
           recursive = TRUE, showWarnings = FALSE)

# 1. Package table
ip <- installed.packages()
write.table(
  data.frame(Package = ip[, "Package"],
             Version = ip[, "Version"]),
  file = "Proteomic_output_txts/package_versions.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 2. Full session info
sink("Proteomic_output_txts/sessionInfo.txt")
print(sessionInfo())
sink()
