##read parquet
##load libraries
library(arrow)
library(reshape2)
library(readxl)
library(tidyverse)
library(readr)
library(dplyr)
#######################

input_dir <- if (dir.exists("PPF_datasets")) "PPF_datasets" else "."
report_path <- file.path(input_dir, "report.parquet")
metadata_path <- file.path(input_dir, "metadata.txt")
filtered_path <- file.path(input_dir, "Protein_filtered.csv")

# report.parquet is not in git; it is deposited on PRIDE (PXD080689). Without it,
# load PPFtable_v4 from Protein_filtered.csv, the tracked output of the steps below.
if (!file.exists(report_path)) {
  if (!file.exists(filtered_path)) {
    stop("Could not find ", report_path, " or ", filtered_path,
         ". Download the timsTOF basal-set report.parquet from PRIDE PXD080689 into ", input_dir, ".")
  }
  message("report.parquet not found; loading PPFtable_v4 from ", filtered_path,
          ". To rebuild it, download report.parquet from PRIDE PXD080689 into ", input_dir, ".")
  # na = "NA" keeps the empty Genes strings as "" rather than NA, matching the parquet build
  PPFtable_v4 <- read_csv(
    filtered_path,
    col_types = cols(.default = col_character(), PG.MaxLFQ = col_double()),
    na = "NA",
    name_repair = "unique_quiet"
  ) |>
    select(-1)
} else {

if (!file.exists(metadata_path)) {
  stop("Could not find metadata file: ", metadata_path)
}

sample_metadata <- read_tsv(
  metadata_path,
  col_names = c("run_number", "sample_name"),
  col_types = cols(
    run_number = col_character(),
    sample_name = col_character()
  )
) |>
  distinct(run_number, .keep_all = TRUE)

PPFtable_v1 <- read_parquet(report_path)
exists("PPFtable_v1")
colnames(PPFtable_v1)

#filter PG based on certain values
##filter for only proteotypic and get rid of "crap" or decoy protein IDs

PPFtable_v2 <- PPFtable_v1 |>
  filter(Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & Channel.Q.Value <= 0.05)|>
  filter(Proteotypic == 1)|>
  filter(!grepl("cRAP", Protein.Ids, ignore.case = TRUE))

#known human contaminants:
contaminant_proteins <- c(
  "Q6E0U4", "P20930", "Q5D862", "Q86YZ3", "P04264", "P13645", "Q99456",
  "P13646", "P02533", "P19012", "P08779", "Q04695", "P05783", "P08727",
  "P35908", "P35900", "Q8N1A0", "Q9C075", "Q2M2I5", "Q7Z3Z0", "Q7Z3Y9",
  "Q7Z3Y8", "Q7Z3Y7", "P12035", "Q15323", "Q14532", "O76009", "Q14525",
  "O76011", "Q92764", "O76013", "O76014", "O76015", "Q6A163", "P19013",
  "Q6A162", "P13647", "P02538", "P04259", "P48668", "Q3KNV1", "P08729",
  "Q3SY84", "Q14CN4", "Q86Y46", "Q7RTS7", "O95678", "Q01546", "Q7Z794",
  "Q8N1N4", "Q5XKE5", "P05787", "Q6KB66", "Q14533", "Q9NSB4", "P78385",
  "Q9NSB2", "P78386", "O43790", "A6NCN2", "P35527", "P60331", "P60014",
  "P60412", "P60413", "P60368", "P60369", "P60372", "P60370", "P60371",
  "P60409", "P60410", "P60411", "Q07627", "Q8IUC1", "P59990", "P59991",
  "P60328", "P60329", "Q8IUG1", "Q8IUC0", "Q52LG2", "Q3SY46", "Q3LI77",
  "P0C5Y4", "Q9BYS1", "Q3LI76", "A8MUX0", "Q9BYP8", "Q8IUB9", "Q3LHN2",
  "Q7Z4W3", "Q3LI73", "Q3LI72", "Q3LI70", "Q3SYF9", "Q3LI54", "Q3LI63",
  "Q3LI61", "Q9BYU5", "Q3LI58", "Q3LI59", "Q3LHN1", "Q9BYT5", "Q3MIV0",
  "P0C7H8", "A1A580", "Q9BYR9", "Q3LI83", "Q3LHN0", "Q6PEX3", "Q3LI81",
  "A8MX34", "Q9BYR8", "Q9BYR7", "Q9BYR6", "Q9BYQ7", "Q9BYQ6", "Q9BQ66",
  "Q9BYR5", "Q9BYR4", "Q9BYR3", "Q9BYR2", "Q9BYQ5", "Q9BYR0", "Q9BYQ9",
  "Q9BYQ8", "Q6L8H4", "Q6L8G5", "Q6L8G4", "Q701N4", "Q6L8H2", "Q6L8H1",
  "Q701N2", "Q6L8G9", "Q6L8G8", "O75690", "P26371", "Q3LI64", "Q3LI66",
  "Q3LI67", "Q8IUC3", "Q8IUC2", "A8MXZ3", "Q9BYQ4", "Q9BYQ3", "Q9BYQ2",
  "A8MVA2", "A8MTY7", "Q9BYQ0", "Q9BYP9"
)

PPFtable_v3 <- PPFtable_v2 |>
  filter(!Protein.Group %in% contaminant_proteins)

PPFtable_v4 <- PPFtable_v3 |>
  mutate(run_number = str_extract(Run, "(?<=Vijay-)\\d+(?=_)")) |>
  left_join(sample_metadata, by = "run_number") |>
  mutate(
    sample_run_name = if_else(
      is.na(sample_name),
      Run,
      paste(sample_name, Run, sep = "_")
    )
  )

unmapped_runs <- PPFtable_v4 |>
  filter(is.na(sample_name)) |>
  distinct(Run, run_number)

if (nrow(unmapped_runs) > 0) {
  warning(
    "Some runs did not map to metadata: ",
    paste(unmapped_runs$Run, collapse = ", ")
  )
}

PPFtable_v4 <- PPFtable_v4 |>
  select(sample_name, sample_run_name, Run, run_number, Protein.Group, Genes, PG.MaxLFQ) |>
  unique()

write.csv(PPFtable_v4, filtered_path)

} # end rebuild from report.parquet
