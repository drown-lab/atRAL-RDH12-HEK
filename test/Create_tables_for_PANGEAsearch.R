library(dplyr)
library(stringr)
library(readr)
library(purrr)

## Get txt files for PANGEA searches

# Call DEP output text file
Paper2_SILAC <- read.csv("Cell_Perturbed/Paper2/output_txts/DEPresults_test2.csv")

# inspect channel values first
unique(Paper2_SILAC$channel)

df <- Paper2_SILAC

#-----------------------------
# 1) Identify contrasts
#-----------------------------
pval_cols <- colnames(df)[str_detect(colnames(df), "_p\\.val$")]
contrasts <- str_remove(pval_cols, "_p\\.val$")

#-----------------------------
# 2) Output directory
#-----------------------------
outdir <- "Cell_Perturbed/Paper2/output_txts/PANGEA_exports"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------
# 3) Function to make tables
#-----------------------------
make_pangea_tables <- function(df, ct, channel_value) {
  pval_col  <- paste0(ct, "_p.val")
  padj_col  <- paste0(ct, "_p.adj")
  ratio_col <- paste0(ct, "_ratio")
  
  out_df <- df %>%
    filter(channel == channel_value) %>%
    transmute(
      name  = name,
      ID    = ID,
      p.val = .data[[pval_col]],
      p.adj = round(.data[[padj_col]], 2),
      ratio = round(.data[[ratio_col]], 1)
    )
  
  list(
    all  = out_df,
    up   = out_df %>% filter(ratio >= 1.5, p.adj <= 0.05),
    down = out_df %>% filter(ratio <= -1.5, p.adj <= 0.05)
  )
}

#-----------------------------
# 4) Build nested list
#-----------------------------
pangea_tables <- setNames(
  lapply(contrasts, function(ct) {
    list(
      Heavy = make_pangea_tables(df, ct, "H"),
      Light = make_pangea_tables(df, ct, "L")
    )
  }),
  contrasts
)

#-----------------------------
# 5) Write files
#-----------------------------
for (ct in names(pangea_tables)) {
  for (ch in names(pangea_tables[[ct]])) {
    for (type in names(pangea_tables[[ct]][[ch]])) {
      
      out_df <- pangea_tables[[ct]][[ch]][[type]]
      
      filename <- paste0(ct, "_", ch, "_", type, ".csv")
      
      write.csv(
        out_df,
        file = file.path(outdir, filename),
        row.names = FALSE
      )
    }
  }
}