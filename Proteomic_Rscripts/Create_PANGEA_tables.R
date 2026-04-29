library(dplyr)
library(stringr)
library(readr)
library(purrr)

# --------------------------------------------------
# 1) Read DEP results
# --------------------------------------------------
Cells <- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")

# Use this as the working dataframe
df <- Cells

# --------------------------------------------------
# 2) Identify contrasts from *_p.val columns
# --------------------------------------------------
pval_cols <- colnames(df)[str_detect(colnames(df), "_p\\.val$")]
contrasts <- str_remove(pval_cols, "_p\\.val$")

# --------------------------------------------------
# 3) Output directory
# --------------------------------------------------
outdir <- "Proteomic_output_txts/PANGEA_exports"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# 4) Thresholds
# --------------------------------------------------
padj_cutoff <- 0.01
lfc_cutoff  <- log2(2)

# --------------------------------------------------
# 5) Function to make tables for one contrast
#    - all proteins with stats
#    - up proteins
#    - down proteins
# --------------------------------------------------
make_pangea_tables <- function(df, ct,
                               padj_cutoff = 0.05,
                               lfc_cutoff = log2(1.5)) {
  
  pval_col  <- paste0(ct, "_p.val")
  padj_col  <- paste0(ct, "_p.adj")
  ratio_col <- paste0(ct, "_ratio")
  
  needed_cols <- c("name", "ID", "Gene", pval_col, padj_col, ratio_col)
  missing_cols <- setdiff(needed_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns for contrast ", ct, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out_df <- df %>%
    transmute(
      name  = name,
      Gene  = Gene,
      ID    = ID,
      p.val = as.numeric(.data[[pval_col]]),
      p.adj = as.numeric(.data[[padj_col]]),
      ratio = as.numeric(.data[[ratio_col]])
    ) %>%
    filter(!is.na(Gene), Gene != "") %>%
    distinct(Gene, .keep_all = TRUE)
  
  list(
    all = out_df,
    up = out_df %>%
      filter(!is.na(ratio), !is.na(p.adj),
             ratio >=  lfc_cutoff,
             p.adj <= padj_cutoff),
    down = out_df %>%
      filter(!is.na(ratio), !is.na(p.adj),
             ratio <= -lfc_cutoff,
             p.adj <= padj_cutoff)
  )
}

# --------------------------------------------------
# 6) Build named list of tables
# --------------------------------------------------
pangea_tables <- setNames(
  lapply(contrasts, function(ct) {
    make_pangea_tables(
      df = df,
      ct = ct,
      padj_cutoff = padj_cutoff,
      lfc_cutoff = lfc_cutoff
    )
  }),
  contrasts
)

# --------------------------------------------------
# 7) Write csv files
# --------------------------------------------------
for (ct in names(pangea_tables)) {
  for (type in names(pangea_tables[[ct]])) {
    
    out_df <- pangea_tables[[ct]][[type]]
    filename <- paste0(ct, "_", type, ".csv")
    
    write.csv(
      out_df,
      file = file.path(outdir, filename),
      row.names = FALSE
    )
  }
}

# --------------------------------------------------
# 8) Optional summary table
# --------------------------------------------------
summary_tbl <- map_dfr(names(pangea_tables), function(ct) {
  tibble(
    contrast = ct,
    n_all  = nrow(pangea_tables[[ct]]$all),
    n_up   = nrow(pangea_tables[[ct]]$up),
    n_down = nrow(pangea_tables[[ct]]$down)
  )
})

write.csv(
  summary_tbl,
  file.path(outdir, "PANGEA_export_summary.csv"),
  row.names = FALSE
)

print(summary_tbl)