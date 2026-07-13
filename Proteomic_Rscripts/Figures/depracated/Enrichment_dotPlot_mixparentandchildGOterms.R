library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(viridis)

# ---------------------------------------------------------
# 1. INPUT
# ---------------------------------------------------------

enrich_df <- read.csv(
  "Proteomic_output_txts/PANGEA_results/Acute_100vsVeh_100uMenrichment_2026-04-06 135232.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

outdir <- "Proteomic_Figs/PANGEA_top_terms_plots/Acute_100vsVeh_100uMenrichment/dotplots/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bh_cutoff <- 0.25

# number of terms to plot after redundancy reduction
n_pathway_terms <- 20
n_bp_terms <- 20
n_cc_terms <- 10

# redundancy settings
remove_near_duplicates <- TRUE
jaccard_cutoff <- 0.9   # terms with highly overlaping enriched genes are treated as redundant

colnames(enrich_df)

# ---------------------------------------------------------
# 2. CLEAN COLUMN TYPES
# ---------------------------------------------------------

enrich_df <- enrich_df %>%
  mutate(
    `Fold Enrichment` = as.numeric(`log2 fold`),
    `Benjamini & Hochberg` = as.numeric(`Benjamini & Hochberg`)
  )

# ---------------------------------------------------------
# 3. FIND OVERLAP GENE COLUMN
# ---------------------------------------------------------

candidate_overlap_cols <- c(
  "Overlap Gene Symbols",
  "Overlaping Gene Symbols",
  "Shared Gene Symbols",
  "Genes in Overlap",
  "Gene Symbols",
  "Input gene(s) in this gene set",
  "Genes",
  "Overlap Genes",
  "overlap_genes"
)

overlap_col <- candidate_overlap_cols[candidate_overlap_cols %in% colnames(enrich_df)][1]

if (is.na(overlap_col) || length(overlap_col) == 0) {
  stop(
    paste0(
      "Could not find the overlap gene column.\n",
      "Available columns are:\n",
      paste(colnames(enrich_df), collapse = "\n")
    )
  )
}

message("Using overlap gene column: ", overlap_col)

# ---------------------------------------------------------
# 4. HELPER FUNCTIONS
# ---------------------------------------------------------

normalize_term_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9 ]+", " ") %>%
    str_squish()
}

parse_gene_set <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  
  genes <- unlist(str_split(x, "\\s*[,;/|]\\s*"))
  genes <- str_trim(genes)
  genes <- genes[genes != ""]
  genes <- unique(genes)
  genes
}

canonical_gene_set_string <- function(x) {
  genes <- parse_gene_set(x)
  if (length(genes) == 0) return(NA_character_)
  paste(sort(genes), collapse = ";")
}

jaccard_similarity <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  inter / union
}

reduce_redundant_terms <- function(df,
                                   n_terms,
                                   bh_cutoff,
                                   overlap_col,
                                   remove_near_duplicates = TRUE,
                                   jaccard_cutoff = 0.85) {
  
  if (nrow(df) == 0) {
    return(df[0, , drop = FALSE])
  }
  
  # prepare a working table
  work_df <- df %>%
    filter(
      !is.na(`Fold Enrichment`),
      !is.na(`Benjamini & Hochberg`),
      `Benjamini & Hochberg` <= bh_cutoff
    ) %>%
    mutate(
      term_norm = normalize_term_name(`Gene Set Name`),
      overlap_string = vapply(.data[[overlap_col]], canonical_gene_set_string, character(1)),
      overlap_count = vapply(.data[[overlap_col]], function(x) length(parse_gene_set(x)), integer(1))
    ) %>%
    arrange(`Benjamini & Hochberg`, desc(`Fold Enrichment`), desc(overlap_count)) %>%
    distinct(`Gene Set Name`, .keep_all = TRUE)
  
  if (nrow(work_df) == 0) {
    return(work_df)
  }
  
  # 1) remove exact duplicates by overlap gene set
  exact_dedup_df <- work_df %>%
    mutate(overlap_string = if_else(is.na(overlap_string), paste0("NO_OVERLAP_", row_number()), overlap_string)) %>%
    distinct(overlap_string, .keep_all = TRUE)
  
  if (!remove_near_duplicates) {
    return(exact_dedup_df %>% slice_head(n = n_terms))
  }
  
  # 2) greedy removal of near-duplicate terms based on Jaccard similarity
  kept_rows <- vector("list", length = 0)
  kept_gene_sets <- list()
  
  for (i in seq_len(nrow(exact_dedup_df))) {
    current_row <- exact_dedup_df[i, , drop = FALSE]
    current_genes <- parse_gene_set(current_row[[overlap_col]])
    
    is_redundant <- FALSE
    
    if (length(kept_gene_sets) > 0) {
      sims <- vapply(kept_gene_sets, function(gs) jaccard_similarity(current_genes, gs), numeric(1))
      if (any(sims >= jaccard_cutoff, na.rm = TRUE)) {
        is_redundant <- TRUE
      }
    }
    
    if (!is_redundant) {
      kept_rows[[length(kept_rows) + 1]] <- current_row
      kept_gene_sets[[length(kept_gene_sets) + 1]] <- current_genes
    }
  }
  
  out_df <- bind_rows(kept_rows)
  
  out_df %>%
    slice_head(n = n_terms)
}

make_enrichment_plot <- function(df,
                                 n_terms,
                                 title_text,
                                 out_file,
                                 bh_cutoff,
                                 overlap_col,
                                 remove_near_duplicates = TRUE,
                                 jaccard_cutoff = 0.85) {
  
  if (nrow(df) == 0) {
    message("No rows found for: ", title_text)
    return(NULL)
  }
  
  plot_df <- reduce_redundant_terms(
    df = df,
    n_terms = n_terms,
    bh_cutoff = bh_cutoff,
    overlap_col = overlap_col,
    remove_near_duplicates = remove_near_duplicates,
    jaccard_cutoff = jaccard_cutoff
  )
  
  if (nrow(plot_df) == 0) {
    message("No significant nonredundant rows after BH filtering for: ", title_text)
    return(NULL)
  }
  
  plot_df <- plot_df %>%
    mutate(
      `Gene Set Name` = fct_reorder(`Gene Set Name`, `Fold Enrichment`)
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = `Fold Enrichment`,
      y = fct_reorder(`Gene Set Name`, `Fold Enrichment`),
      size = `Count Overlap Gene`,
      color = -log10(`Benjamini & Hochberg`)
    )
  ) +
    geom_point() +
    scale_color_viridis_c(option = "cividis") +
    labs(
      title = title_text,
      x = "Fold Enrichment",
      y = NULL,
      color = "-log10(BH)",
      size = "Gene Count"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12)
    )
  
  print(p)
  
  ggsave(
    filename = file.path(outdir, out_file),
    plot = p,
    width = 10,
    height = 7
  )
  
  return(plot_df)
}

# ---------------------------------------------------------
# 5. DEFINE CATEGORIES
# ---------------------------------------------------------

pathway_categories <- c(
  "REACTOME pathway",
  "PANTHER Pathway H.sap",
  "KEGG Pathway H.sap"
)

bp_categories <- c(
  "GO Biological Processes",
  "Direct GO Biological Processes",
  "SLIM2 GO BP"
)

cc_categories <- c(
  "GO Cellular Component",
  "Direct GO Cellular Component",
  "SLIM2 GO CC"
)

# ---------------------------------------------------------
# 6. SUBSET TABLES
# ---------------------------------------------------------

pathway_df <- enrich_df %>%
  filter(`Gene Set Category` %in% pathway_categories)

bp_df <- enrich_df %>%
  filter(`Gene Set Category` %in% bp_categories)

cc_df <- enrich_df %>%
  filter(`Gene Set Category` %in% cc_categories)

# ---------------------------------------------------------
# 7. MAKE PLOTS
# ---------------------------------------------------------

top_pathways <- make_enrichment_plot(
  df = pathway_df,
  n_terms = n_pathway_terms,
  title_text = paste0(
    "Top Nonredundant Enriched Pathways (BH <= ", bh_cutoff, ")"
  ),
  out_file = "top_nonredundant_pathways.pdf",
  bh_cutoff = bh_cutoff,
  overlap_col = overlap_col,
  remove_near_duplicates = remove_near_duplicates,
  jaccard_cutoff = jaccard_cutoff
)

top_bp <- make_enrichment_plot(
  df = bp_df,
  n_terms = n_bp_terms,
  title_text = paste0(
    "Top Nonredundant Enriched Biological Processes (BH <= ", bh_cutoff, ")"
  ),
  out_file = "top_nonredundant_biological_processes.pdf",
  bh_cutoff = bh_cutoff,
  overlap_col = overlap_col,
  remove_near_duplicates = remove_near_duplicates,
  jaccard_cutoff = jaccard_cutoff
)

top_cc <- make_enrichment_plot(
  df = cc_df,
  n_terms = n_cc_terms,
  title_text = paste0(
    "Top Nonredundant Enriched Cellular Components (BH <= ", bh_cutoff, ")"
  ),
  out_file = "top_nonredundant_cellular_components.pdf",
  bh_cutoff = bh_cutoff,
  overlap_col = overlap_col,
  remove_near_duplicates = remove_near_duplicates,
  jaccard_cutoff = jaccard_cutoff
)

# ---------------------------------------------------------
# 8. WRITE TABLES USED FOR EACH PLOT
# ---------------------------------------------------------

if (!is.null(top_pathways)) {
  write.csv(
    top_pathways,
    file.path(outdir, "top_nonredundant_pathways.csv"),
    row.names = FALSE
  )
}

if (!is.null(top_bp)) {
  write.csv(
    top_bp,
    file.path(outdir, "top_nonredundant_biological_processes.csv"),
    row.names = FALSE
  )
}

if (!is.null(top_cc)) {
  write.csv(
    top_cc,
    file.path(outdir, "top_nonredundant_cellular_components.csv"),
    row.names = FALSE
  )
}