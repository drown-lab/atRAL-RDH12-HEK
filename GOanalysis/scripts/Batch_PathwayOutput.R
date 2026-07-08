library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(readr)
library(paletteer)
library(ReactomePA)

# =========================================================
# 1. USER INPUTS
# =========================================================

input_dir <- "Proteomic_output_txts/PANGEA_exports"
output_root <- "PathwayAnalysis/output"

# Choose one or both:
# pathway_databases <- c("Reactome")
# pathway_databases <- c("KEGG")
pathway_databases <- c("Reactome", "KEGG")

# Optional terminal override, for example:
# Sys.setenv(PATHWAY_DATABASES = "Reactome")
# Sys.setenv(PATHWAY_DATABASES = "KEGG")
pathway_database_override <- Sys.getenv("PATHWAY_DATABASES", unset = "")
if (pathway_database_override != "") {
  pathway_databases <- str_split(pathway_database_override, ",")[[1]] %>%
    str_trim() %>%
    .[. != ""]
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Enrichment settings
pvalue_cutoff <- 0.2
qvalue_cutoff <- 0.2
min_gs_size <- 3
max_gs_size <- 700

# Optional manual filtering
use_manual_filter <- TRUE
min_count_keep <- 3
max_terms_to_plot <- 15
max_padj_keep <- 0.10
terms_to_plot <- 15

terms_to_remove_exact <- character(0)
patterns_to_remove <- character(0)

# =========================================================
# 2. HELPER FUNCTIONS
# =========================================================

clean_vector <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    str_remove("-\\d+$") %>%
    .[. != ""] %>%
    unique()
}

detect_id_type <- function(df) {
  if ("Protein" %in% colnames(df)) return("SYMBOL")
  if ("Gene" %in% colnames(df)) return("SYMBOL")
  if ("ID" %in% colnames(df)) return("UNIPROT")
  stop("Could not detect identifier column. Expected Protein, Gene, or ID.")
}

extract_ids <- function(df, id_type) {
  if (id_type == "SYMBOL") {
    if ("Protein" %in% colnames(df)) return(clean_vector(df$Protein))
    if ("Gene" %in% colnames(df)) return(clean_vector(df$Gene))
  }
  if (id_type == "UNIPROT") {
    return(clean_vector(df$ID))
  }
  stop("Unsupported id_type")
}

map_ids_to_entrez <- function(ids, id_type) {
  bitr(
    geneID = ids,
    fromType = id_type,
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db,
    drop = TRUE
  )
}

make_manual_table <- function(enrich_df) {
  enrich_df2 <- enrich_df %>%
    mutate(
      Description_lower = str_to_lower(Description),
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      FoldEnrichment_num = as.numeric(FoldEnrichment)
    )
  
  if (use_manual_filter) {
    enrich_df2 <- enrich_df2 %>%
      filter(Count >= min_count_keep) %>%
      filter(p.adjust <= max_padj_keep) %>%
      filter(!Description %in% terms_to_remove_exact)
    
    if (length(patterns_to_remove) > 0) {
      combined_pattern <- paste(patterns_to_remove, collapse = "|")
      enrich_df2 <- enrich_df2 %>%
        filter(!str_detect(Description_lower, regex(combined_pattern, ignore_case = TRUE)))
    }
  }
  
  enrich_df2 %>%
    arrange(p.adjust, desc(Count), desc(FoldEnrichment_num)) %>%
    slice_head(n = max_terms_to_plot) %>%
    mutate(
      Description = fct_reorder(Description, FoldEnrichment_num)
    )
}

make_removed_terms_table <- function(enrich_df) {
  enrich_df %>%
    mutate(
      Description_lower = str_to_lower(Description),
      removed_exact = Description %in% terms_to_remove_exact,
      removed_pattern = if (length(patterns_to_remove) > 0) {
        str_detect(Description_lower, regex(paste(patterns_to_remove, collapse = "|"), ignore_case = TRUE))
      } else {
        FALSE
      },
      removed_low_count = Count < min_count_keep,
      removed_high_padj = p.adjust > max_padj_keep
    ) %>%
    filter(removed_exact | removed_pattern | removed_low_count | removed_high_padj)
}

make_plot_table <- function(enrich_df, n_terms) {
  enrich_df %>%
    mutate(
      FoldEnrichment_num = as.numeric(FoldEnrichment)
    ) %>%
    arrange(p.adjust, desc(Count), desc(FoldEnrichment_num)) %>%
    slice_head(n = n_terms) %>%
    mutate(
      Description = fct_reorder(Description, FoldEnrichment_num)
    )
}

save_pathway_dotplot <- function(plot_df, title, out_prefix, width = 10, height = 6) {
  if (nrow(plot_df) == 0) return(invisible(NULL))
  
  p_dot <- ggplot(
    plot_df,
    aes(x = FoldEnrichment_num, y = Description, size = Count, color = p.adjust)
  ) +
    geom_point() +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 12)
    ) +
    labs(
      title = title,
      x = "Fold enrichment",
      y = NULL,
      color = "p.adjust"
    ) +
    scale_color_paletteer_c("grDevices::Plasma", direction = -1)
  
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = p_dot,
    width = width,
    height = height
  )
  
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p_dot,
    width = width,
    height = height,
    dpi = 900
  )
}

run_pathway_enrichment <- function(database, sig_entrez, universe_entrez) {
  if (database == "Reactome") {
    return(
      enrichPathway(
        gene = sig_entrez,
        universe = universe_entrez,
        organism = "human",
        pvalueCutoff = pvalue_cutoff,
        pAdjustMethod = "BH",
        qvalueCutoff = qvalue_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = TRUE
      )
    )
  }
  
  if (database == "KEGG") {
    ekegg <- enrichKEGG(
      gene = sig_entrez,
      universe = universe_entrez,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = "BH",
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size
    )
    
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      ekegg <- setReadable(
        ekegg,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID"
      )
    }
    
    return(ekegg)
  }
  
  stop("Unsupported pathway database: ", database)
}

run_one_enrichment <- function(sig_file, all_file, database, output_root) {
  contrast_name <- tools::file_path_sans_ext(basename(sig_file))
  out_dir <- file.path(output_root, paste0(database, "_batch"), contrast_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n====================================================\n")
  cat("Processing:", contrast_name, "\n")
  cat("Database:", database, "\n")
  cat("sig_file:", sig_file, "\n")
  cat("all_file:", all_file, "\n")
  
  df_sig <- read.csv(sig_file, stringsAsFactors = FALSE, check.names = FALSE)
  df_all <- read.csv(all_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  id_type_sig <- detect_id_type(df_sig)
  id_type_all <- detect_id_type(df_all)
  
  if (id_type_sig != id_type_all) {
    if (id_type_sig == "SYMBOL" && ("Gene" %in% colnames(df_all))) {
      id_type_all <- "SYMBOL"
    } else if (id_type_all == "SYMBOL" && ("Gene" %in% colnames(df_sig) || "Protein" %in% colnames(df_sig))) {
      id_type_sig <- "SYMBOL"
    } else {
      stop("ID type mismatch between sig and all files that cannot be reconciled.")
    }
  }
  
  sig_ids <- extract_ids(df_sig, id_type_sig)
  universe_ids <- extract_ids(df_all, id_type_all)
  
  if (id_type_sig == "SYMBOL" && id_type_all == "SYMBOL") {
    sig_ids <- intersect(sig_ids, universe_ids)
  }
  
  cat("Detected ID type:", id_type_sig, "\n")
  cat("Unique sig IDs:", length(sig_ids), "\n")
  cat("Unique universe IDs:", length(universe_ids), "\n")
  
  if (length(sig_ids) == 0) {
    warning("No significant IDs found after cleaning/intersection for ", contrast_name)
    return(NULL)
  }
  
  id_map <- map_ids_to_entrez(universe_ids, id_type_all)
  
  write.csv(
    id_map,
    file.path(out_dir, "id_to_entrez_symbol_map.csv"),
    row.names = FALSE
  )
  
  if (id_type_sig == "SYMBOL") {
    sig_entrez <- id_map %>%
      filter(SYMBOL %in% sig_ids) %>%
      pull(ENTREZID) %>%
      unique()
  } else if (id_type_sig == "UNIPROT") {
    if (!"UNIPROT" %in% colnames(id_map)) {
      id_map_sig <- bitr(
        geneID = sig_ids,
        fromType = "UNIPROT",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db,
        drop = TRUE
      )
      sig_entrez <- unique(id_map_sig$ENTREZID)
    } else {
      sig_entrez <- id_map %>%
        filter(UNIPROT %in% sig_ids) %>%
        pull(ENTREZID) %>%
        unique()
    }
  } else {
    stop("Unsupported sig ID type")
  }
  
  universe_entrez <- unique(id_map$ENTREZID)
  
  cat("Mapped sig Entrez IDs:", length(sig_entrez), "\n")
  cat("Mapped universe Entrez IDs:", length(universe_entrez), "\n")
  
  if (length(sig_entrez) == 0 || length(universe_entrez) == 0) {
    warning("No mapped Entrez IDs for ", contrast_name)
    return(NULL)
  }
  
  pathway_result <- run_pathway_enrichment(database, sig_entrez, universe_entrez)
  
  if (is.null(pathway_result) || nrow(as.data.frame(pathway_result)) == 0) {
    cat("No", database, "pathways enriched for", contrast_name, "\n")
    writeLines(
      paste("No", database, "pathways enriched."),
      file.path(out_dir, "NO_ENRICHMENT.txt")
    )
    return(
      data.frame(
        database = database,
        contrast = contrast_name,
        status = "No enrichment",
        n_sig_ids = length(sig_ids),
        n_universe_ids = length(universe_ids),
        n_sig_entrez = length(sig_entrez),
        n_universe_entrez = length(universe_entrez),
        n_raw_terms = NA_integer_,
        n_manual_terms = NA_integer_
      )
    )
  }
  
  pathway_df <- as.data.frame(pathway_result)
  pathway_manual_df <- make_manual_table(pathway_df)
  removed_terms_df <- make_removed_terms_table(pathway_df)
  pathway_plot_df <- make_plot_table(pathway_df, terms_to_plot)
  
  write.csv(
    pathway_df,
    file.path(out_dir, paste0(database, "_raw.csv")),
    row.names = FALSE
  )
  
  write.csv(
    pathway_manual_df,
    file.path(out_dir, paste0(database, "_manual_filtered.csv")),
    row.names = FALSE
  )
  
  write.csv(
    removed_terms_df,
    file.path(out_dir, paste0(database, "_removed_terms.csv")),
    row.names = FALSE
  )
  
  save_pathway_dotplot(
    plot_df = pathway_plot_df,
    title = paste0(contrast_name, ": ", database, " pathways"),
    out_prefix = file.path(out_dir, paste0(database, "_dotplot"))
  )
  
  save_pathway_dotplot(
    plot_df = pathway_manual_df,
    title = paste0(contrast_name, ": ", database, " pathways (manual filtered)"),
    out_prefix = file.path(out_dir, paste0(database, "_manual_filtered_dotplot"))
  )
  
  cat("Finished:", contrast_name, database, "\n")
  
  data.frame(
    database = database,
    contrast = contrast_name,
    status = "OK",
    n_sig_ids = length(sig_ids),
    n_universe_ids = length(universe_ids),
    n_sig_entrez = length(sig_entrez),
    n_universe_entrez = length(universe_entrez),
    n_raw_terms = nrow(pathway_df),
    n_manual_terms = nrow(pathway_manual_df)
  )
}

# =========================================================
# 3. FIND ALL FILE PAIRS
# =========================================================

all_files <- list.files(
  input_dir,
  pattern = "_all\\.csv$",
  full.names = TRUE
)

sig_files <- list.files(
  input_dir,
  pattern = "_(up|down)\\.csv$",
  full.names = TRUE
)

if (length(all_files) == 0) stop("No _all.csv files found.")
if (length(sig_files) == 0) stop("No _up.csv or _down.csv files found.")

find_matching_all <- function(sig_file, all_files) {
  sig_base <- basename(sig_file) %>%
    str_remove("_(up|down)\\.csv$")
  candidates <- all_files[basename(all_files) == paste0(sig_base, "_all.csv")]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

pair_df <- data.frame(
  sig_file = sig_files,
  all_file = vapply(sig_files, find_matching_all, character(1), all_files = all_files),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(all_file))

if (nrow(pair_df) == 0) {
  stop("No matching _all.csv files found for the _up/_down files.")
}

cat("\nFound", nrow(pair_df), "matched input contrasts.\n")

# =========================================================
# 4. RUN ALL JOBS
# =========================================================

results_list <- list()

for (database in pathway_databases) {
  if (!database %in% c("Reactome", "KEGG")) {
    stop("Unsupported database in pathway_databases: ", database)
  }
  
  database_results <- vector("list", nrow(pair_df))
  
  for (i in seq_len(nrow(pair_df))) {
    res <- tryCatch(
      run_one_enrichment(
        sig_file = pair_df$sig_file[i],
        all_file = pair_df$all_file[i],
        database = database,
        output_root = output_root
      ),
      error = function(e) {
        data.frame(
          database = database,
          contrast = tools::file_path_sans_ext(basename(pair_df$sig_file[i])),
          status = paste("ERROR:", conditionMessage(e))
        )
      }
    )
    database_results[[i]] <- res
  }
  
  database_summary_df <- bind_rows(database_results)
  
  write.csv(
    database_summary_df,
    file.path(output_root, paste0(database, "_batch"), paste0(database, "_batch_summary.csv")),
    row.names = FALSE
  )
  
  results_list[[database]] <- database_summary_df
}

summary_df <- bind_rows(results_list)

write.csv(
  summary_df,
  file.path(output_root, "pathway_batch_summary.csv"),
  row.names = FALSE
)

cat("\nPathway batch run complete.\n")
print(summary_df)
