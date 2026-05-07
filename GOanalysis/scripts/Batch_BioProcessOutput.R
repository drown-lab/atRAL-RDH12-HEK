library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(enrichplot)
library(ggplot2)
library(forcats)
library(GOSemSim)
library(readr)
library(scales)
library(paletteer)

# =========================================================
# 1. USER INPUTS
# =========================================================

input_dir <- "Proteomic_output_txts/PANGEA_exports"
output_root <- "GOanalysis/output/clusterProfiler_batch"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Enrichment settings
pvalue_cutoff <- 0.3
qvalue_cutoff <- 0.3
min_gs_size <- 3
max_gs_size <- 700

# Simplify settings
simplify_cutoff <- 0.9

# Optional manual filtering
use_manual_filter <- TRUE
min_count_keep <- 3
max_terms_to_plot <- 15
max_padj_keep <- 0.20
terms_to_plot <- 10

terms_to_remove_exact <- character(0)
patterns_to_remove <- character(0)

use_go_level_filter <- FALSE
go_level_to_keep <- 4

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

make_manual_table <- function(ego_df) {
  ego_df2 <- ego_df
  
  if (use_manual_filter) {
    ego_df2 <- ego_df2 %>%
      mutate(
        Description_lower = str_to_lower(Description),
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
      ) %>%
      filter(Count >= min_count_keep) %>%
      filter(p.adjust <= max_padj_keep) %>%
      filter(!Description %in% terms_to_remove_exact)
    
    if (length(patterns_to_remove) > 0) {
      combined_pattern <- paste(patterns_to_remove, collapse = "|")
      ego_df2 <- ego_df2 %>%
        filter(!str_detect(Description_lower, regex(combined_pattern, ignore_case = TRUE)))
    }
    
    ego_df2 <- ego_df2 %>%
      arrange(p.adjust, desc(Count), desc(GeneRatio_num)) %>%
      slice_head(n = max_terms_to_plot) %>%
      mutate(
        Description = fct_reorder(Description, GeneRatio_num)
      )
  } else {
    ego_df2 <- ego_df2 %>%
      mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
        Description = fct_reorder(Description, GeneRatio_num)
      ) %>%
      arrange(p.adjust, desc(Count), desc(GeneRatio_num)) %>%
      slice_head(n = max_terms_to_plot)
  }
  
  ego_df2
}

run_one_enrichment <- function(sig_file, all_file, output_root) {
  contrast_name <- tools::file_path_sans_ext(basename(sig_file))
  out_dir <- file.path(output_root, contrast_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n====================================================\n")
  cat("Processing:", contrast_name, "\n")
  cat("sig_file:", sig_file, "\n")
  cat("all_file:", all_file, "\n")
  
  df_sig <- read.csv(sig_file, stringsAsFactors = FALSE, check.names = FALSE)
  df_all <- read.csv(all_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  id_type_sig <- detect_id_type(df_sig)
  id_type_all <- detect_id_type(df_all)
  
  # Prefer same identifier space; if mismatch, use SYMBOL if available
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
  
  ego_bp <- enrichGO(
    gene = sig_entrez,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    readable = TRUE
  )
  
  if (is.null(ego_bp) || nrow(as.data.frame(ego_bp)) == 0) {
    cat("No GO BP terms enriched for", contrast_name, "\n")
    writeLines("No GO BP terms enriched.", file.path(out_dir, "NO_ENRICHMENT.txt"))
    return(
      data.frame(
        contrast = contrast_name,
        status = "No enrichment",
        n_sig_ids = length(sig_ids),
        n_universe_ids = length(universe_ids),
        n_sig_entrez = length(sig_entrez),
        n_universe_entrez = length(universe_entrez)
      )
    )
  }
  
  hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP")
  
  ego_bp_s <- simplify(
    ego_bp,
    cutoff = simplify_cutoff,
    by = "p.adjust",
    select_fun = min,
    semData = hsGO
  )
  
  if (use_go_level_filter) {
    ego_bp_for_manual <- gofilter(ego_bp_s, level = go_level_to_keep)
  } else {
    ego_bp_for_manual <- ego_bp_s
  }
  
  ego_bp_df <- as.data.frame(ego_bp)
  ego_bp_s_df <- as.data.frame(ego_bp_s)
  ego_bp_manual_df <- make_manual_table(as.data.frame(ego_bp_for_manual))
  
  write.csv(ego_bp_df, file.path(out_dir, "enrichGO_BP_raw.csv"), row.names = FALSE)
  write.csv(ego_bp_s_df, file.path(out_dir, "enrichGO_BP_simplified.csv"), row.names = FALSE)
  write.csv(ego_bp_manual_df, file.path(out_dir, "enrichGO_BP_manual_filtered.csv"), row.names = FALSE)
  
  p_dot_simplified <- dotplot(ego_bp_s, showCategory = terms_to_plot) +
    ggtitle(paste0(contrast_name, ": GO BP (simplified)")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(out_dir, "enrichGO_BP_simplified_dotplot.pdf"),
    plot = p_dot_simplified,
    width = 9,
    height = 5
  )
  
  ggsave(
    filename = file.path(out_dir, "enrichGO_BP_simplified_dotplot.png"),
    plot = p_dot_simplified,
    width = 9,
    height = 5,
    dpi = 900
  )
  
  if (nrow(ego_bp_manual_df) > 0) {
    p_dot_manual <- ggplot(
      ego_bp_manual_df,
      aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)
    ) +
      geom_point() +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14)
      ) +
      labs(
        title = paste0(contrast_name, ": GO BP (manual filtered)"),
        x = "GeneRatio",
        y = NULL,
        color = "p.adjust"
      ) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_color_paletteer_c("grDevices::Plasma", direction = -1)
    
    ggsave(
      filename = file.path(out_dir, "enrichGO_BP_manual_filtered_dotplot.pdf"),
      plot = p_dot_manual,
      width = 9,
      height = 5
    )
    
    ggsave(
      filename = file.path(out_dir, "enrichGO_BP_manual_filtered_dotplot.png"),
      plot = p_dot_manual,
      width = 9,
      height = 5,
      dpi = 900
    )
  }
  
  cat("Finished:", contrast_name, "\n")
  
  data.frame(
    contrast = contrast_name,
    status = "OK",
    n_sig_ids = length(sig_ids),
    n_universe_ids = length(universe_ids),
    n_sig_entrez = length(sig_entrez),
    n_universe_entrez = length(universe_entrez),
    n_raw_terms = nrow(ego_bp_df),
    n_simplified_terms = nrow(ego_bp_s_df),
    n_manual_terms = nrow(ego_bp_manual_df)
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

cat("\nFound", nrow(pair_df), "matched enrichment jobs.\n")

# =========================================================
# 4. RUN ALL JOBS
# =========================================================

results_list <- vector("list", nrow(pair_df))

for (i in seq_len(nrow(pair_df))) {
  res <- tryCatch(
    run_one_enrichment(
      sig_file = pair_df$sig_file[i],
      all_file = pair_df$all_file[i],
      output_root = output_root
    ),
    error = function(e) {
      data.frame(
        contrast = tools::file_path_sans_ext(basename(pair_df$sig_file[i])),
        status = paste("ERROR:", conditionMessage(e))
      )
    }
  )
  results_list[[i]] <- res
}

summary_df <- bind_rows(results_list)

write.csv(
  summary_df,
  file.path(output_root, "clusterProfiler_batch_summary.csv"),
  row.names = FALSE
)

cat("\nBatch run complete.\n")
print(summary_df)