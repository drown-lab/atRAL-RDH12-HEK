library(dplyr)
library(purrr)
library(stringr)
library(tibble)

# ----------------------------
# helper: safe filenames
# ----------------------------
safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

# ----------------------------
# function: compute all Venn regions for N sets
# ----------------------------
compute_venn_regions <- function(sets_named) {
  stopifnot(is.list(sets_named), length(sets_named) >= 2)
  
  set_names <- names(sets_named)
  if (is.null(set_names) || any(set_names == "")) {
    stop("prot_sets must be a named list")
  }
  
  # universe of all IDs
  universe <- unique(unlist(sets_named, use.names = FALSE))
  
  # membership matrix: rows = IDs, cols = sets
  mem <- sapply(sets_named, function(s) universe %in% s)
  mem <- as.data.frame(mem, check.names = FALSE)  # keep exact names
  mem$id <- universe
  
  # key per ID like "1010" (present/absent pattern)
  key <- apply(mem[, set_names, drop = FALSE], 1, function(v) paste0(as.integer(v), collapse = ""))
  
  # all possible non-empty patterns (exclude all-zeros)
  patterns <- expand.grid(rep(list(c(0, 1)), length(set_names)))
  colnames(patterns) <- set_names                      # <-- FIX
  patterns <- patterns[apply(patterns, 1, sum) > 0, , drop = FALSE]
  patterns$key <- apply(patterns[, set_names, drop = FALSE], 1, paste0, collapse = "")
  
  # label patterns like "A&B" or "A_only"
  pattern_label <- function(row) {
    present <- set_names[as.logical(as.integer(row))]
    if (length(present) == 1) paste0(present, "_only") else paste(present, collapse = "&")
  }
  patterns$label <- apply(patterns[, set_names, drop = FALSE], 1, pattern_label)
  
  # build region lists
  region_lists <- lapply(seq_len(nrow(patterns)), function(i) {
    mem$id[key == patterns$key[i]]
  })
  names(region_lists) <- patterns$label
  
  # counts table
  region_counts <- tibble(
    region = names(region_lists),
    n_ids  = lengths(region_lists)
  ) %>%
    arrange(desc(n_ids))
  
  list(regions = region_lists, counts = region_counts)
}

# ----------------------------
# run + export
# ----------------------------
out_dir <- "RetinaOrganoids/Combined_HekandOrganoid_scripts/venn_regions"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

venn_obj <- compute_venn_regions(prot_sets)
regions  <- venn_obj$regions
counts   <- venn_obj$counts

# write a summary of region sizes
write.csv(counts, file.path(out_dir, "Venn_region_counts.csv"), row.names = FALSE)

# write each region’s IDs
walk(names(regions), function(rg) {
  out <- tibble(ID = regions[[rg]])
  fn <- paste0("IDs_", safe_name(rg), ".csv")
  write.csv(out, file.path(out_dir, fn), row.names = FALSE)
})

message("Exported ", length(regions), " region files + region count table to: ", out_dir)
print(head(counts, 10))
