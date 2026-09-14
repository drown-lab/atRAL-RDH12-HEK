# Plot DIA-NN extracted ion chromatograms (XICs) for selected protein groups
# across all 27 Exploris runs: one page per protein group, 9 conditions (rows) x
# 3 biological replicates (columns). Fragment traces are overlaid per panel, the
# DIA-NN peak boundaries are shaded, and each panel is labelled with the
# identification status in the ORIGINAL (deposited) report.
#
# Inputs
#   xic_dir      : DIANN_xic_rerun/report_xic/  (one <run>.xic.parquet per run; long format:
#                  pr = Precursor.Id, feature = "ms1" | "y6^1" | "b3^1" | "index", rt (min), value)
#   report_orig  : original report.parquet (analysis of record, 1% FDR filtered)
#   report_new   : DIANN_xic_rerun/report.parquet (--qvalue 1 -> every precursor's best candidate, q-value, RT.Start/Stop)
#   peptide_tbl  : SI_Table_Protein_PeptideCounts_Exploris.csv (to pick single-peptide protein groups)
#
# Usage: edit the paths, then
#   source("plot_xics_single_peptide.R")
#   plot_curated()   # 13 single-peptide + 4 two-peptide named proteins -> XIC_curated.pdf
#   plot_bulk()      # all single-peptide protein groups -> XIC_all_single_peptide_PGs.pdf (479 pages)

suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(purrr)
})

xic_dir     <- "F:/MS_Temp/atRAL_manuscript/proteomics/DIANN_xic_rerun/report_xic"
report_orig <- "F:/MS_Temp/atRAL_manuscript/proteomics/report.parquet"
report_new  <- "F:/MS_Temp/atRAL_manuscript/proteomics/DIANN_xic_rerun/report.parquet"
peptide_tbl <- "C:/Users/bryon/atRAL-Stress-Rams-Collab/manuscript/SI_materials/SI_Tables/SI_Table_Protein_PeptideCounts_Exploris.csv"
out_dir     <- "C:/Users/bryon/atRAL-Stress-Rams-Collab/manuscript/SI_materials/XICs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- sample annotation from run names ---------------------------------------
sample_from_run <- function(run) {
  run <- basename(run) |> str_remove("\\.raw$")
  acute <- str_match(run, "Hek_(Control|100uMatRAL|200uMatRAL)_Rep([ABC])")
  recov <- str_match(run, "Hek(GFP|RDH12)_(EtOH|100uM|200uM)_Rep(\\d)")
  tibble(run = run,
    genotype  = ifelse(!is.na(acute[,1]), "RDH12", ifelse(recov[,2] == "GFP", "WT", "RDH12")),
    phase     = ifelse(!is.na(acute[,1]), "acute", "recovery"),
    treatment = ifelse(!is.na(acute[,1]),
                       c(Control = "Veh", `100uMatRAL` = "100 µM", `200uMatRAL` = "200 µM")[acute[,2]],
                       c(EtOH = "Veh", `100uM` = "100 µM", `200uM` = "200 µM")[recov[,3]]),
    rep       = ifelse(!is.na(acute[,1]), match(acute[,3], c("A","B","C")), as.integer(recov[,4]))) |>
    mutate(condition = paste(genotype, treatment, phase),
           condition = factor(condition, levels = c(
             "RDH12 Veh acute", "RDH12 100 µM acute", "RDH12 200 µM acute",
             "RDH12 Veh recovery", "RDH12 100 µM recovery", "RDH12 200 µM recovery",
             "WT Veh recovery", "WT 100 µM recovery", "WT 200 µM recovery")))
}

xic_files <- list.files(xic_dir, pattern = "\\.xic\\.parquet$", full.names = TRUE)
stopifnot(length(xic_files) == 27)
runs <- sample_from_run(str_remove(basename(xic_files), "\\.xic\\.parquet$")) |> mutate(file = xic_files)

# ---- identification status --------------------------------------------------
cols <- c("Run", "Precursor.Id", "Stripped.Sequence", "Proteotypic", "Protein.Group", "Protein.Ids", "Genes", "Q.Value", "PG.Q.Value",
          "Lib.Q.Value", "Lib.PG.Q.Value", "Channel.Q.Value", "RT", "RT.Start", "RT.Stop", "Precursor.Quantity")
orig <- read_parquet(report_orig, col_select = all_of(cols)) |> mutate(Run = basename(Run) |> str_remove("\\.raw$"))
newr <- read_parquet(report_new,  col_select = all_of(cols)) |> mutate(Run = basename(Run) |> str_remove("\\.raw$"))

# precursor -> protein group map (from the original report), PROTEOTYPIC precursors only -
# the same filter the DEP pipeline applied, so the page shows exactly the evidence that was quantified.
# Shared (non-proteotypic) peptides that DIA-NN lists under the group are deliberately excluded.
# The full V12 filter chain (00_Process_DIANNparquetfile.R) is applied so peptide counts match SI_Table_Protein_Support.csv.
pr_map <- orig |>
  filter(Q.Value <= 0.01, PG.Q.Value <= 0.05, Lib.Q.Value <= 0.01, Lib.PG.Q.Value <= 0.01, Channel.Q.Value <= 0.05,
         Proteotypic == 1, !grepl("cRAP", Protein.Ids, ignore.case = TRUE)) |>
  distinct(Precursor.Id, Stripped.Sequence, Protein.Group, Genes)

# ---- core: fetch XICs for a set of precursors across all runs ---------------
fetch_xics <- function(prs) {
  map_dfr(seq_len(nrow(runs)), function(i) {
    open_dataset(runs$file[i]) |>
      filter(pr %in% prs, feature != "index") |>
      select(pr, feature, rt, value) |>
      collect() |>
      mutate(run = runs$run[i])
  })
}

# ---- one page per PRECURSOR (a two-peptide protein group therefore gets two pages) ----
plot_precursor_page <- function(pg, pr_id, xic, k, n, npep_note = "") {
  gene <- pr_map |> filter(Protein.Group == pg) |> pull(Genes) |> unique() |> paste(collapse = ";")
  d <- xic |> filter(pr == pr_id) |> left_join(runs, by = "run")
  if (nrow(d) == 0) return(NULL)
  # per-run status of this precursor
  st <- runs |> select(run, condition, rep) |>
    mutate(Precursor.Id = pr_id) |>
    left_join(orig |> select(Run, Precursor.Id, q_orig = Q.Value, s0 = RT.Start, e0 = RT.Stop), by = c(run = "Run", "Precursor.Id")) |>
    left_join(newr |> select(Run, Precursor.Id, q_new = Q.Value, s1 = RT.Start, e1 = RT.Stop),  by = c(run = "Run", "Precursor.Id")) |>
    mutate(identified = !is.na(q_orig),
           start = ifelse(identified, s0, s1), stop = ifelse(identified, e0, e1),
           label = ifelse(identified, sprintf("q = %.1e", q_orig),
                          ifelse(is.na(q_new), "no candidate", sprintf("not ID'd (q = %.2f)", q_new))))
  # y-scale each panel to the fragment maximum INSIDE the DIA-NN peak boundaries (+/- 0.1 min), so the scored peak
  # is legible even when a single interfering fragment dominates elsewhere in the window (common for weak IDs).
  # Traces above the limit are clipped; absent runs therefore show baseline noise at the candidate position.
  d <- d |> mutate(ms1 = feature == "ms1")
  frag <- d |> filter(!ms1)
  win  <- st |> select(run, start, stop)
  ymax <- frag |> left_join(win, by = "run") |>
    group_by(run) |>
    summarise(ymax = { inw <- !is.na(start) & rt >= start - 0.1 & rt <= stop + 0.1
                       m <- if (any(inw)) max(value[inw]) else max(value); max(m, 1) * 1.15 }, .groups = "drop")
  frag <- frag |> left_join(ymax, by = "run") |> mutate(value = pmin(value, ymax))
  ms1  <- d |> filter(ms1) |> left_join(ymax, by = "run") |> group_by(run) |>
    mutate(value = value / max(value, 1) * ymax * 0.8) |> ungroup()
  st   <- st |> left_join(ymax, by = "run") |> mutate(ymax = coalesce(ymax, 1))

  ggplot() +
    geom_rect(data = st |> filter(!is.na(start)),
              aes(xmin = start, xmax = stop, ymin = 0, ymax = Inf, fill = identified), alpha = 0.15) +
    geom_line(data = ms1, aes(rt, value), colour = "grey55", linetype = "22", linewidth = 0.3) +
    geom_line(data = frag, aes(rt, value, colour = feature), linewidth = 0.35) +
    geom_text(data = st, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.4, size = 2.3,
              colour = ifelse(st$identified, "black", "firebrick")) +
    facet_grid(condition ~ rep, scales = "free", labeller = labeller(rep = function(x) paste("rep", x))) +
    scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "firebrick"), guide = "none") +
    scale_y_continuous(labels = scales::label_scientific(digits = 1), expand = expansion(mult = c(0, 0.25))) +
    labs(title = sprintf("%s (%s)  -  precursor %d of %d: %s", gene, pg, k, n, pr_id),
         subtitle = str_wrap(paste0("Fragment XICs (coloured; y-axis scaled to the in-boundary maximum, off-peak interference clipped); MS1 (grey dashed, rescaled). ",
                           "Shaded = DIA-NN peak boundaries: blue, identified at 1% FDR in the deposited report; red, best candidate in a run where it was not identified. ", npep_note), width = 150),
         x = "Retention time (min)", y = "Intensity", colour = "Fragment") +
    theme_minimal(base_size = 8) +
    theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
          strip.text.y = element_text(angle = 0), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey80"))
}

# all pages for a protein group: one per proteotypic precursor, ordered by peptide then charge
plot_protein_page <- function(pg, xic, npep_note = "") {
  prs <- pr_map |> filter(Protein.Group == pg) |> arrange(Stripped.Sequence, Precursor.Id) |> pull(Precursor.Id) |> unique()
  npep <- pr_map |> filter(Protein.Group == pg) |> pull(Stripped.Sequence) |> n_distinct()
  note <- paste0(sprintf("Protein group has %d proteotypic peptide%s (%d precursor%s). ", npep, ifelse(npep == 1, "", "s"),
                         length(prs), ifelse(length(prs) == 1, "", "s")), npep_note)
  pages <- lapply(seq_along(prs), function(k) plot_precursor_page(pg, prs[k], xic, k, length(prs), note))
  Filter(Negate(is.null), pages)
}

render_pages <- function(pgs, file, notes = NULL) {
  prs <- pr_map |> filter(Protein.Group %in% pgs) |> pull(Precursor.Id) |> unique()
  message("Fetching XICs for ", length(prs), " precursors / ", length(pgs), " protein groups across 27 runs ...")
  xic <- fetch_xics(prs)
  pdf(file.path(out_dir, file), width = 11, height = 14)
  on.exit(dev.off())
  for (pg in pgs) {
    note <- if (!is.null(notes) && pg %in% names(notes)) notes[[pg]] else ""
    pages <- plot_protein_page(pg, xic, note)
    if (length(pages)) for (p in pages) print(p) else message("no XIC rows for ", pg)
  }
  message("Wrote ", file.path(out_dir, file))
}

# ---- curated set: single-peptide + sparse two-peptide proteins named in the text/figures
curated <- c(
  # single-peptide, named in text
  MGST2 = "Q99735", SLC11A2 = "P49281", NFE2L2 = "Q16236", SLC39A7 = "Q92504",
  AGPAT3 = "Q9NRZ7", BCL2L12 = "Q9HB09", PIDD1 = "Q9HB75",
  # single-peptide, Figure 3H only
  TUBB2A = "Q13885", TUBB4B = "P68371", SLC30A5 = "Q8TAD4", SLC31A1 = "O15431", SLC39A1 = "Q9NY26", ALG3 = "Q92685",
  # two-peptide, sparsely detected, named in abstract/discussion
  SLC7A11 = "Q9UPY5", GSTA4 = "O15217", CHKA = "P35790", CYBA = "P13498")

plot_curated <- function() render_pages(unname(curated), "XIC_curated_named_proteins.pdf")

# ---- bulk set: every single-peptide protein group in the peptide-count table
plot_bulk <- function() {
  pc  <- read.csv(peptide_tbl)
  pgs <- pc$Protein_Group[pc$Single_peptide_protein_group %in% c(TRUE, "TRUE", "true")]
  # chunk so each fetch stays modest in memory
  chunks <- split(pgs, ceiling(seq_along(pgs) / 60))
  pdf(file.path(out_dir, "XIC_all_single_peptide_PGs.pdf"), width = 11, height = 14)
  on.exit(dev.off())
  for (ch in chunks) {
    prs <- pr_map |> filter(Protein.Group %in% ch) |> pull(Precursor.Id) |> unique()
    xic <- fetch_xics(prs)
    for (pg in ch) for (p in plot_protein_page(pg, xic)) print(p)
    message("done ", length(ch), " protein groups")
  }
}
