# One-time environment setup for this project (R 4.4.1 + renv)
# Run from the project root; safe to re-run if interrupted.

# Pin Bioconductor to 3.20 — the release matching R 4.4.1 and the versions
# in dependciesAndPackages_info.txt (Biobase 2.66 / GenomeInfoDb 1.42).
renv::settings$bioconductor.version("3.20")

renv::install("BiocManager", prompt = FALSE)

# 'Deriv' (a dependency via doBy/rstatix) was archived from CRAN, which makes
# renv's resolver crash with "failed to create directory at path ''".
# Install it up front so renv sees it already satisfied. A Windows binary is
# still on CRAN; fall back to the archived source tarball (pure R, no
# compilation) if the binary disappears.
if (!requireNamespace("Deriv", quietly = TRUE)) {
  tryCatch(
    install.packages("Deriv", type = "binary"),
    warning = function(w) install.packages(
      "https://cran.r-project.org/src/contrib/Archive/Deriv/Deriv_4.2.6.tar.gz",
      repos = NULL, type = "source"
    )
  )
}

pkgs <- c(
  # CRAN
  "tidyverse", "readxl", "writexl", "Peptides", "UpSetR", "arrow",
  "ggthemes", "reshape2", "ggridges", "ggpubr", "outliers", "paletteer",
  "circlize", "GGally", "ComplexUpset", "ggVennDiagram", "ggrepel",
  "ggsci", "janitor", "patchwork", "pheatmap", "viridis", "RColorBrewer",
  "scales", "minpack.lm", "colorBlindness", "jsonlite",
  # Bioconductor (3.20)
  "bioc::DEP", "bioc::PRONE", "bioc::SummarizedExperiment",
  "bioc::ComplexHeatmap", "bioc::vsn", "bioc::sva",
  "bioc::clusterProfiler", "bioc::org.Hs.eg.db", "bioc::AnnotationDbi",
  "bioc::GOSemSim", "bioc::enrichplot", "bioc::ReactomePA"
)

renv::install(pkgs, prompt = FALSE)
