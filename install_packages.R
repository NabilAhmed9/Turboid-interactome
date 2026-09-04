# install_packages.R
# Install all dependencies for turboid_interactome.R.
# Run once:  source("install_packages.R")

cran <- c("tidyverse", "matrixStats", "cowplot", "ggpubr", "ggrepel",
          "DT", "circlize", "gprofiler2", "tidytext")
bioc <- c("limma", "SummarizedExperiment", "ComplexHeatmap",
          "org.Hs.eg.db", "GO.db", "AnnotationDbi")

installed <- rownames(installed.packages())

need_cran <- setdiff(cran, installed)
if (length(need_cran)) install.packages(need_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
need_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(need_bioc)) BiocManager::install(need_bioc, update = FALSE, ask = FALSE)

missing <- setdiff(c(cran, bioc), rownames(installed.packages()))
if (length(missing))
  stop("Not installed: ", paste(missing, collapse = ", "))
message("All dependencies present.")
