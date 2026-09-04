# turboid_interactome.R
# ============================================================================
# TurboID / BioID proximity-labeling interactome pipeline, implemented A-to-Z
# in R.
#
# Independent, fully scripted R reanalysis of TurboID
# proximity-proteomics dataset.

# The implementation that goes from the raw MaxQuant proteinGroups.txt all the way
# to differential proximity, evidence tiering, and GO / pathway enrichment.

#
# Evidence framing: measured quantitation and
# reproducible on/off detection are the two PRIMARY evidence streams. The
# imputed factorial model is a sensitive screen that defines the candidate set
# and a sensitivity analysis; imputation-only significance is the weakest tier.
# (The evidence_tier logic below ranks: measured > on/off > exploratory >
# imputation-only.)
#
# Inputs (only two):
#   data_raw/proteinGroups.txt   raw MaxQuant output. ;
#                                download from any raw data from PRIDE  and place here.
#   metadata/sample_sheet.csv    shows the LFQ columns "there is an example how to design the sample sheet".
#
# Outputs:
#   results/tables/*.csv
#   results/figures/*.{pdf,png,tiff}
#   results/*_TurboID_IRE1.summarizedExperiment.RDS
#   results/sessionInfo.txt
#
# Run:   Rscript turboid_interactome.R      (or source() interactively)
# Deps:  see install_packages.R and the library() calls in the setup section.
# Reuse: to run on another TurboID/BioID experiment, edit the parameters marked
#        "EDIT FOR YOUR DESIGN" in params_local and the contrast/on-off blocks.
# License: MIT (see LICENSE).
# ============================================================================

# chunk: setup
library(tidyverse)
library(matrixStats)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(DT)
library(limma)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
# Enrichment / GO plotting packages are loaded in the enrichment section:
#   pathways: gprofiler2 ; GO (direct annotations): org.Hs.eg.db, GO.db,
#   AnnotationDbi ; faceted GO bars: tidytext.
#   install.packages(c("gprofiler2", "tidytext"))
#   BiocManager::install(c("org.Hs.eg.db", "GO.db", "AnnotationDbi"))

set.seed(20230710)

# Bioconductor packages (SummarizedExperiment / ComplexHeatmap / org.Hs.eg.db via
# AnnotationDbi and S4Vectors) mask several dplyr verbs. Bind the dplyr versions
# into the global environment so they win over any later masking.
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
slice  <- dplyr::slice;  count  <- dplyr::count;  desc   <- dplyr::desc

# Everything the analysis depends on is defined here. Nothing downstream
# invents a path or a threshold. The single quantitative input is the raw
# MaxQuant proteinGroups.txt with an explicit sample sheet
params_local <- list(
  proteinGroups       = "data_raw/proteinGroups.txt",
  out_dir             = "results",
  min_unique_peptides = 2,
  # valid-value rule for the tested universe:
  #   "group3"  = valid in all 3 replicates of >=1 (construct x stimulation) group
  #   "merged6" = valid in all 6 replicates of >=1 construct   (reproduces paper)
  valid_rule          = "group3",
  impute_downshift    = 1.8,   # downshifted-normal, standard default
  impute_width        = 0.3,   # downshifted-normal, standard default
  fdr_cutoff          = 0.05,
  lfc_cutoff          = 1.0,
  min_obs_per_group   = 2,     # measured quantitative call needs >=2 observed values per group
  baits               = c("ERN1", "ERN2"),  # EDIT FOR YOUR DESIGN: bait gene symbols (IRE1a=ERN1, IRE1b=ERN2)
  biotin_carboxylases = c("PC","PCCA","PCCB","MCCC1","MCCC2","ACACA","ACACB","HLCS"),  # EDIT FOR YOUR DESIGN: endogenously biotinylated background
  # NULL for the published pipeline. If set, the appendix compares against it but
  # nothing in the analysis depends on it.
  meta_validation_csv = NULL,
  label_by            = "gene",              # plot labels: "gene" (short), "protein" (full), "uniprot" ** It's an example on IRE1
  aliases             = c(ERN1 = "IRE1\u03b1", ERN2 = "IRE1\u03b2"))  # display ERN1/2 as IRE1a/b

# Fixed output structure: every table goes to results/tables, every figure to
# results/figures. Both are created on every run.
params_local$tab_dir <- file.path(params_local$out_dir, "tables")
params_local$fig_dir <- file.path(params_local$out_dir, "figures")
invisible(lapply(c(params_local$out_dir, params_local$tab_dir, params_local$fig_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# Save a table to results/tables (returns the data invisibly so it can be piped).
save_table <- function(x, name) {
  readr::write_csv(x, file.path(params_local$tab_dir, paste0(name, ".csv")))
  invisible(x)
}

# Save a ggplot to results/figures at publication quality: vector PDF plus a
# 300-dpi PNG and a 300-dpi TIFF (LZW-compressed) for journal submission.
save_fig <- function(plot, name, width = 8, height = 6) {
  base <- file.path(params_local$fig_dir, name)
  ggplot2::ggsave(paste0(base, ".pdf"),  plot, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(paste0(base, ".png"),  plot, width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(base, ".tiff"), plot, width = width, height = height, dpi = 300,
                  compression = "lzw")
  invisible(plot)
}

# Save a ComplexHeatmap (needs an open device rather than ggsave).
save_heatmap <- function(ht, name, width = 8, height = 9) {
  base <- file.path(params_local$fig_dir, name)
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = width, height = height)
  ComplexHeatmap::draw(ht); grDevices::dev.off()
  grDevices::png(paste0(base, ".png"), width = width, height = height, units = "in", res = 300)
  ComplexHeatmap::draw(ht); grDevices::dev.off()
  grDevices::tiff(paste0(base, ".tiff"), width = width, height = height, units = "in",
                  res = 300, compression = "lzw")
  ComplexHeatmap::draw(ht); grDevices::dev.off()
}

# Save a base-graphics plot (e.g. a dendrogram) at publication quality.
save_base_plot <- function(fun, name, width = 8, height = 6) {
  base <- file.path(params_local$fig_dir, name)
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = width, height = height); fun(); grDevices::dev.off()
  grDevices::png(paste0(base, ".png"), width = width, height = height, units = "in", res = 300); fun(); grDevices::dev.off()
}

dt_table <- function(x) {
  DT::datatable(x, filter = "top",
                options = list(autoWidth = FALSE, scrollX = TRUE, pageLength = 15),
                class = "compact hover row-border stripe dt-left cell-border nowrap")
}

# Choose a display label for a set of UniProt accessions, per params_local$label_by.
# Default is the short gene symbol, with friendly aliases (ERN1 -> IRE1a). Full
# protein names are only used when label_by == "protein", and are truncated so
# plot labels stay legible.
label_for <- function(uniprotids) {
  i <- match(uniprotids, prot.annotation$uniprotid)
  gene <- prot.annotation$gene[i]
  gene <- ifelse(is.na(gene) | gene == "", uniprotids, gene)
  if (length(params_local$aliases)) {
    hit <- gene %in% names(params_local$aliases)
    gene[hit] <- params_local$aliases[gene[hit]]
  }
  switch(params_local$label_by,
         gene    = gene,
         uniprot = uniprotids,
         protein = {
           prot <- prot.annotation$protein[i]
           stringr::str_trunc(ifelse(is.na(prot) | prot == "", gene, prot), 45)
         },
         gene)
}

# Introduction ---------------------------------------------------------------
# This pipline processes the TurboID/BioID proximity-labeling proteomics 
# IRE1α and IRE1β are just examples from the raw MaxQuant that can be replaced with your own data
# `proteinGroups.txt` through to differential proximity analysis and GO
# enrichment. Everything is computed from the raw file: log2 transform, per-
# sample normalization, missing-value assessment, imputation and testing are
# all reconstructed here rather than taken from any pre-processed export.
# The two IRE1 constructs are both baits, not controls. IRE1β is a
# biologically related IRE1 isoform reference that, alongside the matched
# V5-TurboID background control, lets the analysis separate proteins shared
# by both isoforms from those selective for IRE1α or for IRE1β. The
# V5-TurboID construct is the background control; IRE1α and IRE1β are the
# proteins of interest.
# Two design facts govern the choices below and are stated up front.
# The first is that the cytosolic V5-TurboID control has the highest
# fraction of missing values, and the baits themselves illustrate why. IRE1α
# (ERN1) and IRE1β (ERN2) are detected at high intensity in their own
# construct and are entirely absent from the control. A protein absent from
# one arm cannot be tested from measured values alone, so a measured-only
# analysis silently discards exactly the most bait-specific proteins,
# including the baits. This is the reason proximity-labeling data are
# imputed with a downshifted normal: it gives on/off proteins a finite fold
# change. In this script the imputed matrix is used only as a sensitive
# screen that defines the candidate universe and recovers on/off bait-
# specific proteins; it is not the primary evidence. The two primary
# evidence streams are measured quantitation (limma on the measured matrix,
# >=2 observed values per group) and reproducible on/off detection (bait
# 3/3, matched control <=1/3). Imputation-only significance is a supporting
# sensitivity tier.
# The second is that the stimulation effect is small. DMSO and TM do not
# resolve into separate clusters. That is a reason to expect few
# stimulation-driven hits, not a reason to delete the stimulation factor, so
# the merged analysis (the manuscript's primary) and the factorial analysis
# (which keeps stimulation) are both run.

# Data carpentry and annotation ----------------------------------------------

# chunk: import
pg <- read.delim(params_local$proteinGroups, sep = "\t",
                 check.names = FALSE, stringsAsFactors = FALSE)

prot.annotation <- tibble(
  protein_ids     = pg[["Majority protein IDs"]],
  uniprotid       = sub(";.*", "", pg[["Majority protein IDs"]]),
  gene            = sub(";.*", "", pg[["Gene names"]]),
  protein         = pg[["Protein names"]],
  peptides        = pg[["Peptides"]],
  unique_peptides = pg[["Unique peptides"]],
  coverage        = pg[["Sequence coverage [%]"]],
  qvalue          = pg[["Q-value"]],
  score           = pg[["Score"]],
  n_accessions    = lengths(strsplit(ifelse(is.na(pg[["Majority protein IDs"]]), "",
                                            pg[["Majority protein IDs"]]), ";")),
  # a protein group with more than one accession cannot be assigned to a single
  # gene product from shared peptides alone; it is flagged, not dropped
  protein_group_ambiguous = lengths(strsplit(ifelse(is.na(pg[["Majority protein IDs"]]), "",
                                            pg[["Majority protein IDs"]]), ";")) > 1)

lfq_cols <- grep("^LFQ intensity ", colnames(pg), value = TRUE)
stopifnot(length(lfq_cols) == 18)
message(nrow(pg), " protein groups; ", length(lfq_cols), " LFQ columns.")

# Sample sheet ---------------------------------------------------------------
# Sample identity is read from a standalone file,
# `metadata/sample_sheet.csv`, which is a required "to be designed" input alongside
# `proteinGroups.txt`. It is the explicit bridge between the anonymous
# MaxQuant LFQ column names and their biological meaning, so a reviewer can
# inspect the experimental design without reading the code, and it is the
# single authoritative source for sample identity rather than a hardcoded
# object. The file has one row per sample with columns `sample_id`,
# `mq_lfq_column`, `construct`, `stimulation`, and `replicate`; the script
# derives `group_merged`, `group_full`, and the control flag internally, and
# validates that every `mq_lfq_column` exists in `proteinGroups.txt` before
# proceeding.

# chunk: sample-sheet
SS_PATH <- "metadata/sample_sheet.csv"
if (!file.exists(SS_PATH))
  stop("Required input not found: ", SS_PATH,
       "\nProvide a CSV with columns: sample_id, mq_lfq_column, construct, stimulation, replicate.")
sample_sheet <- readr::read_csv(SS_PATH, show_col_types = FALSE)

# validate structure
required_cols <- c("sample_id", "mq_lfq_column", "construct", "stimulation", "replicate")
missing_cols <- setdiff(required_cols, colnames(sample_sheet))
if (length(missing_cols))
  stop("sample_sheet.csv is missing columns: ", paste(missing_cols, collapse = ", "))

# validate the mapping against the actual MaxQuant columns
bad <- setdiff(sample_sheet$mq_lfq_column, lfq_cols)
if (length(bad))
  stop("These mq_lfq_column values are not present in proteinGroups.txt:\n  ",
       paste(bad, collapse = "\n  "),
       "\nColumn names must match the MaxQuant 'LFQ intensity ...' headers exactly.")
stopifnot(nrow(sample_sheet) == length(lfq_cols),
          !any(duplicated(sample_sheet$sample_id)),
          all(c("DMSO", "TM") %in% sample_sheet$stimulation))

# derive grouping internally (not stored in the file)
sample_sheet <- sample_sheet %>%
  mutate(group_merged = construct,
         group_full   = paste(construct, stimulation, sep = "_"),
         control      = construct == "con")
message("Sample sheet: ", nrow(sample_sheet), " samples, ",
        dplyr::n_distinct(sample_sheet$group_full), " construct x stimulation groups.")
sample_sheet %>% dt_table()

# Annotation  ------------------------------------------------------
# The identifier table: leading UniProt accession, gene, full protein name,
# and the peptide and score evidence for each identification. Nothing
# proceeds until these look right.

# chunk: annotation-checkpoint
prot.annotation %>%
  select(uniprotid, gene, protein, unique_peptides, peptides, coverage, qvalue, score) %>%
  dt_table()

# Filtering ------------------------------------------------------------------
# Filtering is staged and audited so the path from raw identifications to
# the tested universe is visible.

# chunk: filtering
is_blank <- function(x) is.na(x) | x == ""
audit <- tibble(step = "raw MaxQuant", n = nrow(pg))

pg.f <- pg %>%
  filter(is_blank(Reverse), is_blank(`Potential contaminant`),
         is_blank(`Only identified by site`))
audit <- add_row(audit, step = "remove reverse / contaminant / only-by-site", n = nrow(pg.f))

pg.f <- pg.f %>% filter(`Unique peptides` >= params_local$min_unique_peptides)
audit <- add_row(audit, step = paste0(">= ", params_local$min_unique_peptides,
                                       " unique peptides"), n = nrow(pg.f))

measured <- as.matrix(pg.f[, sample_sheet$mq_lfq_column])
colnames(measured) <- sample_sheet$sample_id
rownames(measured) <- sub(";.*", "", pg.f[["Majority protein IDs"]])
measured[measured == 0] <- NA          # zeros are missing, not zero abundance
measured <- log2(measured)

valid_in <- function(cols) rowSums(!is.na(measured[, cols, drop = FALSE])) == length(cols)
by_group     <- split(sample_sheet$sample_id, sample_sheet$group_full)
by_construct <- split(sample_sheet$sample_id, sample_sheet$construct)
keep_group3  <- Reduce(`|`, lapply(by_group,     valid_in))
keep_merged6 <- Reduce(`|`, lapply(by_construct, valid_in))
audit <- add_row(audit, step = "valid in all 3 reps of >=1 group (group3)",      n = sum(keep_group3))
audit <- add_row(audit, step = "valid in all 6 reps of >=1 construct (merged6)", n = sum(keep_merged6))

keep <- if (params_local$valid_rule == "merged6") keep_merged6 else keep_group3
measured <- measured[keep, ]
audit %>% dt_table()

# Normalization --------------------------------------------------------------
# MaxQuant LFQ already applies a between-sample normalization, so this step is a light per-sample median centering. It is shown per sample, not just in aggregate, so any single channel that is off-centre is visible.

# Per-sample distributions ---------------------------------------------------
# Each panel is one MS sample. The red line is the sample mean; the green
# line is mean, params_local$impute_downshift SD, which is where the imputation in the next section draws missing values from. Reading the two
# together shows both the shape of each sample's intensity distribution and where imputed values will sit relative to the measured bulk.

# chunk: per-sample-hist
col_med <- matrixStats::colMedians(measured, na.rm = TRUE)
measured.norm <- sweep(measured, 2, col_med - median(col_med, na.rm = TRUE))

sample_stats <- as_tibble(measured.norm, rownames = "uniprotid") %>%
  pivot_longer(-uniprotid, names_to = "sample", values_to = "log2") %>%
  group_by(sample) %>%
  summarise(mean = mean(log2, na.rm = TRUE),
            sd   = sd(log2, na.rm = TRUE), .groups = "drop") %>%
  mutate(downshift = mean - params_local$impute_downshift * sd,
         sample = factor(sample, levels = sample_sheet$sample_id))

p_persample <- as_tibble(measured.norm, rownames = "uniprotid") %>%
  pivot_longer(-uniprotid, names_to = "sample", values_to = "log2") %>%
  mutate(sample = factor(sample, levels = sample_sheet$sample_id)) %>%
  ggplot(aes(log2)) +
  geom_histogram(bins = 45, fill = "#9ecae1", colour = "white", linewidth = 0.1) +
  geom_vline(data = sample_stats, aes(xintercept = mean),
             colour = "red", linetype = 2) +
  geom_vline(data = sample_stats, aes(xintercept = downshift),
             colour = "darkgreen", linetype = 2) +
  facet_wrap(~ sample, ncol = 6, scales = "free_y") +
  cowplot::theme_minimal_grid(11) +
  labs(x = "log2 LFQ (median-centred)", y = "Count")
save_fig(p_persample, "QC_per_sample_histograms", width = 14, height = 8)
p_persample

sample_stats %>%
  transmute(sample, mean = round(mean, 2), sd = round(sd, 2),
            `mean-1.8SD` = round(downshift, 2)) %>%
  save_table("QC_per_sample_stats") %>%
  dt_table()

# Before and after centering the data -------------------------------------------------

# chunk: norm-box
p_normbox <- bind_rows(
  as_tibble(measured,      rownames = "uniprotid") %>%
    pivot_longer(-uniprotid, names_to = "sample", values_to = "log2") %>% mutate(stage = "before"),
  as_tibble(measured.norm, rownames = "uniprotid") %>%
    pivot_longer(-uniprotid, names_to = "sample", values_to = "log2") %>% mutate(stage = "after")) %>%
  left_join(sample_sheet, by = c("sample" = "sample_id")) %>%
  mutate(stage = factor(stage, c("before", "after")),
         sample = factor(sample, levels = sample_sheet$sample_id)) %>%
  ggplot(aes(sample, log2, fill = construct)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~ stage, ncol = 1) +
  cowplot::theme_minimal_grid() + cowplot::panel_border() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "", y = "log2 LFQ")
save_fig(p_normbox, "QC_normalization_boxplots", width = 9, height = 6)
p_normbox
# Sample medians are aligned after centering and no channel is an intensity outlier, so no sample is dropped and no further normalization is applied.

# Missing values -------------------------------------------------------------
# The question is whether missingness is intensity-dependent (left-censored) and what fraction is acceptable to impute rather than filter out.

# What is acceptable, and the imputation decision ----------------------------
# Retention keeps the group-wise valid-value rule: a protein is kept if it
# is detected in all three replicates of at least one construct-by- stimulation group. In proximity labeling a protein present in a bait and absent in the control is exactly the kind of hit we must not discard
# so no strict global missingness cut is applied at retention. Instead, each
# retained protein's global missingness drives an explicit imputation  decision rather than a silent blanket fill:
# | Global missing | Imputation policy                   |
# |----------------|-------------------------------------|
# | ≤ 20%          | impute directly                     |
# | \> 20 to 30%   | impute only after checking (ask)    |
# | \> 30 to 40%   | imputation discouraged, seek advice |
# | \> 40 to 50%   | imputation strongly discouraged     |
# | \> 50%         | do not impute                       |
# These bands annotate every protein and feed three analysis tiers: a
# primary high-confidence tier (≤ 20% missing, most robust), a sensitivity, tier (≤ 50% missing, more proteins), and an exploratory tier (strong bait-specific on/off proteins kept by the group-wise rule regardless of global missingness). The message printed below is the actual calculation on this dataset.

# chunk: missing-acceptable
missing_by_protein <- tibble(
  uniprotid = rownames(measured.norm),
  gene = prot.annotation$gene[match(rownames(measured.norm), prot.annotation$uniprotid)],
  n_missing = rowSums(is.na(measured.norm)),
  pct_missing = 100 * rowSums(is.na(measured.norm)) / ncol(measured.norm)) %>%
  mutate(impute_policy = cut(pct_missing, c(-Inf, 20, 30, 40, 50, Inf),
           labels = c("impute", "ask", "discouraged", "strongly_discouraged", "no_imputation")),
         missing_tier = cut(pct_missing, c(-Inf, 20, 50, Inf),
           labels = c("high_confidence_<=20", "sensitivity_<=50", "exploratory_>50")))
save_table(missing_by_protein, "protein_missingness")

# proteins retained at a range of ceilings
tibble(max_missing_pct = c(50, 40, 30, 20, 10, 0)) %>%
  mutate(n_proteins = map_int(max_missing_pct, ~ sum(missing_by_protein$pct_missing <= .x))) %>%
  dt_table()

# the calculation, printed per band, so the imputation decision is explicit
band_counts <- missing_by_protein %>% count(impute_policy, .drop = FALSE)
message("Missing-value imputation decision (", nrow(measured.norm), " retained proteins):")
for (i in seq_len(nrow(band_counts)))
  message(sprintf("  %-22s %4d proteins", as.character(band_counts$impute_policy[i]), band_counts$n[i]))
message("  Policy: <=20% impute; 20-30% ask; 30-40% discouraged; 40-50% strongly discouraged; >50% none.")
message("  ", sum(missing_by_protein$pct_missing > 50),
        " proteins exceed 50% missing -> flagged 'no_imputation' (kept, excluded from high-confidence).")

p_missing <- missing_by_protein %>%
  ggplot(aes(pct_missing, fill = impute_policy)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = c(20, 50), linetype = 2) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  cowplot::theme_minimal_grid() + cowplot::panel_border() +
  labs(x = "Percent of samples missing", y = "proteins", fill = "policy")
save_fig(p_missing, "QC_missingness_policy", width = 8, height = 5)
p_missing


# Is the missingness left-censored? ------------------------------------------

# chunk: missing-leftcensored
tibble(n_missing = rowSums(is.na(measured.norm)),
       mean_log2 = rowMeans(measured.norm, na.rm = TRUE)) %>%
  filter(n_missing > 0) %>%
  ggplot(aes(n_missing, mean_log2)) +
  geom_point(alpha = 0.3) + stat_smooth(method = "lm") + ggpubr::stat_cor() +
  cowplot::theme_minimal_grid() + cowplot::panel_border() +
  labs(x = "missing samples [n]", y = "mean log2 where measured")

# Missingness across the design ----------------------------------------------

# chunk: missing-by-group
map_dfr(names(by_group), ~ tibble(group = .x,
        pct_missing = 100 * mean(is.na(measured.norm[, by_group[[.x]]])))) %>%
  separate(group, c("construct", "stimulation"), sep = "_", remove = FALSE) %>%
  ggplot(aes(group, pct_missing, fill = construct)) +
  geom_col() +
  cowplot::theme_minimal_grid() + cowplot::panel_border() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "missing [%]")
# Mean intensity falls as the number of missing samples rises, so the  missingness is left-censored: proteins drop out near the detection limit,
# not at random. That makes a downshifted-normal imputation appropriate. The missingness is also higher in the control, which is expected: bait-
# specific proteins are genuinely absent there. Imputing those control gaps with low values is what lets a bait-only protein receive a finite, large fold change rather than being dropped.

# Detection per group --------------------------------------------------------
# For every protein the number of replicates it was measured in is recorded, per construct-by-stimulation group and per construct. This is what lets a  hit be reported as "detected 3/3 in IRE1α, 0/3 in control" alongside its statistics, and it feeds the imputation-dependence flag used in candidate calling.

# chunk: detection
det_group <- sapply(by_group, function(cols)
  rowSums(!is.na(measured.norm[, cols, drop = FALSE])))
det_group <- as_tibble(det_group, rownames = "uniprotid")

det_construct <- sapply(by_construct, function(cols)
  rowSums(!is.na(measured.norm[, cols, drop = FALSE])))
colnames(det_construct) <- paste0("det_", colnames(det_construct))
det_construct <- as_tibble(det_construct, rownames = "uniprotid")

detection <- det_group %>%
  rowwise() %>%
  mutate(detection_pattern = sprintf(
    "IRE1a %d/3,%d/3 | IRE1b %d/3,%d/3 | con %d/3,%d/3",
    IRE1a_DMSO, IRE1a_TM, IRE1b_DMSO, IRE1b_TM, con_DMSO, con_TM)) %>%
  ungroup() %>%
  left_join(det_construct, by = "uniprotid")
save_table(detection, "protein_detection_counts")
detection %>% select(uniprotid, detection_pattern, det_IRE1a, det_IRE1b, det_con) %>% dt_table()

# Imputation -----------------------------------------------------------------
# Missing values are imputed per sample from a normal shifted down by params_local$impute_downshift SD with width params_local$impute_width SD, a standard downshifted-normal imputation. Measured and imputed values are kept separate matrices so both are testable.

# chunk: imputation
impute_downshift_normal <- function(x, downshift, width) {
  ok <- !is.na(x); if (!any(ok)) return(x)
  mu <- mean(x[ok]); sdv <- sd(x[ok]); if (is.na(sdv) || sdv == 0) sdv <- 0.25
  n <- sum(!ok)
  if (n) x[!ok] <- rnorm(n, mu - downshift * sdv, width * sdv)
  x
}
set.seed(20230710)
imputed <- apply(measured.norm, 2, impute_downshift_normal,
                 downshift = params_local$impute_downshift,
                 width = params_local$impute_width)
rownames(imputed) <- rownames(measured.norm)
was_imputed <- is.na(measured.norm)
message("Imputed cells: ", sum(was_imputed), " of ", length(was_imputed),
        " (", round(100 * mean(was_imputed), 1), "%).")

# Principal component analysis -----------------------------------------------
# PCA on the imputed matrix, proteins centred not scaled.

# chunk: pca
pca <- prcomp(t(imputed), center = TRUE, scale. = FALSE)
ve  <- round(100 * summary(pca)$importance[2, 1:2], 1)
p_pca <- as_tibble(pca$x[, 1:2], rownames = "sample_id") %>%
  left_join(sample_sheet, by = "sample_id") %>%
  ggplot(aes(PC1, PC2, colour = construct, shape = stimulation)) +
  geom_point(size = 3.5) +
  cowplot::theme_minimal_grid() + cowplot::panel_border() +
  labs(x = paste0("PC1, ", ve[1], "%"), y = paste0("PC2, ", ve[2], "%"))
save_fig(p_pca, "QC_PCA", width = 7, height = 5)
p_pca

# Differential proximity analysis --------------------------------------------
# Primary evidence comes from two parallel streams, neither subordinate to the other: measured quantitative enrichment and reproducible qualitative bait-versus-control detection. The quantitative stream is the factorial `~ construct * stimulation` model fitted on the **measured** matrix, with the within-stimulation bait-versus-control contrasts as the hits; a protein qualifies only when the comparison has at least two observed values in each group, so no quantitative call rests on a 1-versus-2 comparison. The qualitative stream captures proteins reproducibly detected in a bait (3/3  replicates in at least one stimulation) and absent or rare in the matched control (≤ 1/3), which a measured test cannot score because the control side is essentially empty.

# The same factorial model fitted on the **imputed** matrix is retained as a  supporting sensitivity analysis. Proteins significant only after  imputation are labelled `imputation_supported` and do not carry the same evidential weight as measured or on/off candidates, so no biological conclusion rests on imputed values alone. A pooled model merging DMSO and TM reproduces the manuscript's approach, and a one-way ANOVA across constructs is a global diagnostic. All limma fits use moderated statistics (`trend`, `robust`) with Benjamini-Hochberg FDR; the moderated F is the ANOVA and each moderated contrast is a two-sample test.

# One-way ANOVA across constructs (diagnostic) -------------------------------

# chunk: anova
ss <- sample_sheet[match(colnames(imputed), sample_sheet$sample_id), ]

grp <- factor(ss$group_merged)
design_anova <- model.matrix(~ grp)                 # intercept + 2 dummies
fit_anova <- eBayes(lmFit(imputed, design_anova), trend = TRUE, robust = TRUE)
anova_tab <- topTable(fit_anova, coef = 2:3, number = Inf, sort.by = "none") %>%
  as_tibble(rownames = "uniprotid") %>%
  rename(FDR = adj.P.Val, P = P.Value) %>%
  left_join(select(prot.annotation, uniprotid, gene, protein), by = "uniprotid") %>%
  arrange(FDR)
save_table(anova_tab, "oneway_anova_merged")

message("One-way ANOVA (merged): ", sum(anova_tab$FDR < params_local$fdr_cutoff),
        " proteins at FDR < ", params_local$fdr_cutoff)
anova_tab %>%
  filter(FDR < params_local$fdr_cutoff) %>%
  select(gene, uniprotid, protein, F, P, FDR) %>%
  dt_table()

# Pairwise two-sample tests --------------------------------------------------

# chunk: pairwise
fit_pairs <- function(mat, grp, contrast_defs) {
  grp <- factor(grp)
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  cm <- makeContrasts(contrasts = contrast_defs, levels = design)
  colnames(cm) <- names(contrast_defs)
  eBayes(contrasts.fit(lmFit(mat, design), cm), trend = TRUE, robust = TRUE)
}
tidy_pairs <- function(fit, tag) {
  map_dfr(colnames(fit$contrasts), function(cn)
    topTable(fit, coef = cn, number = Inf, sort.by = "none", confint = TRUE) %>%
      as_tibble(rownames = "uniprotid") %>%
      mutate(contrast = cn, analysis = tag) %>%
      rename(FDR = adj.P.Val, P = P.Value)) %>%
    left_join(select(prot.annotation, uniprotid, gene, protein), by = "uniprotid")
}

# every pair among the three merged groups
merged_pairs <- c(IRE1a_vs_con   = "IRE1a - con",
                  IRE1b_vs_con   = "IRE1b - con",
                  IRE1a_vs_IRE1b = "IRE1a - IRE1b")
# factorial contrasts (kept as supplement)
full_pairs <- c(
  IRE1a_vs_con_DMSO = "IRE1a_DMSO - con_DMSO", IRE1a_vs_con_TM = "IRE1a_TM - con_TM",
  IRE1b_vs_con_DMSO = "IRE1b_DMSO - con_DMSO", IRE1b_vs_con_TM = "IRE1b_TM - con_TM",
  IRE1a_vs_IRE1b_DMSO = "IRE1a_DMSO - IRE1b_DMSO", IRE1a_vs_IRE1b_TM = "IRE1a_TM - IRE1b_TM",
  TM_vs_DMSO_IRE1a = "IRE1a_TM - IRE1a_DMSO", TM_vs_DMSO_IRE1b = "IRE1b_TM - IRE1b_DMSO",
  TM_vs_DMSO_con   = "con_TM - con_DMSO",
  interaction_IRE1a = "(IRE1a_TM - IRE1a_DMSO) - (con_TM - con_DMSO)",
  interaction_IRE1b = "(IRE1b_TM - IRE1b_DMSO) - (con_TM - con_DMSO)")

merged_imp  <- tidy_pairs(fit_pairs(imputed,       ss$group_merged, merged_pairs), "imputed")
merged_meas <- tidy_pairs(fit_pairs(measured.norm, ss$group_merged, merged_pairs), "measured_sensitivity")
full_imp    <- tidy_pairs(fit_pairs(imputed,       ss$group_full,   full_pairs),   "imputed")
full_meas   <- tidy_pairs(fit_pairs(measured.norm, ss$group_full,   full_pairs),   "measured_sensitivity")

save_table(bind_rows(merged_imp, merged_meas, full_imp, full_meas), "pairwise_all_contrasts")

bind_rows(merged_imp, merged_meas, full_imp, full_meas) %>%
  group_by(analysis, contrast) %>%
  summarise(n_tested = sum(!is.na(FDR)),
            n_enriched = sum(FDR < params_local$fdr_cutoff & logFC > params_local$lfc_cutoff, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(analysis, contrast) %>%
  dt_table()

# Volcano plots --------------------------------------------------------------
# Primary volcanoes are the four within-stimulation bait-versus-control contrasts on the **measured** matrix. The imputed volcanoes follow as supporting. Baits are always labelled.

# chunk: volcano
bvc <- c("IRE1a_vs_con_DMSO", "IRE1a_vs_con_TM",
         "IRE1b_vs_con_DMSO", "IRE1b_vs_con_TM")
make_volcano <- function(res, cn, tag) {
  d <- res %>% filter(contrast == cn) %>%
    mutate(sig = FDR < params_local$fdr_cutoff & logFC > params_local$lfc_cutoff,
           is_bait = gene %in% params_local$baits,
           flag_biotin = gene %in% params_local$biotin_carboxylases,
           label = label_for(uniprotid))
  if (nrow(d) == 0 || all(is.na(d$FDR))) return(invisible())
  p <- ggplot(d, aes(logFC, -log10(FDR))) +
    geom_point(aes(colour = sig, shape = flag_biotin), alpha = 0.7, size = 1.6) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#b2182b")) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey50") +
    geom_hline(yintercept = -log10(params_local$fdr_cutoff), linetype = 2, colour = "grey50") +
    ggrepel::geom_text_repel(
      data = bind_rows(slice_max(filter(d, sig), logFC, n = 15), filter(d, is_bait)),
      aes(label = label), size = 3, max.overlaps = 25) +
    cowplot::theme_minimal_grid() + cowplot::panel_border() +
    labs(title = paste0(gsub("_", " ", cn), "  (", tag, ")"),
         x = "log2 fold change", y = "-log10 FDR", colour = "enriched", shape = "endog. biotin")
  save_fig(p, paste0("volcano_", tag, "_", cn), width = 8, height = 6)
  print(p)
}
# primary: measured
for (cn in bvc) make_volcano(full_meas, cn, "measured")

# Supporting volcanoes (imputed) {#supporting-volcanoes} ---------------------

# chunk: volcano-imputed
for (cn in c(bvc, "interaction_IRE1a", "interaction_IRE1b"))
  make_volcano(full_imp, cn, "imputed")

# Candidate classification and heatmap ---------------------------------------
# A candidate is a protein either enriched over the matched control in the imputed factorial model (FDR \< params_local$fdr_cutoff, log2FC \> params_local$lfc_cutoff, in either stimulation) or showing a clean on/off  detection pattern. Each candidate is placed in an evidence tier, with the first two tiers as parallel primary evidence and imputation as supporting:

# - **high_confidence**: significant in the measured-data analysis, with at least params_local$min_obs_per_group observed values in each group of the contrast, so it never rests on a sparse comparison; -
# **qualitative_on_off**: reproducibly detected in a bait (3/3 replicates in at least one stimulation) and absent or rare in the matched control (≤ 1/3), a strong biological pattern that a measured limma test cannot score because the control side is essentially empty. This tier is assigned before the high-missingness downgrade, so bait-specific on/off proteins are preserved; - **imputation_supported**: significant only after imputation, supporting evidence, explicitly labelled; - **exploratory**: a candidate with high missingness (\> 50%) and neither measured nor on/off support.

# Endogenous-biotin background is flagged, and protein groups with more than one accession are flagged as ambiguous rather than presented as a single gene product. Classification uses the statistical contrasts and reproducible detection, not a detection Venn diagram.

# chunk: candidates
hit_flag <- function(lfc, fdr)
  !is.na(fdr) & fdr < params_local$fdr_cutoff & lfc > params_local$lfc_cutoff
get_i <- function(cn) full_imp  %>% filter(contrast == cn) %>% select(uniprotid, logFC, FDR)
get_m <- function(cn) full_meas %>% filter(contrast == cn) %>% select(uniprotid, logFC, FDR)
mo <- params_local$min_obs_per_group

candidates <- prot.annotation %>%
  filter(uniprotid %in% rownames(imputed)) %>%
  select(uniprotid, gene, protein, n_accessions, protein_group_ambiguous) %>%
  # imputed within-stimulation bait-vs-control (sensitive; defines the candidate universe)
  left_join(get_i("IRE1a_vs_con_DMSO") %>% rename(a_DMSO_lfc = logFC, a_DMSO_fdr = FDR), by = "uniprotid") %>%
  left_join(get_i("IRE1a_vs_con_TM")   %>% rename(a_TM_lfc = logFC,   a_TM_fdr = FDR),   by = "uniprotid") %>%
  left_join(get_i("IRE1b_vs_con_DMSO") %>% rename(b_DMSO_lfc = logFC, b_DMSO_fdr = FDR), by = "uniprotid") %>%
  left_join(get_i("IRE1b_vs_con_TM")   %>% rename(b_TM_lfc = logFC,   b_TM_fdr = FDR),   by = "uniprotid") %>%
  # measured within-stimulation bait-vs-control (primary quantitative evidence)
  left_join(get_m("IRE1a_vs_con_DMSO") %>% rename(ma_DMSO_lfc = logFC, ma_DMSO_fdr = FDR), by = "uniprotid") %>%
  left_join(get_m("IRE1a_vs_con_TM")   %>% rename(ma_TM_lfc = logFC,   ma_TM_fdr = FDR),   by = "uniprotid") %>%
  left_join(get_m("IRE1b_vs_con_DMSO") %>% rename(mb_DMSO_lfc = logFC, mb_DMSO_fdr = FDR), by = "uniprotid") %>%
  left_join(get_m("IRE1b_vs_con_TM")   %>% rename(mb_TM_lfc = logFC,   mb_TM_fdr = FDR),   by = "uniprotid") %>%
  left_join(detection %>% select(uniprotid, detection_pattern,
                                 IRE1a_DMSO, IRE1a_TM, IRE1b_DMSO, IRE1b_TM, con_DMSO, con_TM,
                                 det_IRE1a, det_IRE1b, det_con), by = "uniprotid") %>%
  left_join(missing_by_protein %>% select(uniprotid, pct_missing, impute_policy, missing_tier),
            by = "uniprotid") %>%
  mutate(
    # on/off qualitative pattern: bait 3/3 in a stimulation AND matched control <=1 there
    a_on_off = (IRE1a_DMSO == 3 & con_DMSO <= 1) | (IRE1a_TM == 3 & con_TM <= 1),
    b_on_off = (IRE1b_DMSO == 3 & con_DMSO <= 1) | (IRE1b_TM == 3 & con_TM <= 1),
    # candidate universe = enriched in imputed limma OR clean on/off detection
    cand_IRE1a = hit_flag(a_DMSO_lfc, a_DMSO_fdr) | hit_flag(a_TM_lfc, a_TM_fdr) | (a_on_off %in% TRUE),
    cand_IRE1b = hit_flag(b_DMSO_lfc, b_DMSO_fdr) | hit_flag(b_TM_lfc, b_TM_fdr) | (b_on_off %in% TRUE),
    class = case_when(cand_IRE1a & cand_IRE1b ~ "shared",
                      cand_IRE1a ~ "IRE1a_specific",
                      cand_IRE1b ~ "IRE1b_specific",
                      TRUE ~ "not_enriched"),
    best_lfc = pmax(a_DMSO_lfc, a_TM_lfc, b_DMSO_lfc, b_TM_lfc, na.rm = TRUE),
    flag_biotin = gene %in% params_local$biotin_carboxylases,
    qualitative_on_off = (cand_IRE1a & a_on_off %in% TRUE) | (cand_IRE1b & b_on_off %in% TRUE),
    # measured quantitative support: significant AND >= min_obs observed in BOTH groups of the contrast
    ok_a_DMSO = hit_flag(ma_DMSO_lfc, ma_DMSO_fdr) & IRE1a_DMSO >= mo & con_DMSO >= mo,
    ok_a_TM   = hit_flag(ma_TM_lfc,   ma_TM_fdr)   & IRE1a_TM   >= mo & con_TM   >= mo,
    ok_b_DMSO = hit_flag(mb_DMSO_lfc, mb_DMSO_fdr) & IRE1b_DMSO >= mo & con_DMSO >= mo,
    ok_b_TM   = hit_flag(mb_TM_lfc,   mb_TM_fdr)   & IRE1b_TM   >= mo & con_TM   >= mo,
    measured_supported = ok_a_DMSO %in% TRUE | ok_a_TM %in% TRUE |
                         ok_b_DMSO %in% TRUE | ok_b_TM %in% TRUE,
    measured_min_FDR = suppressWarnings(pmin(
      ifelse(ok_a_DMSO %in% TRUE, ma_DMSO_fdr, NA_real_),
      ifelse(ok_a_TM   %in% TRUE, ma_TM_fdr,   NA_real_),
      ifelse(ok_b_DMSO %in% TRUE, mb_DMSO_fdr, NA_real_),
      ifelse(ok_b_TM   %in% TRUE, mb_TM_fdr,   NA_real_), na.rm = TRUE)),
    measured_min_FDR = ifelse(is.infinite(measured_min_FDR), NA_real_, measured_min_FDR),
    imputation_dependent = class != "not_enriched" &
                           !(measured_supported %in% TRUE) & !(qualitative_on_off %in% TRUE),
    # tier order: measured first, then on/off (protected from the missingness downgrade),
    # then exploratory, then imputation-only
    evidence_tier = case_when(
      class == "not_enriched"           ~ "not_enriched",
      flag_biotin                       ~ "flagged_background",
      measured_supported %in% TRUE      ~ "high_confidence",
      qualitative_on_off %in% TRUE      ~ "qualitative_on_off",
      missing_tier == "exploratory_>50" ~ "exploratory",
      TRUE                              ~ "imputation_supported"))
save_table(candidates, "candidate_confidence_table")
candidates %>% filter(class != "not_enriched") %>% count(class, evidence_tier) %>% dt_table()

# Candidates by evidence tier ------------------------------------------------
# Two parallel primary streams, then the supporting imputation-only set.

# chunk: tier-tables
# 1. high_confidence: measured quantitative evidence (>= min_obs per group)
candidates %>%
  filter(evidence_tier == "high_confidence") %>%
  arrange(measured_min_FDR, desc(best_lfc)) %>%
  transmute(gene, uniprotid, protein, class, detection_pattern,
            measured_min_FDR, best_log2FC = best_lfc, pct_missing) %>%
  dt_table()

# chunk: tier-on/off
# 2. qualitative_on_off: reproducible bait detection, absent/rare in control
candidates %>%
  filter(evidence_tier == "qualitative_on_off") %>%
  arrange(desc(det_IRE1a + det_IRE1b), det_con) %>%
  transmute(gene, uniprotid, protein, class, detection_pattern, pct_missing,
            measured_supported) %>%
  dt_table()

# chunk: tier-imputation
# 3. imputation_supported: significant only after imputation (supporting evidence)
candidates %>%
  filter(evidence_tier == "imputation_supported") %>%
  arrange(desc(best_lfc)) %>%
  transmute(gene, uniprotid, protein, class, detection_pattern,
            best_log2FC = best_lfc, pct_missing) %>%
  dt_table()

# Supporting: pooled reproduction (imputed, merged DMSO + TM) ----------------

# chunk: pooled-reproduction
bind_rows(
  merged_imp %>% filter(contrast == "IRE1a_vs_con"),
  merged_imp %>% filter(contrast == "IRE1b_vs_con")) %>%
  filter(FDR < params_local$fdr_cutoff, logFC > params_local$lfc_cutoff) %>%
  arrange(contrast, desc(logFC)) %>%
  select(contrast, gene, uniprotid, protein, logFC, CI.L, CI.R, P, FDR) %>%
  dt_table()

# chunk: heatmap
top_ids <- candidates %>%
  filter(class != "not_enriched") %>%
  group_by(class) %>% slice_max(best_lfc, n = 20) %>% ungroup() %>%
  pull(uniprotid)
# always include the baits
bait_ids <- prot.annotation$uniprotid[prot.annotation$gene %in% params_local$baits]
top_ids <- unique(c(bait_ids, top_ids))
top_ids <- intersect(top_ids, rownames(imputed))

hm <- imputed[top_ids, sample_sheet$sample_id]
hm <- (hm - rowMeans(hm)) / rowSds(hm)
amb <- prot.annotation$protein_group_ambiguous[match(rownames(hm), prot.annotation$uniprotid)]
rownames(hm) <- ifelse(amb %in% TRUE,                      # * marks ambiguous groups
                       paste0(label_for(rownames(hm)), " *"),
                       label_for(rownames(hm)))

ht <- ComplexHeatmap::Heatmap(
  hm, name = "z-score",
  col = circlize::colorRamp2(c(-2.5, 0, 2.5), c("#2166ac", "white", "#b2182b")),
  column_split = factor(sample_sheet$construct, c("IRE1a", "IRE1b", "con")),
  top_annotation = HeatmapAnnotation(
    construct = sample_sheet$construct, stimulation = sample_sheet$stimulation,
    col = list(construct = c(IRE1a = "#E41A1C", IRE1b = "#377EB8", con = "#4DAF4A"),
               stimulation = c(DMSO = "grey80", TM = "grey30"))),
  row_names_gp = gpar(fontsize = 7), show_column_names = FALSE)
save_heatmap(ht, "heatmap_top_candidates", width = 9, height = 10)
ht
# ERN1 (IRE1α) and ERN2 (IRE1β) are forced into the heatmap and are the strongest rows in their own construct blocks, which is the internal positive control for the whole experiment. Rows marked with an asteriskare ambiguous protein groups (more than one accession) and should not be read as a single uniquely identified gene product.

# Additional quality control -------------------------------------------------
# MA plots, replicate correlation, hierarchical clustering and coefficient- of- variation diagnostics, on the imputed matrix unless noted.

# MA plots -------------------------------------------------------------------
# For each bait-versus-control comparison, M (log2 fold change) against A (mean log2 abundance). A funnel that widens at low A is the expected mean- variance pattern; systematic curvature would indicate a normalization problem.

# chunk: ma-plots
ma_one <- function(bait) {
  a_cols <- by_construct[[bait]]; c_cols <- by_construct[["con"]]
  tibble(uniprotid = rownames(imputed),
         A = 0.5 * (rowMeans(imputed[, a_cols]) + rowMeans(imputed[, c_cols])),
         M = rowMeans(imputed[, a_cols]) - rowMeans(imputed[, c_cols]),
         comparison = paste0(bait, " vs con"))
}
ma_df <- bind_rows(ma_one("IRE1a"), ma_one("IRE1b")) %>%
  left_join(bind_rows(
    merged_imp %>% filter(contrast == "IRE1a_vs_con") %>% mutate(comparison = "IRE1a vs con"),
    merged_imp %>% filter(contrast == "IRE1b_vs_con") %>% mutate(comparison = "IRE1b vs con")) %>%
    select(uniprotid, comparison, FDR), by = c("uniprotid", "comparison")) %>%
  mutate(sig = !is.na(FDR) & FDR < params_local$fdr_cutoff & M > params_local$lfc_cutoff)

p_ma <- ggplot(ma_df, aes(A, M, colour = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#b2182b")) +
  geom_hline(yintercept = 0, linetype = 1, colour = "grey50") +
  geom_hline(yintercept = c(-1, 1), linetype = 2, colour = "grey50") +
  facet_wrap(~ comparison) +
  cowplot::theme_minimal_grid() + cowplot::panel_border() +
  labs(x = "A (mean log2 abundance)", y = "M (log2 fold change)")
save_fig(p_ma, "QC_MA_plots", width = 10, height = 5)
p_ma

# Replicate correlation and hierarchical clustering --------------------------

# chunk: corr-heatmap
cor_m <- cor(imputed, use = "pairwise.complete.obs")
save_table(as_tibble(cor_m, rownames = "sample"), "sample_correlation_matrix")

ann <- HeatmapAnnotation(
  construct = sample_sheet$construct, stimulation = sample_sheet$stimulation,
  col = list(construct = c(IRE1a = "#E41A1C", IRE1b = "#377EB8", con = "#4DAF4A"),
             stimulation = c(DMSO = "grey80", TM = "grey30")))
ht_cor <- ComplexHeatmap::Heatmap(
  cor_m, name = "Pearson r",
  col = circlize::colorRamp2(c(min(cor_m), 1), c("white", "#08519c")),
  top_annotation = ann,
  clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7))
save_heatmap(ht_cor, "QC_replicate_correlation", width = 9, height = 8)
ht_cor

# chunk: hclust
hc <- hclust(as.dist(1 - cor_m), method = "average")
draw_dend <- function()
  plot(hc, main = "Sample hierarchical clustering (1 - Pearson r)", xlab = "", sub = "")
save_base_plot(draw_dend, "QC_hierarchical_clustering", width = 9, height = 6)
draw_dend()

# Replicate CV ---------------------------------------------------------------
# Within-group coefficient of variation on the linear scale, from measured values only (imputed values are excluded so the CV reflects real reproducibility), for proteins with at least two measured replicates in the group.

# chunk: cv
lin <- 2^measured.norm
cv_long <- map_dfr(names(by_group), function(g) {
  m <- lin[, by_group[[g]], drop = FALSE]
  tibble(uniprotid = rownames(lin), group = g,
         n_valid = rowSums(!is.na(m)),
         cv = 100 * matrixStats::rowSds(m, na.rm = TRUE) / rowMeans(m, na.rm = TRUE)) %>%
    filter(n_valid >= 2)
})
cv_summary <- cv_long %>% group_by(group) %>%
  summarise(n = n(), median_CV = median(cv, na.rm = TRUE),
            mean_CV = mean(cv, na.rm = TRUE), .groups = "drop")
save_table(cv_summary, "replicate_CV_summary")
cv_summary %>% dt_table()

p_cv <- cv_long %>%
  separate(group, c("construct", "stimulation"), sep = "_", remove = FALSE) %>%
  ggplot(aes(group, cv, fill = construct)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white") +
  coord_cartesian(ylim = c(0, quantile(cv_long$cv, 0.98, na.rm = TRUE))) +
  cowplot::theme_minimal_grid() + cowplot::panel_border() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "CV [%] (within-group, linear scale)")
save_fig(p_cv, "QC_replicate_CV", width = 8, height = 5)
p_cv

# GO and pathway enrichment --------------------------------------------------
# Enrichment is split by the kind of database, because the two need different handling. Pathways (KEGG and Reactome) are run with `gprofiler2::gost()`, which gives specific pathway terms and a proper enrichment p-value against a background; those are the pathway bar and dot plots. GO (biological process, molecular function, cellular component) is instead computed from direct GO annotations in `org.Hs.eg.db`, because  :Profiler and clusterProfiler propagate each gene's annotations up the  hierarchy and the highest-count terms then become abstract parents ("localization", "transport"). Direct annotations keep the specific terms, matching a DAVID-style GO summary, and are covered in the GO section further down. For the pathway call, the query gene set defaults to the union of proteins enriched over control in either bait; set `GO_QUERY_CLASS` to a single class to restrict it. The `GO_BACKGROUND` switch chooses the annotated- genome background (default, populates both databases) or a custom background of the detected proteins (stricter, sparser).

# chunk: go-run
gostres <- NULL
if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  message("Install gprofiler2: install.packages('gprofiler2')")
} else {
  library(gprofiler2)

  GO_QUERY_CLASS <- "any_enriched"     # or IRE1a_specific / IRE1b_specific / shared
  GO_BACKGROUND  <- "annotated"        # "annotated" (genome, default) or "custom" (tested proteins)

  query_genes <- if (GO_QUERY_CLASS == "any_enriched")
    candidates$gene[candidates$class != "not_enriched"] else
    candidates$gene[candidates$class == GO_QUERY_CLASS]
  query_genes <- unique(na.omit(query_genes))
  background  <- unique(na.omit(prot.annotation$gene[
    match(rownames(imputed), prot.annotation$uniprotid)]))
  message("g:Profiler query: ", length(query_genes), " genes; background mode: ", GO_BACKGROUND)

  gost_args <- list(query = query_genes, organism = "hsapiens",
                    sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),
                    correction_method = "g_SCS", significant = TRUE, evcodes = TRUE)
  if (GO_BACKGROUND == "custom")
    gost_args <- c(gost_args, list(custom_bg = background, domain_scope = "custom"))

  gost_raw <- do.call(gprofiler2::gost, gost_args)

  if (is.null(gost_raw) || is.null(gost_raw$result) || nrow(gost_raw$result) == 0) {
    message("g:Profiler returned no significant terms for this query/background.")
  } else {
    gostres <- as_tibble(gost_raw$result) %>%
      mutate(genes = vapply(intersection, function(x) paste(x, collapse = ", "),
                            character(1), USE.NAMES = FALSE)) %>%
      select(source, term_id, term_name, p_value, term_size,
             query_size, intersection_size, genes) %>%
      mutate(pct = 100 * intersection_size / query_size)
    save_table(gostres, "gprofiler_enrichment_all")
  }
}

# What came back per source --------------------------------------------------
# A checkpoint before plotting: how many significant terms each source returned. If a source is missing here, its panel below will say so rather than draw an empty box.

# chunk: go-source-counts
if (!is.null(gostres)) {
  gostres %>% count(source, name = "n_terms") %>% arrange(desc(n_terms)) %>% dt_table()
}

# chunk: go-palettes
pal_go     <- c(BP = "#1F4E9B", CC = "#E8000B", MF = "#33A02C")
pal_go_dk  <- c(`Biological Process` = "#8B1A1A",
                `Molecular Function` = "#1B6B1B",
                `Cellular component` = "#1B2E6B")
pal_source <- c(KEGG = "#2E4A52", REAC = "#E69F00")
pal_db     <- c(KEGG = "black",   REACTOME = "#8B0000")
msg_empty  <- function(what) message("No ", what, " terms to plot for this query/background.")

# Significant pathways (KEGG + Reactome) -------------------------------------
# Image-1 style: bars are `-log10(p)`, coloured by database.

# chunk: fig-pathways-bar
if (!is.null(gostres)) {
  d <- gostres %>% filter(source %in% c("KEGG", "REAC"))
  if (nrow(d) == 0) msg_empty("KEGG/Reactome") else {
    p_path_bar <- d %>% slice_min(p_value, n = 30) %>%
      mutate(term_name = fct_reorder(term_name, -log10(p_value))) %>%
      ggplot(aes(-log10(p_value), term_name, fill = source)) +
      geom_col() + scale_fill_manual(values = pal_source) +
      labs(x = expression(-log[10](p)), y = NULL, fill = "source") +
      theme_classic(base_size = 13) +
      theme(legend.position = "top", axis.text.y = element_text(size = 11))
    save_fig(p_path_bar, "Significant_Pathways", width = 11, height = 7)
    print(p_path_bar)
  }
}

# Pathway enrichment dot plot ------------------------------------------------
# Image-2 style: x is the number of detected proteins, dot size is the  adjusted p, colour is the database.

# chunk: fig-pathways-dot
if (!is.null(gostres)) {
  d <- gostres %>% filter(source %in% c("KEGG", "REAC"))
  if (nrow(d) == 0) msg_empty("KEGG/Reactome") else {
    p_path_dot <- d %>% slice_min(p_value, n = 28) %>%
      mutate(Database = recode(source, KEGG = "KEGG", REAC = "REACTOME"),
             term_name = fct_reorder(term_name, intersection_size)) %>%
      ggplot(aes(intersection_size, term_name)) +
      geom_point(aes(size = p_value, colour = Database)) +
      scale_size("p.adj", trans = "reverse", range = c(1.5, 7)) +
      scale_colour_manual(values = pal_db) +
      labs(x = "# of detected proteins", y = "Pathway") +
      theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())
    save_fig(p_path_dot, "Enrichment_analysis", width = 10, height = 9)
    print(p_path_dot)
  }
}

# GO enrichment (direct annotations) -----------------------------------------
# The GO figures below use *direct* GO annotations from `org.Hs.eg.db`, not the propagated (ancestor) annotations that g:Profiler and clusterProfiler return by default. Propagated annotations push abstract parent terms ("localization", "transport", "cellular process") to the top by count and bury the specific biology. Direct annotations assign each gene only its most specific terms, which is what DAVID does and what gives interpretable terms such as "regulation of cellular response to heat" or "Hsp90 protein  binding". Each term also carries a hypergeometric enrichment p-value and  BH-adjusted FDR against the tested-protein background, so the figure is descriptive and testable.

# chunk: go-direct
go_direct <- NULL
gene_go_direct <- tibble(gene = character())
have_go <- all(vapply(c("org.Hs.eg.db", "GO.db", "AnnotationDbi"),
                      requireNamespace, logical(1), quietly = TRUE))
if (!have_go) {
  message("Install org.Hs.eg.db, GO.db and AnnotationDbi for direct-annotation GO.")
} else {
  suppressPackageStartupMessages({library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi)})

  GO_QUERY_CLASS <- "any_enriched"   # or IRE1a_specific / IRE1b_specific / shared
  q_genes <- if (GO_QUERY_CLASS == "any_enriched")
    candidates$gene[candidates$class != "not_enriched"] else
    candidates$gene[candidates$class == GO_QUERY_CLASS]
  q_genes  <- unique(as.character(na.omit(q_genes)))
  bg_genes <- unique(as.character(na.omit(prot.annotation$gene[
    match(rownames(imputed), prot.annotation$uniprotid)])))

  # symbols -> Entrez, always returning a plain character vector
  sym2entrez <- function(s) {
    s <- unique(as.character(s)); s <- s[!is.na(s) & s != ""]
    if (!length(s)) return(character(0))
    ez <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = s, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
    unique(as.character(na.omit(unname(ez))))
  }
  q_ez  <- sym2entrez(q_genes)
  bg_ez <- unique(c(sym2entrez(bg_genes), q_ez))

  if (!length(q_ez)) {
    message("No query genes mapped to Entrez; skipping direct GO.")
  } else {
    # DIRECT annotations: the "GO" column (not "GOALL") is non-propagated
    get_direct <- function(ez) {
      ez <- unique(as.character(ez)); ez <- ez[!is.na(ez) & ez != ""]
      if (!length(ez)) return(tibble(ENTREZID = character(), GO = character(), ONTOLOGY = character()))
      suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db, keys = ez, keytype = "ENTREZID", columns = c("GO", "ONTOLOGY"))) %>%
        as_tibble() %>% filter(!is.na(GO)) %>% distinct(ENTREZID, GO, ONTOLOGY)
    }
    qa <- get_direct(q_ez); ba <- get_direct(bg_ez)

    ez2sym <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = q_ez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
    n_query <- length(q_ez); n_bg <- length(bg_ez)

    go_ids <- unique(as.character(qa$GO)); go_ids <- go_ids[!is.na(go_ids)]
    term_names <- suppressMessages(AnnotationDbi::select(
      GO.db, keys = go_ids, keytype = "GOID", columns = "TERM")) %>%
      as_tibble() %>% rename(GO = GOID, term = TERM)

    go_direct <- qa %>%
      group_by(GO, ONTOLOGY) %>%
      summarise(Count = n_distinct(ENTREZID),
                genes = paste(sort(unique(ez2sym[ENTREZID])), collapse = ", "),
                .groups = "drop") %>%
      left_join(ba %>% count(GO, name = "bg_count"), by = "GO") %>%
      left_join(term_names, by = "GO") %>%
      mutate(pct = 100 * Count / n_query,
             p = phyper(Count - 1, bg_count, n_bg - bg_count, n_query, lower.tail = FALSE)) %>%
      group_by(ONTOLOGY) %>% mutate(FDR = p.adjust(p, "BH")) %>% ungroup() %>%
      mutate(ontology = recode(ONTOLOGY, BP = "Biological Process",
                               MF = "Molecular Function", CC = "Cellular component")) %>%
      filter(!is.na(term))
    save_table(go_direct %>% select(ontology, GO, term, Count, pct, p, FDR, genes),
               "GO_direct_enrichment_all")

    # gene -> direct terms per ontology, for the significant-protein table
    gene_go_direct <- qa %>%
      left_join(term_names, by = "GO") %>%
      mutate(gene = ez2sym[ENTREZID]) %>% filter(!is.na(gene), !is.na(term)) %>%
      group_by(gene, ONTOLOGY) %>%
      summarise(terms = paste(head(unique(term), 6), collapse = "; "), .groups = "drop") %>%
      pivot_wider(names_from = ONTOLOGY, values_from = terms) %>%
      rename(GO_BP = any_of("BP"), GO_MF = any_of("MF"), GO_CC = any_of("CC"))

    message("Direct GO terms: ", nrow(go_direct), " across ", n_query, " query genes.")
  }
}

# GO, all three ontologies ---------------------------------------------------
# Image-3 style: counts per GO term, grouped and coloured by ontology, top terms per ontology, using direct annotations.

# chunk: fig-go-all
if (!is.null(go_direct)) {
  N_PER_ONT <- 8
  go_all <- go_direct %>%
    mutate(source = recode(ONTOLOGY, BP = "BP", CC = "CC", MF = "MF")) %>%
    group_by(source) %>% slice_max(Count, n = N_PER_ONT, with_ties = FALSE) %>% ungroup() %>%
    mutate(source = factor(source, levels = c("MF", "CC", "BP")),
           row = paste(source, term, sep = "|")) %>%
    arrange(source, Count) %>%
    mutate(row = factor(row, levels = row))
  p_go_all <- ggplot(go_all, aes(Count, row, fill = source)) +
    geom_col() +
    scale_fill_manual(values = pal_go, breaks = c("BP", "CC", "MF")) +
    scale_y_discrete(labels = function(x) sub("^[A-Z]+\\|", "", x)) +
    labs(x = "Count", y = "Goterm", fill = "source") +
    theme_classic(base_size = 13) +
    theme(legend.position = "top", axis.text.y = element_text(size = 11))
  save_fig(p_go_all, "GO_All", width = 11, height = 8)
  p_go_all
}

# GO summary panels and count table ------------------------------------------
# Image-4 style: one panel per ontology (top 10 by protein count) plus a table giving the count, the percent of query proteins, the FDR, and the genes in each term.

# chunk: fig-go-faceted
if (!is.null(go_direct)) {
  if (!requireNamespace("tidytext", quietly = TRUE))
    message("Install tidytext for the faceted GO panels: install.packages('tidytext')")
  go_facet <- go_direct %>%
    mutate(ontology = factor(ontology,
             c("Biological Process", "Molecular Function", "Cellular component"))) %>%
    group_by(ontology) %>% slice_max(Count, n = 10, with_ties = FALSE) %>% ungroup()
  p_go_facet <- go_facet %>%
    mutate(term = tidytext::reorder_within(term, Count, ontology)) %>%
    ggplot(aes(Count, term, fill = ontology)) +
    geom_col() + tidytext::scale_y_reordered() +
    scale_fill_manual(values = pal_go_dk, guide = "none") +
    facet_wrap(~ ontology, ncol = 1, scales = "free_y") +
    labs(x = "Protein count", y = NULL) +
    theme_bw(base_size = 12) + theme(strip.text = element_text(face = "bold"))
  save_fig(p_go_facet, "GoPlotInfo", width = 9, height = 9)
  p_go_facet
}

# chunk: tbl-go-info
if (!is.null(go_direct)) {
  go_info <- go_direct %>%
    group_by(ontology) %>% slice_max(Count, n = 10, with_ties = FALSE) %>% ungroup() %>%
    transmute(Ontology = ontology, Goterm = term, Count,
              `%` = sprintf("%.1f%%", pct), FDR = signif(FDR, 3), Genes = genes) %>%
    arrange(factor(Ontology, c("Biological Process","Molecular Function","Cellular component")),
            desc(Count))
  save_table(go_info, "GO_count_percent_table")
  go_info %>% dt_table()
}
# Enrichment is reported as proximity, not physical interaction: GO:CC terms in particular reflect where the labeled proteins localize (endoplasmic reticulum, membrane, cytosol) and are expected to be dominated by the bait compartment.

# Significant proteins with roles and GO terms -------------------------------
# Every protein significant over control in the primary (imputed) analysis, in any of the three merged contrasts, is collected here with its statistics and its functional annotation: the descriptive protein name serves as its role, and its direct GO terms (biological process, molecular function, cellular component) are attached from `org.Hs.eg.db`. This is the master per-protein table for the candidate list, written to the tables folder.

# Every enriched protein (factorial model, either stimulation) is collected here with its statistics, its role (the descriptive protein name), its direct GO terms from `org.Hs.eg.db`, and the evidence needed to judge it:the detection pattern (for example 3/3 in a bait, 0/3 in control), whether it is supported by the measured-data analysis, whether it is imputation-only, its evidence tier (high_confidence, qualitative_on_off, imputation_supported, or exploratory), and whether its protein group is ambiguous. This is the master per-protein table, sorted so measured-  supported and on/off proteins come first.

# chunk: sig-proteins-go
sig_proteins <- candidates %>%
  filter(class != "not_enriched") %>%
  left_join(gene_go_direct, by = "gene")
for (col in c("GO_BP", "GO_MF", "GO_CC"))
  if (!col %in% names(sig_proteins)) sig_proteins[[col]] <- NA_character_

sig_proteins <- sig_proteins %>%
  transmute(uniprotid, gene, role = protein, class, evidence_tier,
            detection_pattern, best_log2FC = best_lfc, pct_missing, impute_policy,
            measured_supported, measured_min_FDR, imputation_dependent,
            qualitative_on_off, protein_group_ambiguous, n_accessions,
            GO_BP, GO_MF, GO_CC) %>%
  arrange(factor(evidence_tier,
                 c("high_confidence","qualitative_on_off","imputation_supported",
                   "exploratory","flagged_background")),
          desc(best_log2FC))

save_table(sig_proteins, "significant_proteins_with_GO_roles")
message(nrow(sig_proteins), " enriched proteins | ",
        sum(sig_proteins$evidence_tier == "high_confidence", na.rm = TRUE), " high-confidence, ",
        sum(sig_proteins$evidence_tier == "qualitative_on_off", na.rm = TRUE), " on/off, ",
        sum(sig_proteins$evidence_tier == "imputation_supported", na.rm = TRUE), " imputation-only, ",
        sum(sig_proteins$protein_group_ambiguous, na.rm = TRUE), " ambiguous groups.")
sig_proteins %>% dt_table()

# Export ---------------------------------------------------------------------

# chunk: export
# standalone reproducibility tables (the sample sheet is an input, not re-exported)
save_table(prot.annotation %>% filter(uniprotid %in% rownames(imputed)), "protein_annotation")
save_table(as_tibble(measured.norm, rownames = "uniprotid"), "abundance_measured_wide")
save_table(as_tibble(imputed,       rownames = "uniprotid"), "abundance_imputed_wide")

abundance_long <- as_tibble(imputed, rownames = "uniprotid") %>%
  pivot_longer(-uniprotid, names_to = "sample", values_to = "imputed") %>%
  left_join(as_tibble(measured.norm, rownames = "uniprotid") %>%
              pivot_longer(-uniprotid, names_to = "sample", values_to = "measured"),
            by = c("uniprotid", "sample")) %>%
  mutate(was_imputed = is.na(measured)) %>%
  left_join(sample_sheet %>% select(sample = sample_id, construct, stimulation, replicate),
            by = "sample")
save_table(abundance_long, "abundance_long")

# SummarizedExperiment holding both assays
se <- SummarizedExperiment(
  assays = list(measured = measured.norm, imputed = imputed),
  colData = column_to_rownames(as.data.frame(sample_sheet), "sample_id")[colnames(imputed), ],
  rowData = prot.annotation[match(rownames(imputed), prot.annotation$uniprotid), ])
saveRDS(se, file.path(params_local$out_dir,
                      paste0(Sys.Date(), "_TurboID_IRE1.summarizedExperiment.RDS")))
message("Final object: ", nrow(se), " proteins x ", ncol(se), " samples, assays: ",
        paste(assayNames(se), collapse = ", "))
# Everything above is computed from `proteinGroups.txt`. The imputed assay is the sensitive screen that recovers the baits and the on/off bait- specific proteins and defines the candidate universe; the primary evidence is measured quantitation and reproducible on/off detection, with imputation-only significance retained as a supporting sensitivity tier. No protein was removed on biological grounds, and measured and imputed values are kept distinguishable throughout.

# Session info ---------------------------------------------------------------

# chunk: session-info
.si <- sessionInfo()
print(.si)
writeLines(capture.output(print(.si)),
           file.path(params_local$out_dir, "sessionInfo.txt"))

