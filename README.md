# turboid-interactome

**A reproducible R pipeline for TurboID (and BioID) proximity-labeling proteomics.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.3390%2Fcells13090747-blue)](https://doi.org/10.3390/cells13090747)

This repository takes a raw MaxQuant `proteinGroups.txt` from A to Z: filtering, log2
transformation, per-sample normalization, a tiered missing-value and imputation
strategy, quality control, factorial `~ construct * stimulation` limma modeling,
an evidence-tiered candidate classification, and GO / pathway enrichment with
publication figures.

![Pipeline overview](docs/pipeline_overview.jpg)

> **Evidence framing.** Measured quantitation and reproducible on/off detection are the
> two **primary** evidence streams. The imputed factorial model is a sensitive screen
> that defines the candidate set and a sensitivity analysis; imputation-only
> significance is the weakest tier and never the sole basis for a conclusion.

It is demonstrated on the human IRE1α / IRE1β TurboID interactome in the mast
cell line HMC-1.2. The published study analyzed these data in Perseus; this
repository is an independent, fully scripted reimplementation in R that does not
use Perseus or any pre-processed intermediate. The only two inputs are the raw
MaxQuant file and a sample sheet.

## Biology in brief

IRE1 is an endoplasmic reticulum sensor of the unfolded protein response. Humans
have two isoforms, IRE1α (`ERN1`) and IRE1β (`ERN2`). Each was fused to TurboID,
a biotin ligase that tags proteins within a few nanometres, so streptavidin
capture followed by mass spectrometry reports each isoform's molecular
neighbourhood. A cytosolic V5-TurboID construct is the background control.

> ⚠️ **TurboID reports proximity, not physical binding:** a labeled protein is near
> the bait. Any direct-interaction claim needs an orthogonal experiment.

## Experimental design

Three constructs × two stimulation conditions × three biological replicates = **18
mass-spectrometry samples**.

| Construct | Gene | Role |
|-----------|------|------|
| IRE1α-TurboID | `ERN1` | bait |
| IRE1β-TurboID | `ERN2` | bait / isoform reference |
| V5-TurboID | — | cytosolic background control |

Stimulation is DMSO (vehicle) or tunicamycin (TM, ER stress). The stimulation
effect is small in these data, so DMSO and TM are analyzed both merged and as a
factorial with stimulation retained.

## Repository layout

```
turboid-interactome/
├── turboid_interactome.R      # the entire pipeline, one script, top to bottom
├── install_packages.R         # CRAN + Bioconductor dependency installer
├── metadata/
│   └── sample_sheet.csv       # the 18 LFQ columns and their biological meaning
├── data_raw/
│   └── proteinGroups.txt      # raw MaxQuant output (download, see below)
├── docs/
│   └── pipeline_overview.jpg  # pipeline schematic
└── results/                   # created on every run
    ├── tables/                # every table as CSV
    ├── figures/               # every figure as PDF + 300-dpi PNG + TIFF
    └── sessionInfo.txt
```

## Inputs

Everything is computed from two files:

- **`data_raw/proteinGroups.txt`** — raw MaxQuant output. This is **not** included you should added it by yourself in the repository.  and place the file at `data_raw/proteinGroups.txt`.
- **`metadata/sample_sheet.csv`** — the explicit, human-readable bridge between the
  MaxQuant `LFQ intensity ...` column names and their biological meaning. It is
  tracked in the repository and has one row per sample:

| column | meaning |
|--------|---------|
| `sample_id` | readable sample name, e.g. `IRE1a_DMSO_1` |
| `mq_lfq_column` | exact `LFQ intensity ...` header in `proteinGroups.txt` |
| `construct` | `IRE1a`, `IRE1b`, or `con` |
| `stimulation` | `DMSO` or `TM` |
| `replicate` | biological replicate number |

The script validates that every `mq_lfq_column` exists in `proteinGroups.txt` and
stops with an informative error if not, then derives the grouping and the control
flag internally.

## How to run

```r
# 1. install dependencies (CRAN + Bioconductor)
source("install_packages.R")

# 2. run the whole pipeline
#    from a shell:
Rscript turboid_interactome.R
#    or interactively in RStudio: open turboid_interactome.R and source it
```

The script creates `results/tables/` and `results/figures/` on every run and
writes all outputs there, plus a `SummarizedExperiment` `.RDS` and a
`sessionInfo.txt` for the record.

## How the pipeline works

1. **Filtering** — contaminants, reverse hits, only-by-site identifications, and
   protein groups with < 2 unique peptides are removed (staged and audited).
2. **Normalization** — zeros → NA, log2 transform, per-sample median centering
   (QC figures show before/after).
3. **Missing values** — missingness is shown to be left-censored, so a
   downshifted-normal imputation is applied per sample. Imputed and measured
   matrices are kept separate for the whole analysis.
4. **Differential analysis** — `limma` with moderated statistics (`trend`,
   `robust`): factorial `~ construct * stimulation` models on both the measured
   and the imputed matrix, plus a merged (DMSO + TM) reproduction of the
   manuscript's approach and a one-way ANOVA diagnostic.
5. **Candidate classification** — evidence-tiered (see below).
6. **Enrichment** — KEGG/Reactome via `g:Profiler`, and GO via *direct*
   annotations from `org.Hs.eg.db` (DAVID-style, avoiding abstract parent terms).
7. **Outputs** — all tables and all figures (vector PDF + 300-dpi PNG/TIFF).

## How to read the evidence tiers

Measured quantitation and reproducible on/off detection are the two primary
evidence streams; the imputed factorial model is a sensitive screen that defines
the candidate set and a sensitivity analysis. Each candidate is placed in one
tier:

| Tier | Meaning |
|------|---------|
| `high_confidence` | Significant in the **measured-data** limma analysis with ≥ 2 observed values in each group of the contrast — no call rests on a sparse comparison. |
| `qualitative_on_off` | Detected 3/3 replicates in a bait (in at least one stimulation) and absent or rare in the matched control (≤ 1/3). A measured test cannot score this because the control side is essentially empty. Assigned *before* the high-missingness downgrade so bait-specific on/off proteins are preserved. |
| `imputation_supported` | Significant only after imputation. Supporting evidence, explicitly labelled — never the sole basis for a conclusion. |
| `exploratory` | High missingness (> 50%) with neither measured nor on/off support. |
| `flagged_background` | Endogenously biotinylated carboxylases (e.g. `PC`, `ACACA`); always flagged, never promoted to candidate. |

Protein groups with more than one accession are flagged as **ambiguous** rather
than presented as a single gene product (marked `*` in the heatmap).

## Outputs

`results/tables/` holds every table as CSV (filter audit, per-protein missingness
and detection, all pairwise contrasts, the candidate confidence table, the master
significant-protein table with GO roles, and the enrichment results).

`results/figures/` holds every figure as a vector PDF plus a 300-dpi PNG and a
300-dpi LZW-compressed TIFF: per-sample distributions, normalization boxplots,
missingness policy, PCA, volcanoes, the candidate heatmap, MA plots, replicate
correlation and clustering, CV, and the pathway and GO figures.

## Adapting to another TurboID / BioID experiment

The engine is general, but several parameters are specific to this experiment. To
run it on your own bait and control design:

1. Edit the items marked **`EDIT FOR YOUR DESIGN`** in `params_local` (the bait
   gene symbols and the endogenous-biotin background list).
2. Provide a `sample_sheet.csv` with your own `construct` levels.
3. Adjust the contrast definitions and the on/off rule (currently written for
   three constructs named `IRE1a` / `IRE1b` / `con` with `n = 3` replicates) in
   the differential-analysis and candidate sections.

The downshifted-normal imputation, the valid-value rule, and the FDR and
fold-change thresholds are all in `params_local` at the top of the script.

## Reproducibility

Package versions are recorded two ways: `results/sessionInfo.txt` is written on
every run, and you can pin an exact environment with
[`renv`](https://rstudio.github.io/renv/) by running `renv::init()` once in the
project and committing the resulting `renv.lock`. R itself is not installed as
part of this repository; install the packages with `install_packages.R`.

## Data availability

Mass spectrometry proteomics data

## Citation

If you use this pipeline, please cite the associated article (see also
`CITATION.cff`):

> Ahmed N, Preisinger C, Wilhelm T, Huber M. TurboID-Based IRE1 Interactome
> Reveals Participants of the Endoplasmic Reticulum-Associated Protein Degradation
> Machinery in the Human Mast Cell Leukemia Cell Line HMC-1.2. *Cells* 2024;
> 13(9):747. https://doi.org/10.3390/cells13090747

## License

[MIT](LICENSE).
