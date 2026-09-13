#!/usr/bin/env Rscript

# ==============================================================================
# piRNA-seq transposable-element analysis
# ==============================================================================
#
# Purpose:
#   Compare antisense piRNA abundance across transposable-element (TE) families
#   and Drosophila libraries/strains, with particular emphasis on Tirant and ZAM.
#
# Privacy / portability:
#   No local or institutional filesystem path is stored in this script.
#   The data directory must be supplied by the user either:
#
#     1) as the first command-line argument:
#          Rscript piRNAseq_commented.R /path/to/data
#
#        or
#
#     2) through an environment variable:
#          export PIRNASEQ_DATA_DIR="/path/to/data"
#          Rscript piRNAseq_commented.R
#
# Expected files inside DATA_DIR:
#
#   TE_RPKM_antisense-bowtie.txt
#   22libs/TE_RPKM_antisense-bowtie.txt
#
# Notes:
#   - The column selection and row range below are specific to the original
#     input tables and should be checked if the table format changes.
#   - RPKM values of zero cannot be displayed on a log10 axis and will be
#     omitted by ggplot2.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(ggplot2)
library(ggbeeswarm)


# ------------------------------------------------------------------------------
# 2. User-supplied data directory
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  DATA_DIR <- args[1]
} else {
  DATA_DIR <- Sys.getenv("PIRNASEQ_DATA_DIR")
}

if (!nzchar(DATA_DIR)) {
  stop(
    paste0(
      "No data directory was provided.\n",
      "Run the script as:\n",
      "  Rscript piRNAseq_commented.R /path/to/data\n",
      "or set the PIRNASEQ_DATA_DIR environment variable."
    )
  )
}

# Build input paths from the user-provided directory.
CURRENT_TE_FILE <- file.path(DATA_DIR, "TE_RPKM_antisense-bowtie.txt")
OLD_TE_FILE <- file.path(DATA_DIR, "22libs", "TE_RPKM_antisense-bowtie.txt")

if (!file.exists(CURRENT_TE_FILE)) {
  stop("Input file not found: TE_RPKM_antisense-bowtie.txt")
}

if (!file.exists(OLD_TE_FILE)) {
  stop("Historical input file not found: 22libs/TE_RPKM_antisense-bowtie.txt")
}


# ------------------------------------------------------------------------------
# 3. Read and prepare the current antisense piRNA table
# ------------------------------------------------------------------------------

# Read the table without converting character columns to factors.
te_data <- read.table(
  CURRENT_TE_FILE,
  header = FALSE,
  stringsAsFactors = FALSE
)

# Remove columns that are not used in this analysis.
# These column numbers correspond to the structure of the original input table.
te_data$V4 <- NULL
te_data$V5 <- NULL
te_data$V6 <- NULL
te_data$V9 <- NULL
te_data$V11 <- NULL
te_data$V12 <- NULL

# Rename the remaining columns to the corresponding libraries/strains.
colnames(te_data) <- c(
  "TE",
  "5906",
  "tjGal4",
  "23230",
  "DB6",
  "86Fb",
  "58A"
)

# Keep the TE rows used in the original analysis.
# This is dataset-specific and should be updated if the input table changes.
te_data <- te_data[2:123, , drop = FALSE]

# Remove suffixes following "_" from TE identifiers.
# Example: "Tirant_something" becomes "Tirant".
te_data$TE <- as.character(te_data$TE)
te_data$TE <- vapply(
  strsplit(te_data$TE, "_"),
  function(x) x[1],
  character(1)
)


# ------------------------------------------------------------------------------
# 4. Define TE clades/groups
# ------------------------------------------------------------------------------

gypsy_clade <- c(
  "gypsy5", "ZAM", "Tirant", "Idefix", "Quasimodo",
  "rover", "17.6", "297",
  "gypsy2", "gypsy3", "gypsy4", "gypsy6", "gypsy10",
  "gypsy", "gtwin", "springer", "HMS-Beagle2"
)

Tabor_group <- c(
  "412", "Tabor", "Stalker", "Stalker2", "mdg1", "blood"
)

Mdg3_group <- c(
  "mdg3", "micropia", "invader1", "invader2",
  "invader3", "invader4", "invader5", "invader6"
)

osvaldo_group <- c(
  "Circe", "gypsy8", "gypsy12", "Osvaldo"
)

bel_group <- c(
  "roo", "rooA", "Max-element", "diver", "diver2"
)

copia_group <- c(
  "copia", "Copia1", "1731"
)


# ------------------------------------------------------------------------------
# 5. Assign each TE to a clade
# ------------------------------------------------------------------------------

te_data$clade <- ifelse(
  te_data$TE %in% gypsy_clade,
  "gypsy_gypsy",
  ifelse(
    te_data$TE %in% Tabor_group,
    "Tabor_group",
    ifelse(
      te_data$TE %in% Mdg3_group,
      "Mdg3_group",
      ifelse(
        te_data$TE %in% osvaldo_group,
        "osvaldo_group",
        ifelse(
          te_data$TE %in% bel_group,
          "bel_group",
          ifelse(
            te_data$TE %in% copia_group,
            "copia_group",
            NA
          )
        )
      )
    )
  )
)

# Save TEs that were not assigned to one of the groups above.
# These are classified further at the end of the script.
te_rest <- te_data[is.na(te_data$clade), , drop = FALSE]

# Keep only the classified LTR-retrotransposon groups for the main analysis.
te_data <- te_data[!is.na(te_data$clade), , drop = FALSE]


# ------------------------------------------------------------------------------
# 6. Read the historical / reference library table
# ------------------------------------------------------------------------------

old_data <- read.table(
  OLD_TE_FILE,
  header = FALSE,
  stringsAsFactors = FALSE
)

# Keep only the TE name and libraries used in the original comparison.
old_data <- old_data[, c("V1", "V11", "V15", "V16")]

colnames(old_data) <- c(
  "TE",
  "Tj_old",
  "DB_old",
  "iso1"
)

# Standardize TE names in the same way as for the current table.
old_data$TE <- vapply(
  strsplit(as.character(old_data$TE), "_"),
  function(x) x[1],
  character(1)
)

# Keep only TEs belonging to the groups included in the main analysis.
old_data <- old_data[
  old_data$TE %in% c(
    gypsy_clade,
    Tabor_group,
    Mdg3_group,
    osvaldo_group,
    bel_group,
    copia_group
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# 7. Match old-library values to the current TE table
# ------------------------------------------------------------------------------

# Match by TE name rather than assuming that rows in the two input files are
# already in exactly the same order.
old_match <- match(te_data$TE, old_data$TE)

te_data$iso1 <- old_data$iso1[old_match]
te_data$Tj_old <- old_data$Tj_old[old_match]
te_data$DB_old <- old_data$DB_old[old_match]

# Convert values to numeric before plotting.
numeric_columns <- c(
  "5906", "tjGal4", "23230", "DB6", "86Fb", "58A",
  "iso1", "Tj_old", "DB_old"
)

te_data[numeric_columns] <- lapply(
  te_data[numeric_columns],
  function(x) as.numeric(as.character(x))
)


# ------------------------------------------------------------------------------
# 8. Optional quality-control comparison
# ------------------------------------------------------------------------------

# Compare the newer Tj-Gal4 measurements with the older Tj-Gal4 dataset.
# This is useful for visually checking consistency between datasets.
ggplot(
  data = te_data,
  aes(x = tjGal4, y = Tj_old)
) +
  geom_point() +
  theme_classic() +
  xlab("Tj-Gal4") +
  ylab("Tj-Gal4 old") +
  ggtitle("Comparison of current and historical Tj-Gal4 TE piRNA levels")


# ------------------------------------------------------------------------------
# 9. Convert the wide table into a long-format summary table
# ------------------------------------------------------------------------------

summary_data <- data.frame(
  TE = te_data$TE,
  data = te_data$iso1,
  clade = te_data$clade,
  lib = "iso1",
  stringsAsFactors = FALSE
)

append_library <- function(summary_df, values, library_name) {
  rbind(
    summary_df,
    data.frame(
      TE = te_data$TE,
      data = values,
      clade = te_data$clade,
      lib = library_name,
      stringsAsFactors = FALSE
    )
  )
}

summary_data <- append_library(summary_data, te_data$`86Fb`, "ZH-86Fb")
summary_data <- append_library(summary_data, te_data$`5906`, "BL-5906")
summary_data <- append_library(summary_data, te_data$`58A`, "ZH-58A")
summary_data <- append_library(summary_data, te_data$`23230`, "BL-23230")
summary_data <- append_library(summary_data, te_data$DB6, "DB")
summary_data <- append_library(summary_data, te_data$tjGal4, "Tj-Gal4")
summary_data <- append_library(summary_data, te_data$Tj_old, "Tj-Gal4_old")
summary_data <- append_library(summary_data, te_data$DB_old, "DB_old")

summary_data$data <- as.numeric(as.character(summary_data$data))


# ------------------------------------------------------------------------------
# 10. Highlight Tirant and ZAM
# ------------------------------------------------------------------------------

# Use category names rather than literal color names in the data.
# Colors are assigned later in scale_color_manual().
summary_data$highlight <- ifelse(
  summary_data$TE == "Tirant",
  "Tirant",
  ifelse(
    summary_data$TE == "ZAM",
    "ZAM",
    "Other"
  )
)

# Historical DB and Tj-Gal4 libraries are excluded from the final plot,
# reproducing the intent of the original analysis.
summary_data <- summary_data[
  !summary_data$lib %in% c("DB_old", "Tj-Gal4_old"),
  ,
  drop = FALSE
]

# Desired ordering of libraries on the x-axis.
level_order <- c(
  "iso1",
  "ZH-86Fb",
  "BL-5906",
  "ZH-58A",
  "BL-23230",
  "DB",
  "Tj-Gal4"
)

summary_data$lib <- factor(
  summary_data$lib,
  levels = level_order
)

# The final figure focuses on the gypsy and Tabor groups.
summary_plot <- summary_data[
  summary_data$clade %in% c("gypsy_gypsy", "Tabor_group"),
  ,
  drop = FALSE
]

# Mark non-highlighted TEs so connecting lines are drawn only for background TEs.
summary_plot$grey_dots <- !summary_plot$TE %in% c("Tirant", "ZAM")


# ------------------------------------------------------------------------------
# 11. Plot antisense piRNA abundance across libraries
# ------------------------------------------------------------------------------

# Background TEs are shown in grey and connected across libraries.
# Tirant and ZAM are highlighted as individual points.
#
# NOTE:
#   scale_y_log10() cannot represent values <= 0. ggplot2 will omit those values.
piRNA_plot <- ggplot(
  data = summary_plot,
  aes(x = lib, y = data)
) +
  geom_line(
    data = subset(summary_plot, grey_dots),
    aes(group = TE),
    linewidth = 0.5,
    color = "grey88",
    linetype = "dashed",
    alpha = 0.9
  ) +
  geom_beeswarm(
    aes(color = highlight),
    cex = 1.5,
    size = 3,
    alpha = 0.99
  ) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 100, 1000, 10000)
  ) +
  scale_color_manual(
    values = c(
      "Tirant" = "red",
      "ZAM" = "black",
      "Other" = "grey88"
    )
  ) +
  theme_classic() +
  xlab("") +
  ylab("Antisense piRNA levels, log10(RPKM)") +
  labs(color = NULL) +
  theme(
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )

print(piRNA_plot)


# ------------------------------------------------------------------------------
# 12. Additional TE classification for unassigned elements
# ------------------------------------------------------------------------------

# Expanded gypsy list retained from the original script for future analyses.
# It is not used to regenerate the main plot above.
gypsy_clade_extended <- c(
  "gypsy5", "ZAM", "Tirant", "accord", "accord2", "Idefix", "Quasimodo",
  "McClintock", "rover", "17.6", "297", "Transpac",
  "gypsy2", "gypsy3", "gypsy4", "gypsy6", "gypsy7", "gypsy9", "gypsy10",
  "Burdock", "gypsy", "gtwin", "springer", "HMS-Beagle",
  "HMS-Beagle2", "opus"
)

LINE_elements <- c(
  "TART-A", "Doc", "Doc3-element", "Doc4-element", "BS", "BS3", "BS4",
  "jockey2", "jockey", "Helena", "Ivk", "G4", "baggins", "G2", "R1-2",
  "R1A1-element", "R2-element", "G3", "Fw2", "Fw3", "G5", "G5A", "G6",
  "G7", "HeT-A", "TAHRE", "I-element"
)

DNA_elements <- c(
  "Transib5", "transib4", "transib3", "transib2", "mariner2", "hobo",
  "pogo", "P-element", "Helitron", "Bari1", "Bari2", "S2", "Tc1-2",
  "S-element", "Tc1", "looper1"
)

te_rest$clade <- ifelse(
  te_rest$TE %in% LINE_elements,
  "LINE",
  ifelse(
    te_rest$TE %in% DNA_elements,
    "DNA",
    NA
  )
)

# Objects available at the end of the script:
#   te_data      - classified TE table used for the main analysis
#   summary_data - long-format table across libraries
#   summary_plot - subset used in the final figure
#   piRNA_plot   - ggplot object for the final figure
#   te_rest      - remaining TEs with LINE/DNA annotations where applicable
