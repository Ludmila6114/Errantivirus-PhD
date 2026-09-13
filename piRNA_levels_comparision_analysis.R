#!/usr/bin/env Rscript

# piRNA-seq TE analysis
#
# Compare antisense piRNA levels across TE families and Drosophila libraries.
# Tirant and ZAM are highlighted in the final plot.
#
# Run with:
#   Rscript piRNAseq.R /path/to/data
#
# or:
#   export PIRNASEQ_DATA_DIR="/path/to/data"
#   Rscript piRNAseq.R
#
# Expected files:
#   TE_RPKM_antisense-bowtie.txt
#   22libs/TE_RPKM_antisense-bowtie.txt
#
# Column selection and row numbers below are specific to these input tables.
# Values <= 0 cannot be shown on the log10 scale.


library(ggplot2)
library(ggbeeswarm)


# ------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------

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
      "Run as:\n",
      "  Rscript piRNAseq.R /path/to/data\n",
      "or set PIRNASEQ_DATA_DIR."
    )
  )
}


CURRENT_TE_FILE <- file.path(
  DATA_DIR,
  "TE_RPKM_antisense-bowtie.txt"
)

OLD_TE_FILE <- file.path(
  DATA_DIR,
  "22libs",
  "TE_RPKM_antisense-bowtie.txt"
)


if (!file.exists(CURRENT_TE_FILE)) {
  stop("Input file not found: TE_RPKM_antisense-bowtie.txt")
}

if (!file.exists(OLD_TE_FILE)) {
  stop("Input file not found: 22libs/TE_RPKM_antisense-bowtie.txt")
}


# ------------------------------------------------------------------
# Current antisense piRNA dataset
# ------------------------------------------------------------------

te_data <- read.table(
  CURRENT_TE_FILE,
  header = FALSE,
  stringsAsFactors = FALSE
)


# remove columns that are not used
te_data$V4 <- NULL
te_data$V5 <- NULL
te_data$V6 <- NULL
te_data$V9 <- NULL
te_data$V11 <- NULL
te_data$V12 <- NULL


# library names
colnames(te_data) <- c(
  "TE",
  "5906",
  "tjGal4",
  "23230",
  "DB6",
  "86Fb",
  "58A"
)


# rows containing the TE table
te_data <- te_data[
  2:123,
  ,
  drop = FALSE
]


# simplify TE names
# e.g. Tirant_something -> Tirant
te_data$TE <- as.character(
  te_data$TE
)

te_data$TE <- vapply(
  strsplit(te_data$TE, "_"),
  function(x) x[1],
  character(1)
)


# ------------------------------------------------------------------
# TE groups
# ------------------------------------------------------------------

gypsy_clade <- c(
  "gypsy5",
  "ZAM",
  "Tirant",
  "Idefix",
  "Quasimodo",
  "rover",
  "17.6",
  "297",
  "gypsy2",
  "gypsy3",
  "gypsy4",
  "gypsy6",
  "gypsy10",
  "gypsy",
  "gtwin",
  "springer",
  "HMS-Beagle2"
)


Tabor_group <- c(
  "412",
  "Tabor",
  "Stalker",
  "Stalker2",
  "mdg1",
  "blood"
)


Mdg3_group <- c(
  "mdg3",
  "micropia",
  "invader1",
  "invader2",
  "invader3",
  "invader4",
  "invader5",
  "invader6"
)


osvaldo_group <- c(
  "Circe",
  "gypsy8",
  "gypsy12",
  "Osvaldo"
)


bel_group <- c(
  "roo",
  "rooA",
  "Max-element",
  "diver",
  "diver2"
)


copia_group <- c(
  "copia",
  "Copia1",
  "1731"
)


# assign TE groups
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


# keep unclassified TEs for later
te_rest <- te_data[
  is.na(te_data$clade),
  ,
  drop = FALSE
]


# use the selected LTR groups for the main analysis
te_data <- te_data[
  !is.na(te_data$clade),
  ,
  drop = FALSE
]


# ------------------------------------------------------------------
# Older/reference dataset
# ------------------------------------------------------------------

old_data <- read.table(
  OLD_TE_FILE,
  header = FALSE,
  stringsAsFactors = FALSE
)


# keep only the libraries needed here
old_data <- old_data[
  ,
  c(
    "V1",
    "V11",
    "V15",
    "V16"
  )
]


colnames(old_data) <- c(
  "TE",
  "Tj_old",
  "DB_old",
  "iso1"
)


# simplify TE names in the same way
old_data$TE <- vapply(
  strsplit(
    as.character(old_data$TE),
    "_"
  ),
  function(x) x[1],
  character(1)
)


# keep the same TE groups
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


# ------------------------------------------------------------------
# Match old and current datasets
# ------------------------------------------------------------------

# match by TE name instead of assuming the row order is identical
old_match <- match(
  te_data$TE,
  old_data$TE
)


te_data$iso1 <- old_data$iso1[
  old_match
]

te_data$Tj_old <- old_data$Tj_old[
  old_match
]

te_data$DB_old <- old_data$DB_old[
  old_match
]


# convert expression columns to numeric
numeric_columns <- c(
  "5906",
  "tjGal4",
  "23230",
  "DB6",
  "86Fb",
  "58A",
  "iso1",
  "Tj_old",
  "DB_old"
)


te_data[numeric_columns] <- lapply(
  te_data[numeric_columns],
  function(x) as.numeric(as.character(x))
)


# quick comparison of old and current Tj-Gal4 data
ggplot(
  te_data,
  aes(
    x = tjGal4,
    y = Tj_old
  )
) +
  geom_point() +
  theme_classic() +
  xlab("Tj-Gal4") +
  ylab("Tj-Gal4 old") +
  ggtitle(
    "Current vs old Tj-Gal4 TE piRNA levels"
  )


# ------------------------------------------------------------------
# Make a long table for plotting
# ------------------------------------------------------------------

summary_data <- data.frame(
  TE = te_data$TE,
  data = te_data$iso1,
  clade = te_data$clade,
  lib = "iso1",
  stringsAsFactors = FALSE
)


# helper for adding another library
append_library <- function(
    summary_df,
    values,
    library_name
) {

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


summary_data <- append_library(
  summary_data,
  te_data$`86Fb`,
  "ZH-86Fb"
)

summary_data <- append_library(
  summary_data,
  te_data$`5906`,
  "BL-5906"
)

summary_data <- append_library(
  summary_data,
  te_data$`58A`,
  "ZH-58A"
)

summary_data <- append_library(
  summary_data,
  te_data$`23230`,
  "BL-23230"
)

summary_data <- append_library(
  summary_data,
  te_data$DB6,
  "DB"
)

summary_data <- append_library(
  summary_data,
  te_data$tjGal4,
  "Tj-Gal4"
)

summary_data <- append_library(
  summary_data,
  te_data$Tj_old,
  "Tj-Gal4_old"
)

summary_data <- append_library(
  summary_data,
  te_data$DB_old,
  "DB_old"
)


summary_data$data <- as.numeric(
  as.character(summary_data$data)
)


# ------------------------------------------------------------------
# Highlight Tirant and ZAM
# ------------------------------------------------------------------

summary_data$highlight <- ifelse(
  summary_data$TE == "Tirant",
  "Tirant",
  ifelse(
    summary_data$TE == "ZAM",
    "ZAM",
    "Other"
  )
)


# old DB and Tj-Gal4 libraries are not used in the final figure
summary_data <- summary_data[
  !summary_data$lib %in% c(
    "DB_old",
    "Tj-Gal4_old"
  ),
  ,
  drop = FALSE
]


# order libraries on the x-axis
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


# plot gypsy and Tabor groups
summary_plot <- summary_data[
  summary_data$clade %in% c(
    "gypsy_gypsy",
    "Tabor_group"
  ),
  ,
  drop = FALSE
]


# background TEs will be connected with grey lines
summary_plot$grey_dots <- !summary_plot$TE %in% c(
  "Tirant",
  "ZAM"
)


# ------------------------------------------------------------------
# Plot antisense piRNA levels
# ------------------------------------------------------------------

# Tirant = red
# ZAM = black
# other TEs = grey
#
# Values <= 0 are not shown on the log10 scale.

piRNA_plot <- ggplot(
  summary_plot,
  aes(
    x = lib,
    y = data
  )
) +

  geom_line(
    data = subset(
      summary_plot,
      grey_dots
    ),
    aes(
      group = TE
    ),
    linewidth = 0.5,
    color = "grey88",
    linetype = "dashed",
    alpha = 0.9
  ) +

  geom_beeswarm(
    aes(
      color = highlight
    ),
    cex = 1.5,
    size = 3,
    alpha = 0.99
  ) +

  scale_y_log10(
    breaks = c(
      0.1,
      1,
      10,
      100,
      1000,
      10000
    )
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

  ylab(
    "Antisense piRNA levels, log10(RPKM)"
  ) +

  labs(
    color = NULL
  ) +

  theme(
    axis.text = element_text(
      size = 16
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    legend.position = "bottom"
  )


print(
  piRNA_plot
)


# ------------------------------------------------------------------
# Other TE classes
# ------------------------------------------------------------------

# longer gypsy list kept for other analyses
gypsy_clade_extended <- c(
  "gypsy5",
  "ZAM",
  "Tirant",
  "accord",
  "accord2",
  "Idefix",
  "Quasimodo",
  "McClintock",
  "rover",
  "17.6",
  "297",
  "Transpac",
  "gypsy2",
  "gypsy3",
  "gypsy4",
  "gypsy6",
  "gypsy7",
  "gypsy9",
  "gypsy10",
  "Burdock",
  "gypsy",
  "gtwin",
  "springer",
  "HMS-Beagle",
  "HMS-Beagle2",
  "opus"
)


LINE_elements <- c(
  "TART-A",
  "Doc",
  "Doc3-element",
  "Doc4-element",
  "BS",
  "BS3",
  "BS4",
  "jockey2",
  "jockey",
  "Helena",
  "Ivk",
  "G4",
  "baggins",
  "G2",
  "R1-2",
  "R1A1-element",
  "R2-element",
  "G3",
  "Fw2",
  "Fw3",
  "G5",
  "G5A",
  "G6",
  "G7",
  "HeT-A",
  "TAHRE",
  "I-element"
)


DNA_elements <- c(
  "Transib5",
  "transib4",
  "transib3",
  "transib2",
  "mariner2",
  "hobo",
  "pogo",
  "P-element",
  "Helitron",
  "Bari1",
  "Bari2",
  "S2",
  "Tc1-2",
  "S-element",
  "Tc1",
  "looper1"
)


# classify remaining TEs as LINE or DNA where possible
te_rest$clade <- ifelse(
  te_rest$TE %in% LINE_elements,
  "LINE",
  ifelse(
    te_rest$TE %in% DNA_elements,
    "DNA",
    NA
  )
)


# main objects:
# te_data      - TE table used in the analysis
# summary_data - long table with all libraries
# summary_plot - subset used for the final plot
# piRNA_plot   - final plot
# te_rest      - remaining LINE/DNA elements
