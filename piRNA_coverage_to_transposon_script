###############################################################################
# Tirant piRNA coverage analysis across Drosophila strains
#
# Purpose
# -------
# This script calculates strand-specific small-RNA/piRNA coverage across a
# transposable-element consensus sequence (here: Tirant), normalizes coverage
# using sample-specific coefficients, summarizes coverage in genomic windows,
# and compares the resulting profiles between strains.
#
# The original analysis logic and sample-specific normalization coefficients
# are retained. Comments have been added to make each processing step explicit.
#
# Path privacy
# ------------
# The original script contained a local absolute path to the TE FASTA library.
# That path has been replaced with the TE_LIBRARY_PATH environment variable so
# the script can be shared publicly without exposing a local filesystem path.
#
# Before running, set for example:
#   export TE_LIBRARY_PATH="/path/to/TE_base.fasta"
#   export DATA_DIR="/path/to/alignment/files"   # optional; default is '.'
#
# NOTE: Rsamtools::scanBam() expects BAM-format input. The original input
# filenames ending in '.sam.out' are preserved below. If these files are
# plain-text SAM rather than BAM-format files with unusual names, convert them
# to BAM before running this script.
###############################################################################

library(Biostrings)
library(Rsamtools)
library(ggplot2)

###############################################################################
# User-configurable paths
###############################################################################

# Path to the FASTA file containing TE consensus/reference sequences.
# It is supplied through an environment variable to avoid hard-coding private
# paths in a script that may be uploaded to GitHub.
TE_LIBRARY_PATH <- Sys.getenv("TE_LIBRARY_PATH")

if (TE_LIBRARY_PATH == "") {
  stop(
    paste0(
      "TE_LIBRARY_PATH is not set. Set it before running this script, e.g. ",
      "export TE_LIBRARY_PATH='/path/to/TE_base.fasta'"
    )
  )
}

# Directory containing the BAM/alignment files. If DATA_DIR is not set, the
# current working directory is used.
DATA_DIR <- Sys.getenv("DATA_DIR", unset = ".")

###############################################################################
# Function: plotcoverage_dots
#
# Calculates strand-specific coverage across one TE reference sequence.
#
# Arguments
# ---------
# TE_name:
#   Name of the TE in the FASTA library (for example, "Tirant").
#
# TE_library_path:
#   Path to the FASTA file containing TE reference sequences.
#
# bam_path:
#   Despite the historical argument name, this function expects a DATA FRAME
#   produced by bam_to_df(), not a filesystem path.
#
# plot_name:
#   Plot label retained from the original function interface. It is currently
#   not used inside the function.
#
# coeff:
#   Sample-specific normalization coefficient. In the original analysis these
#   coefficients were used to scale coverage to 1 million miRNA reads.
#
# Returns
# -------
# A data frame with one row per coverage window and strand:
#   window | strand | value
###############################################################################

plotcoverage_dots <- function(TE_name, TE_library_path, bam_path, plot_name, coeff) {

  # Read the TE reference/consensus FASTA library.
  TE_library <- readDNAStringSet(TE_library_path)

  # Obtain the length of the requested TE sequence.
  # The function assumes that TE_name uniquely identifies a sequence in FASTA.
  TE_length <- TE_library[TE_library@ranges@NAMES == TE_name, ]@ranges@width

  # Initialize per-nucleotide coverage vectors for sense and antisense reads.
  TE_coverage <- data.frame(
    position = 1:TE_length,
    sense = 0,
    antisense = 0
  )

  # Keep only reads aligned to the TE being analyzed.
  intervals <- bam_path
  intervals <- intervals[intervals$rname == TE_name, ]

  # Recover the read count encoded in the query/read name.
  # The original naming convention is expected to contain '=' followed by the
  # count and then ':', e.g. something conceptually like "read=12:...".
  intervals$count <- unlist(
    lapply(strsplit(intervals$qname, "="), function(x) x[2])
  )
  intervals$count <- unlist(
    lapply(strsplit(intervals$count, ":"), function(x) x[1])
  )

  # Add each alignment's weighted read count to every TE position covered by
  # that alignment, keeping sense and antisense coverage separate.
  for (i in 1:nrow(intervals)) {

    # Explicitly define the aligned interval on the TE.
    # Parentheses around the end coordinate are important because ':' has high
    # operator precedence in R. This expresses the intended full read span.
    read_positions <- intervals$pos[i]:(
      intervals$pos[i] + intervals$qwidth[i] - 1
    )

    if (intervals$strand[i] == "+") {
      TE_coverage$sense[read_positions] <-
        TE_coverage$sense[read_positions] +
        as.numeric(as.character(intervals$count[i]))
    } else {
      TE_coverage$antisense[read_positions] <-
        TE_coverage$antisense[read_positions] +
        as.numeric(as.character(intervals$count[i]))
    }

    # Print progress for large alignment tables.
    if ((i %% 1000) == 0) {
      print(paste(100 * i / nrow(intervals), "percent done"))
    }
  }

  # Convert the two coverage columns into a long-format data frame that is
  # convenient for strand-specific plotting and summarization.
  TE_coverage_melted <- data.frame(
    position = rep(TE_coverage$position, 2),
    reads = c(TE_coverage$sense, TE_coverage$antisense),
    label = c(
      rep("sense", nrow(TE_coverage)),
      rep("antisense", nrow(TE_coverage))
    )
  )

  # Normalize coverage with the sample-specific coefficient.
  TE_coverage_melted$reads <- TE_coverage_melted$reads / coeff

  # Plot antisense signal below zero to visually separate the two strands.
  TE_coverage_melted$reads[
    TE_coverage_melted$label == "antisense"
  ] <- -1 * TE_coverage_melted$reads[
    TE_coverage_melted$label == "antisense"
  ]

  # Assign positions to the windowing scheme used in the original analysis.
  # NOTE: round(position / 100) + 1 is preserved from the original script.
  # If strict non-overlapping 100-bp bins are desired, this can instead be
  # changed to ceiling(position / 100) together with the window-count logic.
  TE_coverage_melted$window <- round(TE_coverage_melted$position / 100) + 1

  # Prepare one output table for each strand.
  df_sense <- data.frame(
    window = 1:round(nrow(TE_coverage_melted) / 200 + 1),
    strand = "sense",
    value = NA
  )
  df_antisense <- data.frame(
    window = 1:round(nrow(TE_coverage_melted) / 200 + 1),
    strand = "antisense",
    value = NA
  )

  # Calculate the mean normalized sense coverage within each window.
  for (i in 1:nrow(df_sense)) {
    df_sense$value[i] <- mean(
      TE_coverage_melted[
        as.numeric(as.character(TE_coverage_melted$window)) == i &
          TE_coverage_melted$label == "sense",
      ]$reads
    )
  }

  # Calculate the mean normalized antisense coverage within each window.
  for (i in 1:nrow(df_antisense)) {
    df_antisense$value[i] <- mean(
      TE_coverage_melted[
        as.numeric(as.character(TE_coverage_melted$window)) == i &
          TE_coverage_melted$label == "antisense",
      ]$reads
    )
  }

  # Combine the strand-specific summaries into one output data frame.
  df <- rbind(df_sense, df_antisense)
  return(df)
}

###############################################################################
# Helper: .unlist
#
# Combines fields returned by scanBam() while preserving factor levels.
# do.call(c, ...) can otherwise coerce factors to their underlying integers.
###############################################################################

.unlist <- function(x) {
  x1 <- x[[1L]]

  if (is.factor(x1)) {
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

###############################################################################
# Function: bam_to_df
#
# Reads a BAM-format alignment file with Rsamtools::scanBam() and converts the
# resulting list structure into a standard R data frame.
###############################################################################

bam_to_df <- function(path) {
  bam <- scanBam(path)

  # scanBam() returns a list; collect each BAM field across all chunks/records.
  bam_field <- names(bam[[1]])
  list_bam <- lapply(
    bam_field,
    function(y) .unlist(lapply(bam, "[[", y))
  )

  bam_df <- do.call("DataFrame", list_bam)
  names(bam_df) <- bam_field
  bam_df <- data.frame(bam_df)

  return(bam_df)
}

###############################################################################
# Load alignment data for each strain/library
#
# file.path() keeps the directory separate from the filenames, so DATA_DIR can
# point to a private location without embedding that location in this script.
###############################################################################

attp2_bam <- bam_to_df(
  file.path(DATA_DIR, "attp2_annotated.fa.gz.size23.fasta.sam.out")
)
attp40_bam <- bam_to_df(
  file.path(DATA_DIR, "attp40_annotated.fa.gz.size23.fasta.sam.out")
)
chr2_bam <- bam_to_df(
  file.path(DATA_DIR, "chr2balIfCuO_annotated.fa.gz.size23.fasta.sam.out")
)
chr3_bam <- bam_to_df(
  file.path(DATA_DIR, "chr3balLyTm3Sb_annotated.fa.gz.size23.fasta.sam.out")
)
ZH11_bam <- bam_to_df(
  file.path(DATA_DIR, "ZH11_annotated.fa.gz.size23.fasta.sam.out")
)
iso1_bam <- bam_to_df(
  file.path(DATA_DIR, "iso1_annotated.fa.gz.size23.fasta.sam.out")
)
DB_bam <- bam_to_df(
  file.path(DATA_DIR, "DBifcuosbtm3ser_annotated.fa.gz.size23.fasta.sam.out")
)

###############################################################################
# TE to analyze
###############################################################################

TE_name <- "Tirant"

###############################################################################
# Calculate normalized coverage profiles for all libraries
#
# The numeric values are the sample-specific normalization coefficients from
# the original analysis. The plot_name strings are retained for documentation,
# although plot_name is not currently used by plotcoverage_dots().
###############################################################################

attp2 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  attp2_bam,
  "attp2 strain\nnormalized to 1M miRNA reads",
  0.534589
)

attp40 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  attp40_bam,
  "attp40 strain\nnormalized to 1M miRNA reads",
  0.62268
)

chr2 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  chr2_bam,
  "Chr2 balancer strain\nnormalized to 1M miRNA reads",
  1.08494
)

chr3 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  chr3_bam,
  "Chr3 balancer strain\nnormalized to 1M miRNA reads",
  0.55205
)

ZH <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  ZH11_bam,
  "ZH11 (integrase) strain\nnormalized to 1M miRNA reads",
  0.041453
)

iso1 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  iso1_bam,
  "Iso1 strain\nnormalized to 1M miRNA reads",
  0.782243
)

DB <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  DB_bam,
  "Double balancer strain\nnormalized to 1M miRNA reads",
  0.660626
)

###############################################################################
# Add human-readable library labels before combining the profiles.
###############################################################################

attp2$library <- "attp2"
attp40$library <- "attp40"
ZH$library <- "ZH11 integrase"
chr2$library <- "balancer chr2"
chr3$library <- "balancer chr3"
DB$library <- "double balancer"
iso1$library <- "Iso1"

# Combine all strains into one table for plotting.
summary <- rbind(
  attp2,
  attp40,
  ZH,
  chr2,
  chr3,
  DB,
  iso1
)

###############################################################################
# Plot strand-specific TE coverage profiles.
#
# Sense coverage is positive and antisense coverage is negative, allowing both
# strands to be compared around the zero line. Color identifies the library;
# line type identifies the strand.
###############################################################################

ggplot(
  data = summary,
  aes(x = window, y = value, linetype = strand, col = library)
) +
  geom_line(size = 1.3) +
  theme_bw() +
  xlab("Genomic coordinate (100bp window number)") +
  ylab("Scaled coverage") +
  labs(col = "Library:", linetype = "Strand:") +
  theme(text = element_text(size = 12)) +
  ggtitle(paste("TE:", TE_name))

# ZAM could also be used as an additional example TE in this analysis.
