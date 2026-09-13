# Tirant 25-mer and genomic insertion analysis
#
# Main steps:
# - split the Tirant consensus into overlapping 25-mers
# - use alignment results to find Tirant regions that also match other TEs
# - map the remaining Tirant 25-mers to the OSC genome
# - identify larger Tirant-like regions
# - extract these regions from the genome
# - calculate coverage along Tirant from a BAM file
#
# Alignment itself is done separately. This script reads the resulting
# alignment files.


library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(IRanges)
library(ggplot2)


# ------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------

TE_LIBRARY_FILE <- "/path/to/TE_library.fasta"
OSC_GENOME_FILE <- "/path/to/OSC_genome.fasta"

TE_ALIGNMENT_FILE <- "/path/to/Mapped_Tirant_to_TEbase.out"
OSC_ALIGNMENT_FILE <- "/path/to/Mapped_Tirant_to_OSC_genome.out"

COVERAGE_BAM <- "/path/to/result.bam"

OUTPUT_DIR <- "/path/to/output"


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# analysis settings
TE_NAME <- "Tirant"
WINDOW_SIZE <- 25

# sequence added around clusters of Tirant 25-mer mappings
FLANK_SIZE <- 1000

# minimum size of the final candidate Tirant region
MIN_REGION_SIZE <- 3000


# chromosomes used in the OSC analysis
OSC_CHROMOSOMES <- c(
  "3L_RagTag",
  "3R_RagTag",
  "2L_RagTag",
  "2R_RagTag",
  "X_RagTag"
)


# ------------------------------------------------------------------
# Make overlapping Tirant 25-mers
# ------------------------------------------------------------------

TE_library <- readDNAStringSet(
  TE_LIBRARY_FILE
)


if (!TE_NAME %in% names(TE_library)) {

  stop(
    paste(
      TE_NAME,
      "was not found in the TE library."
    )
  )
}


Tirant <- TE_library[
  TE_NAME
]


Tirant_length <- width(
  Tirant
)


cat(
  "Tirant length:",
  Tirant_length,
  "bp\n"
)


if (Tirant_length < WINDOW_SIZE) {

  stop(
    "Tirant sequence is shorter than the requested window size."
  )
}


# overlapping windows:
# 1-25, 2-26, 3-27, ...
window_starts <- seq_len(
  Tirant_length - WINDOW_SIZE + 1
)


Tirant_25mers <- DNAStringSet(
  lapply(
    window_starts,
    function(i) {

      subseq(
        Tirant,
        start = i,
        width = WINDOW_SIZE
      )[[1]]
    }
  )
)


names(Tirant_25mers) <- paste0(
  "Tirant_seq_",
  window_starts
)


TIRANT_25MER_FASTA <- file.path(
  OUTPUT_DIR,
  "Tirant_25mers.fasta"
)


writeXStringSet(
  Tirant_25mers,
  TIRANT_25MER_FASTA
)


cat(
  "Number of Tirant 25-mers:",
  length(Tirant_25mers),
  "\n"
)


# ------------------------------------------------------------------
# Tirant 25-mers aligned to the TE library
# ------------------------------------------------------------------

# This file should contain the alignment of Tirant_25mers.fasta
# against the TE consensus library.

TE_hits <- read.table(
  TE_ALIGNMENT_FILE,
  fill = TRUE,
  stringsAsFactors = FALSE
)


# position of each 25-mer along Tirant
TE_hits$number <- vapply(

  strsplit(
    as.character(TE_hits$V1),
    "_"
  ),

  function(x) {

    as.numeric(
      tail(x, 1)
    )
  },

  numeric(1)
)


# Original analysis:
# keep non-perfect alignments to TEs other than Tirant.
#
# These positions are treated as potentially non-specific and removed
# from the genome analysis below.

TE_hits_other <- TE_hits[
  TE_hits$V13 != paste0("MD:Z:", WINDOW_SIZE) &
    TE_hits$V3 != TE_NAME,
  ,
  drop = FALSE
]


# Which parts of Tirant also align to other TEs?
ggplot(
  TE_hits_other,
  aes(
    x = number,
    fill = V3
  )
) +
  geom_bar(
    stat = "count"
  ) +
  theme_bw() +
  facet_wrap(
    ~ V3,
    scales = "free"
  ) +
  xlab(
    "Position along Tirant"
  ) +
  ylab(
    "Number of alignments"
  )


# Tirant positions that may give non-specific mappings
suspicious <- unique(
  TE_hits_other$number
)


cat(
  "Potentially non-specific Tirant 25-mers:",
  length(suspicious),
  "\n"
)


# ------------------------------------------------------------------
# Tirant 25-mers aligned to the OSC genome
# ------------------------------------------------------------------

OSC <- read.table(
  OSC_ALIGNMENT_FILE,
  fill = TRUE,
  stringsAsFactors = FALSE
)


# recover the original Tirant position from the 25-mer name
OSC$number <- vapply(

  strsplit(
    as.character(OSC$V1),
    "_"
  ),

  function(x) {

    as.numeric(
      tail(x, 1)
    )
  },

  numeric(1)
)


# mark perfect and non-perfect genome alignments
OSC$matches <- ifelse(
  OSC$V13 == paste0("MD:Z:", WINDOW_SIZE),
  "perfect",
  "mutations"
)


# remove Tirant positions that also map to other TE families
OSC <- OSC[
  !OSC$number %in% suspicious,
  ,
  drop = FALSE
]


# keep the main OSC chromosomes
OSC <- OSC[
  OSC$V3 %in% OSC_CHROMOSOMES,
  ,
  drop = FALSE
]


cat(
  "Genome mappings after filtering:",
  nrow(OSC),
  "\n"
)


# ------------------------------------------------------------------
# Find Tirant-like regions in the OSC genome
# ------------------------------------------------------------------

# each mapping corresponds to a 25-bp Tirant window
Tirant_windows <- GRanges(

  seqnames = as.character(
    OSC$V3
  ),

  ranges = IRanges(

    start = as.numeric(
      OSC$V4
    ),

    width = WINDOW_SIZE
  )
)


# merge overlapping 25-mer mappings
Tirant_windows <- reduce(
  Tirant_windows
)


cat(
  "Initial Tirant-like regions:",
  length(Tirant_windows),
  "\n"
)


# add 1 kb around each region
start(Tirant_windows) <- pmax(
  1,
  start(Tirant_windows) - FLANK_SIZE
)

end(Tirant_windows) <- (
  end(Tirant_windows) +
  FLANK_SIZE
)


# merge regions again after extending them
Tirant_windows <- reduce(
  Tirant_windows
)


# keep larger regions
Tirant_windows <- Tirant_windows[
  width(Tirant_windows) > MIN_REGION_SIZE
]


cat(
  "Candidate Tirant regions >",
  MIN_REGION_SIZE,
  "bp:",
  length(Tirant_windows),
  "\n"
)


# ------------------------------------------------------------------
# Extract candidate Tirant regions from the OSC genome
# ------------------------------------------------------------------

OSC_genome <- readDNAStringSet(
  OSC_GENOME_FILE
)


genome_lengths <- setNames(
  width(OSC_genome),
  names(OSC_genome)
)


# make sure chromosome names are present in the FASTA
missing_chr <- setdiff(
  unique(
    as.character(seqnames(Tirant_windows))
  ),
  names(OSC_genome)
)


if (length(missing_chr) > 0) {

  stop(
    paste(
      "Chromosomes missing from genome FASTA:",
      paste(
        missing_chr,
        collapse = ", "
      )
    )
  )
}


# prevent regions from extending beyond chromosome ends
end(Tirant_windows) <- pmin(
  end(Tirant_windows),
  genome_lengths[
    as.character(
      seqnames(Tirant_windows)
    )
  ]
)


Tirant_windows$sequence <- getSeq(
  OSC_genome,
  Tirant_windows
)


# save coordinates
Tirant_table <- as.data.frame(
  Tirant_windows
)


write.table(
  Tirant_table,
  file = file.path(
    OUTPUT_DIR,
    "candidate_Tirant_regions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# save candidate regions as FASTA
old_Tirant <- DNAStringSet(
  Tirant_windows$sequence
)


names(old_Tirant) <- paste0(
  seq_along(old_Tirant),
  ":",
  seqnames(Tirant_windows),
  ":",
  start(Tirant_windows),
  ":",
  width(Tirant_windows)
)


writeXStringSet(
  old_Tirant,
  file.path(
    OUTPUT_DIR,
    "candidate_Tirant_regions.fasta"
  )
)


cat(
  "Candidate Tirant sequences written:",
  length(old_Tirant),
  "\n"
)


# ------------------------------------------------------------------
# BAM helper functions
# ------------------------------------------------------------------

.unlist <- function(x) {

  x1 <- x[[1L]]

  if (is.factor(x1)) {

    structure(
      unlist(x),
      class = "factor",
      levels = levels(x1)
    )

  } else {

    do.call(
      c,
      x
    )
  }
}


bam_to_df <- function(path) {

  param <- ScanBamParam(
    what = c(
      "qname",
      "rname",
      "strand",
      "pos",
      "qwidth"
    )
  )


  bam <- scanBam(
    path,
    param = param
  )


  bam_fields <- names(
    bam[[1]]
  )


  bam_list <- lapply(

    bam_fields,

    function(field) {

      .unlist(
        lapply(
          bam,
          "[[",
          field
        )
      )
    }
  )


  bam_df <- do.call(
    "DataFrame",
    bam_list
  )


  names(
    bam_df
  ) <- bam_fields


  data.frame(
    bam_df
  )
}


# ------------------------------------------------------------------
# Read alignment used for Tirant coverage
# ------------------------------------------------------------------

tirant_bam <- bam_to_df(
  COVERAGE_BAM
)


# ------------------------------------------------------------------
# Get read abundance
# ------------------------------------------------------------------

# Collapsed small-RNA libraries often store read counts as:
#
# read_name=25
#
# or
#
# read_name=25:...
#
# If no count is present, use 1.

extract_read_count <- function(qname) {

  qname <- as.character(
    qname
  )


  count <- rep(
    1,
    length(qname)
  )


  has_count <- grepl(
    "=",
    qname,
    fixed = TRUE
  )


  if (any(has_count)) {

    encoded <- sub(
      ".*=",
      "",
      qname[has_count]
    )


    encoded <- sub(
      ":.*",
      "",
      encoded
    )


    encoded <- suppressWarnings(
      as.numeric(encoded)
    )


    encoded[
      is.na(encoded)
    ] <- 1


    count[
      has_count
    ] <- encoded
  }


  count
}


# ------------------------------------------------------------------
# Coverage along Tirant
# ------------------------------------------------------------------

plotcoverage_dots <- function(
    TE_name,
    TE_library_path,
    bam_df,
    coeff = 1,
    window_size = 100
) {

  TE_library <- readDNAStringSet(
    TE_library_path
  )


  if (!TE_name %in% names(TE_library)) {

    stop(
      paste(
        TE_name,
        "was not found in the TE library."
      )
    )
  }


  TE_length <- width(
    TE_library[
      TE_name
    ]
  )


  # reads mapping to Tirant
  intervals <- bam_df[
    as.character(bam_df$rname) == TE_name,
    ,
    drop = FALSE
  ]


  if (nrow(intervals) == 0) {

    warning(
      paste(
        "No reads mapped to",
        TE_name
      )
    )

    return(
      NULL
    )
  }


  intervals$count <- extract_read_count(
    intervals$qname
  )


  intervals$start <- as.integer(
    intervals$pos
  )


  intervals$end <- (
    intervals$start +
    as.integer(intervals$qwidth) -
    1
  )


  # keep coordinates inside the TE consensus
  intervals$start <- pmax(
    1,
    intervals$start
  )

  intervals$end <- pmin(
    TE_length,
    intervals$end
  )


  intervals <- intervals[
    intervals$start <= intervals$end,
    ,
    drop = FALSE
  ]


  # calculate coverage separately for both strands
  calculate_strand <- function(strand_value) {

    x <- intervals[
      as.character(intervals$strand) == strand_value,
      ,
      drop = FALSE
    ]


    if (nrow(x) == 0) {

      return(
        numeric(TE_length)
      )
    }


    ranges <- IRanges(
      start = x$start,
      end = x$end
    )


    as.numeric(
      coverage(
        ranges,
        weight = x$count,
        width = TE_length
      )
    )
  }


  sense <- calculate_strand(
    "+"
  )

  antisense <- calculate_strand(
    "-"
  )


  coverage_data <- data.frame(

    position = rep(
      seq_len(TE_length),
      2
    ),

    reads = c(
      sense,
      -antisense
    ),

    strand = rep(
      c(
        "sense",
        "antisense"
      ),
      each = TE_length
    )
  )


  # normalize coverage
  coverage_data$reads <- (
    coverage_data$reads /
    coeff
  )


  # average coverage in fixed-size windows
  coverage_data$window <- ceiling(
    coverage_data$position /
      window_size
  )


  coverage_summary <- aggregate(

    cbind(
      position,
      reads
    ) ~ window + strand,

    data = coverage_data,

    FUN = mean
  )


  return(
    coverage_summary
  )
}


# ------------------------------------------------------------------
# Example: Tirant coverage
# ------------------------------------------------------------------

Tirant_coverage <- plotcoverage_dots(

  TE_name = "Tirant",

  TE_library_path = TE_LIBRARY_FILE,

  bam_df = tirant_bam,

  coeff = 1,

  window_size = 100
)


if (!is.null(Tirant_coverage)) {

  ggplot(
    Tirant_coverage,
    aes(
      x = position,
      y = reads,
      linetype = strand
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3
    ) +
    geom_line(
      linewidth = 0.8
    ) +
    theme_bw() +
    xlab(
      "Position along Tirant"
    ) +
    ylab(
      "Coverage"
    ) +
    ggtitle(
      "Tirant coverage"
    )
}
