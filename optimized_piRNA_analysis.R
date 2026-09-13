# ==============================================================================
# Transposon small-RNA targeting visualization
# ==============================================================================
#
# Purpose:
#   Visualize small-RNA coverage along transposable-element (TE) consensus
#   sequences.
#
#   Sense reads are plotted above zero.
#   Antisense reads are plotted below zero.
#
#   Coverage is normalized using a user-provided normalization coefficient
#   (for example, normalization to 1 million miRNA reads).
#
# Workflow:
#   1. Read small-RNA alignments.
#   2. Select reads mapping to a TE.
#   3. Calculate sense and antisense coverage along the TE.
#   4. Normalize coverage.
#   5. Average coverage in fixed-size windows.
#   6. Compare targeting profiles between libraries.
#
# IMPORTANT:
#   No private/local paths are included in this script.
#   Set your own paths in the USER SETTINGS section.
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

library(Biostrings)
library(Rsamtools)
library(IRanges)
library(ggplot2)


# ==============================================================================
# 2. USER SETTINGS
# ==============================================================================

# Directory containing mapped small-RNA alignment files.
ALIGNMENT_DIR <- "/path/to/alignment/files"

# FASTA file containing TE consensus sequences.
TE_LIBRARY_FILE <- "/path/to/TE_library.fasta"

# Directory where plots should be saved.
OUTPUT_DIR <- "/path/to/output"


# Size of the windows used to smooth coverage.
WINDOW_SIZE <- 100


# Create output directory if necessary.
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. SAMPLE INFORMATION
# ==============================================================================

# Define alignment files, display names and normalization coefficients.
#
# The normalization coefficient should correspond to the normalization strategy
# used for the experiment (for example, reads per 1 million miRNA reads).

samples <- data.frame(

  sample = c(
    "GFP",
    "Panx",
    "Piwi",
    "Nxf2_1",
    "Nxf2_2"
  ),

  label = c(
    "Control KD",
    "Panx KD",
    "Piwi KD",
    "Nxf2 KD (1)",
    "Nxf2 KD (2)"
  ),

  file = c(
    "GFP3_mapped.out",
    "Panx23_mapped.out",
    "Piwi23_mapped.out",
    "Nxf2123_mapped.out",
    "Nxf2223_mapped.out"
  ),

  normalization = c(
    1.06818,
    0.827488,
    2.4986,
    0.98335,
    1.09356
  ),

  stringsAsFactors = FALSE
)


# Construct full file paths.
samples$path <- file.path(
  ALIGNMENT_DIR,
  samples$file
)


# ==============================================================================
# 4. LOAD TE LIBRARY
# ==============================================================================

TE_library <- readDNAStringSet(
  TE_LIBRARY_FILE
)


cat(
  "Loaded",
  length(TE_library),
  "TE consensus sequences\n"
)


# ==============================================================================
# 5. HELPER FUNCTION: CONVERT BAM TO DATA FRAME
# ==============================================================================

# scanBam() returns a nested list.
# This helper preserves factor levels while converting the output into
# a standard R data frame.

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

  if (!file.exists(path)) {
    stop(
      paste(
        "Alignment file not found:",
        basename(path)
      )
    )
  }


  # Only fields required for the coverage analysis are loaded.
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


  bam_df <- data.frame(
    bam_df
  )


  return(
    bam_df
  )
}


# ==============================================================================
# 6. LOAD ALIGNMENT FILES
# ==============================================================================

# Load every small-RNA library once.
#
# The resulting list can then be reused for any number of TEs.

bam_data <- setNames(

  lapply(
    samples$path,
    bam_to_df
  ),

  samples$sample
)


# ==============================================================================
# 7. EXTRACT READ COPY NUMBER
# ==============================================================================

# In this dataset, collapsed reads contain their abundance in the read name.
#
# Example:
#
#   read_name=25:...
#
# means that the sequence represents 25 reads.
#
# If no count is encoded in the read name, the read is treated as one read.

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


  encoded_count <- sub(
    ".*=",
    "",
    qname[has_count]
  )


  encoded_count <- sub(
    ":.*",
    "",
    encoded_count
  )


  encoded_count <- suppressWarnings(
    as.numeric(encoded_count)
  )


  # If a count cannot be parsed, use one read.
  encoded_count[
    is.na(encoded_count)
  ] <- 1


  count[
    has_count
  ] <- encoded_count


  return(
    count
  )
}


# ==============================================================================
# 8. CALCULATE SMALL-RNA COVERAGE ALONG ONE TE
# ==============================================================================

calculate_te_coverage <- function(
    TE_name,
    bam_df,
    normalization = 1,
    window_size = 100
) {

  # ---------------------------------------------------------------------------
  # Check input
  # ---------------------------------------------------------------------------

  if (!TE_name %in% names(TE_library)) {

    stop(
      paste(
        "TE not found in TE library:",
        TE_name
      )
    )
  }


  if (normalization <= 0) {

    stop(
      "Normalization coefficient must be greater than zero."
    )
  }


  # ---------------------------------------------------------------------------
  # Determine TE length
  # ---------------------------------------------------------------------------

  TE_index <- match(
    TE_name,
    names(TE_library)
  )


  TE_length <- width(
    TE_library
  )[TE_index]


  # ---------------------------------------------------------------------------
  # Select reads mapping to the TE
  # ---------------------------------------------------------------------------

  intervals <- bam_df[
    as.character(bam_df$rname) == TE_name,
    ,
    drop = FALSE
  ]


  # If no reads map to this TE, return zero coverage.
  if (nrow(intervals) == 0) {

    positions <- seq_len(
      TE_length
    )


    empty <- data.frame(

      position = rep(
        positions,
        2
      ),

      value = 0,

      strand = rep(
        c(
          "sense",
          "antisense"
        ),
        each = TE_length
      )
    )


    empty$window <- ceiling(
      empty$position /
        window_size
    )


    result <- aggregate(
      cbind(
        position,
        value
      ) ~ window + strand,
      data = empty,
      FUN = mean
    )


    return(
      result
    )
  }


  # ---------------------------------------------------------------------------
  # Extract abundance of collapsed reads
  # ---------------------------------------------------------------------------

  intervals$count <- extract_read_count(
    intervals$qname
  )


  # ---------------------------------------------------------------------------
  # Calculate alignment coordinates
  # ---------------------------------------------------------------------------

  intervals$start <- as.integer(
    intervals$pos
  )


  intervals$end <- (
    intervals$start +
    as.integer(intervals$qwidth) -
    1
  )


  # Ensure coordinates do not extend outside the TE consensus sequence.
  intervals$start <- pmax(
    1,
    intervals$start
  )


  intervals$end <- pmin(
    TE_length,
    intervals$end
  )


  # Remove invalid intervals.
  intervals <- intervals[
    !is.na(intervals$start) &
      !is.na(intervals$end) &
      intervals$start <= intervals$end,
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------------------------
  # Calculate strand-specific coverage
  # ---------------------------------------------------------------------------

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


    coverage_vector <- coverage(

      ranges,

      weight = x$count,

      width = TE_length
    )


    return(
      as.numeric(
        coverage_vector
      )
    )
  }


  sense_coverage <- calculate_strand(
    "+"
  )


  antisense_coverage <- calculate_strand(
    "-"
  )


  # ---------------------------------------------------------------------------
  # Build long-format coverage table
  # ---------------------------------------------------------------------------

  coverage_data <- data.frame(

    position = rep(
      seq_len(TE_length),
      2
    ),

    value = c(
      sense_coverage,
      -antisense_coverage
    ),

    strand = rep(
      c(
        "sense",
        "antisense"
      ),
      each = TE_length
    )
  )


  # Normalize coverage.
  coverage_data$value <- (
    coverage_data$value /
    normalization
  )


  # ---------------------------------------------------------------------------
  # Average coverage in fixed-size windows
  # ---------------------------------------------------------------------------

  # Window 1 = positions 1-100
  # Window 2 = positions 101-200
  # etc.

  coverage_data$window <- ceiling(
    coverage_data$position /
      window_size
  )


  coverage_summary <- aggregate(

    cbind(
      position,
      value
    ) ~ window + strand,

    data = coverage_data,

    FUN = mean
  )


  return(
    coverage_summary
  )
}


# ==============================================================================
# 9. PLOT ONE TE ACROSS MULTIPLE LIBRARIES
# ==============================================================================

plot_te_targeting <- function(
    TE_name,
    sample_ids = samples$sample,
    window_size = WINDOW_SIZE,
    show_legend = TRUE
) {

  # Check that requested sample names exist.
  if (!all(sample_ids %in% samples$sample)) {

    stop(
      "One or more requested samples are not present in the sample table."
    )
  }


  coverage_list <- lapply(

    sample_ids,

    function(sample_id) {

      sample_info <- samples[
        samples$sample == sample_id,
        ,
        drop = FALSE
      ]


      df <- calculate_te_coverage(

        TE_name = TE_name,

        bam_df = bam_data[
          [sample_id]
        ],

        normalization = sample_info$normalization,

        window_size = window_size
      )


      df$library <- sample_info$label


      return(
        df
      )
    }
  )


  summary_data <- do.call(
    rbind,
    coverage_list
  )


  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------

  p <- ggplot(

    summary_data,

    aes(
      x = position,
      y = value,
      color = library,
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

    labs(

      title = paste(
        "TE:",
        TE_name
      ),

      x = paste0(
        "Position along TE (",
        window_size,
        "-bp bins)"
      ),

      y = "Scaled small-RNA coverage",

      color = "Library:",

      linetype = "Strand:"
    ) +

    theme(

      text = element_text(
        size = 12
      ),

      plot.title = element_text(
        face = "bold"
      )
    )


  if (!show_legend) {

    p <- p +

      theme(
        legend.position = "none"
      )
  }


  return(
    p
  )
}


# ==============================================================================
# 10. EXAMPLE: CONTROL VERSUS PANX KD
# ==============================================================================

plot_te_targeting(

  TE_name = "Stalker",

  sample_ids = c(
    "GFP",
    "Panx"
  )
)


# ==============================================================================
# 11. EXAMPLE: COMPARE ALL LIBRARIES FOR ONE TE
# ==============================================================================

Stalker_plot <- plot_te_targeting(

  TE_name = "Stalker",

  sample_ids = c(
    "GFP",
    "Panx",
    "Piwi",
    "Nxf2_1",
    "Nxf2_2"
  )
)


print(
  Stalker_plot
)


# Save single-TE plot.

ggsave(

  filename = file.path(
    OUTPUT_DIR,
    "Stalker_smallRNA_targeting.pdf"
  ),

  plot = Stalker_plot,

  width = 9,

  height = 6,

  units = "in"
)


# ==============================================================================
# 12. LIST OF TRANSPOSONS TO VISUALIZE
# ==============================================================================

TE_LIST <- c(

  "412",
  "gtwin",
  "mdg1",
  "Idefix",
  "rover",
  "gypsy5",
  "gypsy",
  "17.6",
  "Stalker",
  "blood",
  "Tabor",
  "gypsy10_new",
  "springer",
  "Stalker2",
  "gypsy2_new",
  "Quasimodo",
  "gypsy4",
  "HMS-Beagle2",
  "ZAM",
  "gypsy6_new",
  "copia",
  "F-element",
  "roo",
  "297",
  "Doc",
  "HMS-Beagle",
  "baggins",
  "TART-A",
  "flea",
  "gypsy3_new",
  "Helena",
  "Transpac",
  "G6",
  "jockey2",
  "Circe",
  "Copia1",
  "Repbase_DNAREP1",
  "I-element_new",
  "Cr1a",
  "mdg3",
  "rooA",
  "1360",
  "gypsy12",
  "opus",
  "Juan"
)


# ==============================================================================
# 13. CREATE MULTI-PAGE PDF FOR MANY TEs
# ==============================================================================

# Instead of manually creating obj1, obj2, obj3, ...,
# loop over the TE list.
#
# Each TE is written to a separate page of the PDF.
#
# This is much easier to inspect than one extremely tall figure.

pdf(

  file = file.path(
    OUTPUT_DIR,
    "Transposon_smallRNA_targeting_profiles.pdf"
  ),

  width = 8,

  height = 5
)


for (TE_name in TE_LIST) {


  # Skip TE names that are absent from the FASTA library.
  if (!TE_name %in% names(TE_library)) {

    warning(
      paste(
        "Skipping TE not found in library:",
        TE_name
      )
    )

    next
  }


  p <- plot_te_targeting(

    TE_name = TE_name,

    sample_ids = c(
      "GFP",
      "Panx"
    )
  )


  print(
    p
  )
}


dev.off()


# ==============================================================================
# 14. OPTIONAL: ALL LIBRARIES FOR A SELECTED SET OF TEs
# ==============================================================================

# For example, compare Control, Panx, Piwi and Nxf2 depletion
# for selected elements.

SELECTED_TES <- c(
  "Stalker",
  "ZAM",
  "gypsy",
  "mdg1",
  "roo"
)


pdf(

  file = file.path(
    OUTPUT_DIR,
    "Selected_TEs_all_libraries.pdf"
  ),

  width = 9,

  height = 6
)


for (TE_name in SELECTED_TES) {


  if (!TE_name %in% names(TE_library)) {

    warning(
      paste(
        "Skipping TE not found in library:",
        TE_name
      )
    )

    next
  }


  p <- plot_te_targeting(

    TE_name = TE_name,

    sample_ids = c(
      "GFP",
      "Panx",
      "Piwi",
      "Nxf2_1",
      "Nxf2_2"
    )
  )


  print(
    p
  )
}


dev.off()


# ==============================================================================
# End of script
# ==============================================================================
