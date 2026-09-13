# small-RNA targeting along transposable-element consensus sequences
#
# sense reads are plotted above zero
# antisense reads are plotted below zero
#
# coverage is normalized with sample-specific coefficients and
# averaged in fixed-size windows


library(Biostrings)
library(Rsamtools)
library(IRanges)
library(ggplot2)


# ------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------

# folder with mapped small-RNA files
ALIGNMENT_DIR <- "/path/to/alignment/files"

# TE consensus FASTA
TE_LIBRARY_FILE <- "/path/to/TE_library.fasta"

# output folder
OUTPUT_DIR <- "/path/to/output"

# window size for coverage smoothing
WINDOW_SIZE <- 100


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------
# Samples
# ------------------------------------------------------------------

# normalization coefficients correspond to the normalization used
# for each small-RNA library

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


samples$path <- file.path(
  ALIGNMENT_DIR,
  samples$file
)


# ------------------------------------------------------------------
# TE library
# ------------------------------------------------------------------

TE_library <- readDNAStringSet(
  TE_LIBRARY_FILE
)


cat(
  "Loaded",
  length(TE_library),
  "TE consensus sequences\n"
)


# ------------------------------------------------------------------
# BAM helper functions
# ------------------------------------------------------------------

# scanBam() returns lists; this keeps factor levels intact
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


# convert BAM data to a regular data frame
bam_to_df <- function(path) {

  if (!file.exists(path)) {

    stop(
      paste(
        "Alignment file not found:",
        basename(path)
      )
    )
  }


  # only load fields needed for the coverage analysis
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


# load all libraries once
bam_data <- setNames(

  lapply(
    samples$path,
    bam_to_df
  ),

  samples$sample
)


# ------------------------------------------------------------------
# Read counts
# ------------------------------------------------------------------

# collapsed reads contain their abundance in the read name
#
# example:
# read_name=25:...
#
# if no count is present, the read is counted as 1

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


# ------------------------------------------------------------------
# Coverage along one TE
# ------------------------------------------------------------------

calculate_te_coverage <- function(
    TE_name,
    bam_df,
    normalization = 1,
    window_size = 100
) {

  # check TE name
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


  # TE length
  TE_index <- match(
    TE_name,
    names(TE_library)
  )


  TE_length <- width(
    TE_library
  )[TE_index]


  # reads mapping to this TE
  intervals <- bam_df[
    as.character(bam_df$rname) == TE_name,
    ,
    drop = FALSE
  ]


  # return zero coverage if there are no mapped reads
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


  # get abundance of collapsed reads
  intervals$count <- extract_read_count(
    intervals$qname
  )


  # alignment coordinates
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
    !is.na(intervals$start) &
      !is.na(intervals$end) &
      intervals$start <= intervals$end,
    ,
    drop = FALSE
  ]


  # calculate coverage separately for each strand
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


    as.numeric(
      coverage_vector
    )
  }


  sense_coverage <- calculate_strand(
    "+"
  )


  antisense_coverage <- calculate_strand(
    "-"
  )


  # combine sense and antisense coverage
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


  # normalize coverage
  coverage_data$value <- (
    coverage_data$value /
    normalization
  )


  # 100-bp windows:
  # 1-100, 101-200, 201-300, ...
  coverage_data$window <- ceiling(
    coverage_data$position /
      window_size
  )


  # average coverage within each window
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


# ------------------------------------------------------------------
# Plot one TE
# ------------------------------------------------------------------

plot_te_targeting <- function(
    TE_name,
    sample_ids = samples$sample,
    window_size = WINDOW_SIZE,
    show_legend = TRUE
) {

  # check sample names
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


# ------------------------------------------------------------------
# Example: control vs Panx KD
# ------------------------------------------------------------------

plot_te_targeting(

  TE_name = "Stalker",

  sample_ids = c(
    "GFP",
    "Panx"
  )
)


# ------------------------------------------------------------------
# Example: compare all libraries for Stalker
# ------------------------------------------------------------------

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


# ------------------------------------------------------------------
# TE list
# ------------------------------------------------------------------

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


# ------------------------------------------------------------------
# Plot many TEs
# ------------------------------------------------------------------

# one TE per page
pdf(

  file = file.path(
    OUTPUT_DIR,
    "Transposon_smallRNA_targeting_profiles.pdf"
  ),

  width = 8,
  height = 5
)


for (TE_name in TE_LIST) {

  # skip TEs that are not present in the FASTA
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


# ------------------------------------------------------------------
# Selected TEs with all libraries
# ------------------------------------------------------------------

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
