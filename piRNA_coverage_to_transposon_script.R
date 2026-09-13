#!/usr/bin/env Rscript

# Tirant piRNA coverage across Drosophila strains
#
# Calculates strand-specific small-RNA coverage across a TE consensus,
# normalizes each library, averages coverage in 100-bp windows,
# and compares the profiles between strains.
#
# Run after setting:
#
# export TE_LIBRARY_PATH="/path/to/TE_base.fasta"
# export DATA_DIR="/path/to/alignment/files"
#
# DATA_DIR is optional and defaults to the current directory.
#
# scanBam() expects BAM-format input. If the files below are plain SAM files,
# convert them to BAM first.


library(Biostrings)
library(Rsamtools)
library(ggplot2)


# ------------------------------------------------------------------
# Input paths
# ------------------------------------------------------------------

TE_LIBRARY_PATH <- Sys.getenv("TE_LIBRARY_PATH")

if (TE_LIBRARY_PATH == "") {
  stop(
    "TE_LIBRARY_PATH is not set. Set it before running the script."
  )
}

DATA_DIR <- Sys.getenv(
  "DATA_DIR",
  unset = "."
)


# ------------------------------------------------------------------
# Coverage function
# ------------------------------------------------------------------

# Calculate sense and antisense coverage across one TE.
#
# bam_path is a data frame produced by bam_to_df().
# coeff is the normalization coefficient for the library.

plotcoverage_dots <- function(
    TE_name,
    TE_library_path,
    bam_path,
    plot_name,
    coeff
) {

  # load TE consensus library
  TE_library <- readDNAStringSet(
    TE_library_path
  )


  # TE length
  TE_length <- TE_library[
    names(TE_library) == TE_name
  ] |>
    width()


  TE_coverage <- data.frame(
    position = 1:TE_length,
    sense = 0,
    antisense = 0
  )


  # keep reads mapping to this TE
  intervals <- bam_path[
    bam_path$rname == TE_name,
    ,
    drop = FALSE
  ]


  # read abundance is encoded in the read name after "="
  intervals$count <- unlist(
    lapply(
      strsplit(intervals$qname, "="),
      function(x) x[2]
    )
  )

  intervals$count <- unlist(
    lapply(
      strsplit(intervals$count, ":"),
      function(x) x[1]
    )
  )


  # add each read to all positions covered by the alignment
  for (i in 1:nrow(intervals)) {

    read_positions <- intervals$pos[i]:(
      intervals$pos[i] +
        intervals$qwidth[i] -
        1
    )


    if (intervals$strand[i] == "+") {

      TE_coverage$sense[read_positions] <-
        TE_coverage$sense[read_positions] +
        as.numeric(intervals$count[i])

    } else {

      TE_coverage$antisense[read_positions] <-
        TE_coverage$antisense[read_positions] +
        as.numeric(intervals$count[i])

    }


    if ((i %% 1000) == 0) {

      print(
        paste(
          100 * i / nrow(intervals),
          "percent done"
        )
      )

    }
  }


  # long format
  TE_coverage_melted <- data.frame(

    position = rep(
      TE_coverage$position,
      2
    ),

    reads = c(
      TE_coverage$sense,
      TE_coverage$antisense
    ),

    label = c(
      rep(
        "sense",
        nrow(TE_coverage)
      ),
      rep(
        "antisense",
        nrow(TE_coverage)
      )
    )
  )


  # normalize coverage
  TE_coverage_melted$reads <-
    TE_coverage_melted$reads /
    coeff


  # antisense coverage is plotted below zero
  TE_coverage_melted$reads[
    TE_coverage_melted$label == "antisense"
  ] <-
    -TE_coverage_melted$reads[
      TE_coverage_melted$label == "antisense"
    ]


  # original 100-bp windowing scheme
  TE_coverage_melted$window <-
    round(
      TE_coverage_melted$position /
        100
    ) +
    1


  df_sense <- data.frame(
    window = 1:round(
      nrow(TE_coverage_melted) /
        200 +
        1
    ),
    strand = "sense",
    value = NA
  )


  df_antisense <- data.frame(
    window = 1:round(
      nrow(TE_coverage_melted) /
        200 +
        1
    ),
    strand = "antisense",
    value = NA
  )


  # mean coverage per window
  for (i in 1:nrow(df_sense)) {

    df_sense$value[i] <- mean(
      TE_coverage_melted[
        TE_coverage_melted$window == i &
          TE_coverage_melted$label == "sense",
      ]$reads
    )
  }


  for (i in 1:nrow(df_antisense)) {

    df_antisense$value[i] <- mean(
      TE_coverage_melted[
        TE_coverage_melted$window == i &
          TE_coverage_melted$label == "antisense",
      ]$reads
    )
  }


  rbind(
    df_sense,
    df_antisense
  )
}


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


# convert BAM output to a regular data frame
bam_to_df <- function(path) {

  bam <- scanBam(
    path
  )

  bam_field <- names(
    bam[[1]]
  )

  list_bam <- lapply(
    bam_field,
    function(y) {
      .unlist(
        lapply(
          bam,
          "[[",
          y
        )
      )
    }
  )

  bam_df <- do.call(
    "DataFrame",
    list_bam
  )

  names(
    bam_df
  ) <- bam_field

  data.frame(
    bam_df
  )
}


# ------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------

attp2_bam <- bam_to_df(
  file.path(
    DATA_DIR,
    "attp2_annotated.fa.gz.size23.fasta.sam.out"
  )
)

attp40_bam <- bam_to_df(
  file.path(
    DATA_DIR,
    "attp40_annotated.fa.gz.size23.fasta.sam.out"
  )
)

chr2_bam <- bam_to_df(
  file.path(
    DATA_DIR,
    "chr2balIfCuO_annotated.fa.gz.size23.fasta.sam.out"
  )
)

chr3_bam <- bam_to_df(
  file.path(
    DATA_DIR,
    "chr3balLyTm3Sb_annotated.fa.gz.size23.fasta.sam.out"
  )
)

ZH11_bam <- bam_to_df(
  file.path(
    DATA_DIR,
    "ZH11_annotated.fa.gz.size23.fasta.sam.out"
  )
)

iso1_bam <- bam_to_df(
  file.path(
    DATA_DIR,
    "iso1_annotated.fa.gz.size23.fasta.sam.out"
  )
)

DB_bam <- bam_to_df(
  file.path(
    DATA_DIR,
    "DBifcuosbtm3ser_annotated.fa.gz.size23.fasta.sam.out"
  )
)


# TE to analyze
TE_name <- "Tirant"


# ------------------------------------------------------------------
# Calculate coverage for each strain
# ------------------------------------------------------------------

# coefficients normalize the libraries to 1M miRNA reads

attp2 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  attp2_bam,
  "attp2 strain",
  0.534589
)


attp40 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  attp40_bam,
  "attp40 strain",
  0.62268
)


chr2 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  chr2_bam,
  "Chr2 balancer strain",
  1.08494
)


chr3 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  chr3_bam,
  "Chr3 balancer strain",
  0.55205
)


ZH <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  ZH11_bam,
  "ZH11 integrase strain",
  0.041453
)


iso1 <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  iso1_bam,
  "Iso1 strain",
  0.782243
)


DB <- plotcoverage_dots(
  TE_name,
  TE_LIBRARY_PATH,
  DB_bam,
  "Double balancer strain",
  0.660626
)


# add strain names
attp2$library <- "attp2"
attp40$library <- "attp40"
ZH$library <- "ZH11 integrase"
chr2$library <- "balancer chr2"
chr3$library <- "balancer chr3"
DB$library <- "double balancer"
iso1$library <- "Iso1"


# combine all libraries
summary <- rbind(
  attp2,
  attp40,
  ZH,
  chr2,
  chr3,
  DB,
  iso1
)


# ------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------

# sense is positive, antisense is negative
# color = library
# line type = strand

ggplot(
  summary,
  aes(
    x = window,
    y = value,
    linetype = strand,
    color = library
  )
) +
  geom_line(
    linewidth = 1.3
  ) +
  theme_bw() +
  xlab(
    "Genomic coordinate (100-bp window number)"
  ) +
  ylab(
    "Scaled coverage"
  ) +
  labs(
    color = "Library:",
    linetype = "Strand:"
  ) +
  theme(
    text = element_text(
      size = 12
    )
  ) +
  ggtitle(
    paste(
      "TE:",
      TE_name
    )
  )


# ZAM can be analyzed in the same way by changing:
# TE_name <- "ZAM"
