# ==============================================================================
# RepeatMasker and Roo transposable-element analysis
# ==============================================================================
#
# This script:
#
#   1. Reads RepeatMasker annotations from the OSC genome
#   2. Calculates the percentage of each TE consensus represented by each hit
#   3. Identifies Roo insertions
#   4. Tests whether genomic regions around Roo insertions overlap other TEs
#   5. Compares TE copy numbers between the OSC genome and dm6
#   6. Examines Roo nucleotide composition in 50-bp windows
#   7. Extracts Roo insertion coordinates and sequences
#   8. Generates BED files
#   9. Summarises sequencing reads mapping to Roo insertion loci
#
# IMPORTANT:
# No local filesystem paths are stored in this public version.
# Set the paths below before running the analysis.
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

library(Biostrings)
library(GenomicRanges)
library(Rsamtools)

library(ggplot2)
library(dplyr)
library(tidyr)


# ==============================================================================
# 2. USER SETTINGS
# ==============================================================================

# --------------------------------------------------------------------------
# Replace these placeholders with paths on your own computer/server.
#
# Do NOT commit your private paths to GitHub.
# --------------------------------------------------------------------------

OSC_REPEATMASKER_FILE <- "/path/to/OSC_repeatmasker.out"

DM6_REPEATMASKER_FILE <- "/path/to/dm6_repeatmasker.out"

TE_LIBRARY_FILE <- "/path/to/TE_library.fasta"

OSC_GENOME_FILE <- "/path/to/OSC_genome.fasta"

INSERTION_ANALYSIS_DIR <- "/path/to/insertion_analysis"

OUTPUT_DIR <- "/path/to/output"


# Create output directory if it does not already exist.

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. ANALYSIS PARAMETERS
# ==============================================================================

# TE of primary interest.

TARGET_TE <- "roo"


# Minimum percentage of TE consensus sequence represented by a RepeatMasker hit.
#
# The original OSC analysis used 40%.

OSC_MIN_PERCENT <- 40


# Threshold used to define approximately full-length TE copies.

FULL_LENGTH_PERCENT <- 90


# Number of bp surrounding each Roo insertion used for overlap analysis.

FLANK_SIZE <- 5000


# Window size for nucleotide-composition analysis.

SEQUENCE_WINDOW <- 50


# ==============================================================================
# 4. READ TE CONSENSUS LIBRARY
# ==============================================================================

TE_library <- readDNAStringSet(
  TE_LIBRARY_FILE
)


# Check that the target TE exists.

if (!TARGET_TE %in% names(TE_library)) {

  stop(
    paste(
      "TE",
      TARGET_TE,
      "was not found in the TE library."
    )
  )
}


# ==============================================================================
# 5. HELPER FUNCTION FOR REPEATMASKER OUTPUT
# ==============================================================================

read_repeatmasker <- function(path) {

  # RepeatMasker .out files normally contain header rows.
  # This reproduces the structure used in the original analysis.

  RM <- read.table(
    path,
    fill = TRUE,
    stringsAsFactors = FALSE
  )


  # Remove the first two non-data rows.

  if (nrow(RM) >= 3) {

    RM <- RM[
      3:nrow(RM),
      ,
      drop = FALSE
    ]

  }


  # Remove incomplete rows.

  RM <- na.omit(RM)


  return(RM)
}


# ==============================================================================
# 6. REPEATMASKER ANALYSIS — OSC GENOME
# ==============================================================================

RM_OSC <- read_repeatmasker(
  OSC_REPEATMASKER_FILE
)


# ------------------------------------------------------------------------------
# Determine TE consensus length
# ------------------------------------------------------------------------------

# RepeatMasker column V10 contains the TE family name.
#
# Match each RepeatMasker hit to the corresponding sequence in the
# TE consensus library.

RM_OSC$TE_length <- width(TE_library)[
  match(
    as.character(RM_OSC$V10),
    names(TE_library)
  )
]


# Remove RepeatMasker hits whose TE is absent from the TE library.

RM_OSC <- RM_OSC[
  !is.na(RM_OSC$TE_length),
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# Calculate percentage of TE consensus represented by each hit
# ------------------------------------------------------------------------------

RM_OSC$range <- (
  as.numeric(as.character(RM_OSC$V7)) -
  as.numeric(as.character(RM_OSC$V6))
)


RM_OSC$percent <- (
  RM_OSC$range /
  RM_OSC$TE_length *
  100
)


# Keep hits covering at least the requested fraction of the TE consensus.

RM_OSC <- RM_OSC[
  RM_OSC$percent >= OSC_MIN_PERCENT,
  ,
  drop = FALSE
]


# ==============================================================================
# 7. EXTRACT ROO REPEATMASKER HITS
# ==============================================================================

roo <- RM_OSC[
  RM_OSC$V10 == TARGET_TE,
  ,
  drop = FALSE
]


cat(
  "Number of Roo hits with >=",
  OSC_MIN_PERCENT,
  "% consensus coverage:",
  nrow(roo),
  "\n"
)


# Save the RepeatMasker records.

write.table(
  roo,
  file = file.path(
    OUTPUT_DIR,
    "Roo_insertions_repeatmasker.out"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ==============================================================================
# 8. TE COPY-NUMBER SUMMARY IN OSC GENOME
# ==============================================================================

TEs_OSC <- data.frame(
  table(RM_OSC$V10)
)


# Sort from most abundant to least abundant.

TEs_OSC <- TEs_OSC[
  order(
    TEs_OSC$Freq,
    decreasing = TRUE
  ),
]


# Plot TE abundance.

ggplot(
  data = TEs_OSC,
  aes(
    x = reorder(Var1, Freq),
    y = Freq
  )
) +

  geom_bar(
    stat = "identity",
    fill = "red"
  ) +

  coord_flip() +

  theme_bw() +

  xlab("TE name") +

  ylab(
    paste0(
      "Number of copies with at least ",
      OSC_MIN_PERCENT,
      "% length"
    )
  ) +

  theme(
    text = element_text(
      size = 12
    )
  ) +

  ggtitle(
    "OSC genome — TE statistics"
  )


# ==============================================================================
# 9. CREATE GENOMICRANGES FOR ROO
# ==============================================================================

Roo_GR <- GRanges(

  seqnames = as.character(
    roo$V5
  ),

  ranges = IRanges(

    start = as.numeric(
      as.character(roo$V6)
    ),

    end = as.numeric(
      as.character(roo$V7)
    )

  ),

  strand = ifelse(
    as.character(roo$V9) == "+",
    "+",
    "-"
  ),

  Element = as.character(
    roo$V10
  ),

  percent = roo$percent
)


# Keep approximately full-length Roo copies.

Roo_GR <- Roo_GR[
  Roo_GR$percent > FULL_LENGTH_PERCENT
]


cat(
  "Number of Roo copies >",
  FULL_LENGTH_PERCENT,
  "% length:",
  length(Roo_GR),
  "\n"
)


# ==============================================================================
# 10. CREATE GENOMICRANGES FOR ALL OTHER TEs
# ==============================================================================

All_TEs <- GRanges(

  seqnames = as.character(
    RM_OSC$V5
  ),

  ranges = IRanges(

    start = as.numeric(
      as.character(RM_OSC$V6)
    ),

    end = as.numeric(
      as.character(RM_OSC$V7)
    )

  ),

  strand = ifelse(
    as.character(RM_OSC$V9) == "+",
    "+",
    "-"
  ),

  Element = as.character(
    RM_OSC$V10
  ),

  percent = RM_OSC$percent
)


# Exclude Roo itself.

All_TEs <- All_TEs[
  All_TEs$Element != TARGET_TE
]


cat(
  "Number of non-Roo TE hits:",
  length(All_TEs),
  "\n"
)


# ==============================================================================
# 11. ANALYSE REGIONS SURROUNDING ROO INSERTIONS
# ==============================================================================

# Create +/- 5-kb regions around the Roo loci.
#
# pmax(1, ...) prevents genomic coordinates from becoming negative.

Roo_flanks <- GRanges(

  seqnames = as.character(
    roo$V5
  ),

  ranges = IRanges(

    start = pmax(
      1,
      as.numeric(as.character(roo$V6)) -
        FLANK_SIZE
    ),

    end = (
      as.numeric(as.character(roo$V7)) +
      FLANK_SIZE
    )

  ),

  strand = ifelse(
    as.character(roo$V9) == "+",
    "+",
    "-"
  ),

  Element = as.character(
    roo$V10
  ),

  percent = roo$percent
)


cat(
  "Number of Roo flanking regions:",
  length(Roo_flanks),
  "\n"
)


# Identify Roo regions whose surrounding sequence overlaps another TE.

Roo_overlap_regions <- subsetByOverlaps(

  Roo_flanks,
  All_TEs,

  ignore.strand = TRUE
)


cat(
  "Roo regions overlapping another TE:",
  length(Roo_overlap_regions),
  "\n"
)


# Interactive inspection in RStudio if desired:
#
# View(
#   as.data.frame(Roo_overlap_regions)
# )


# ==============================================================================
# 12. REPEATMASKER ANALYSIS — dm6
# ==============================================================================

RM_dm6 <- read_repeatmasker(
  DM6_REPEATMASKER_FILE
)


# Match RepeatMasker TE names to TE consensus sequences.

RM_dm6$TE_length <- width(TE_library)[
  match(
    as.character(RM_dm6$V10),
    names(TE_library)
  )
]


# Remove elements not represented in the TE library.

RM_dm6 <- RM_dm6[
  !is.na(RM_dm6$TE_length),
  ,
  drop = FALSE
]


# Calculate TE hit length.

RM_dm6$range <- (
  as.numeric(as.character(RM_dm6$V7)) -
  as.numeric(as.character(RM_dm6$V6))
)


# Percentage of TE consensus represented.

RM_dm6$percent <- (
  RM_dm6$range /
  RM_dm6$TE_length *
  100
)


# Keep approximately full-length copies.

RM_dm6_full <- RM_dm6[
  RM_dm6$percent >= FULL_LENGTH_PERCENT,
  ,
  drop = FALSE
]


TEs_dm6 <- data.frame(
  table(
    RM_dm6_full$V10
  )
)


# ==============================================================================
# 13. FULL-LENGTH TE COPY NUMBER IN OSCs
# ==============================================================================

# For the dm6-versus-OSC comparison, use the same 90% cutoff for both genomes.
#
# This avoids comparing >=90% dm6 copies against >=40% OSC copies.

RM_OSC_full <- RM_OSC[
  RM_OSC$percent >= FULL_LENGTH_PERCENT,
  ,
  drop = FALSE
]


TEs_OSC_full <- data.frame(
  table(
    RM_OSC_full$V10
  )
)


# ==============================================================================
# 14. COMPARE TE COPY NUMBER BETWEEN dm6 AND OSC GENOMES
# ==============================================================================

transposon_summary <- data.frame(

  transposon = names(
    TE_library
  ),

  dm6 = NA,

  OSCs = NA,

  stringsAsFactors = FALSE
)


# Add dm6 counts.

transposon_summary$dm6 <- TEs_dm6$Freq[
  match(
    transposon_summary$transposon,
    TEs_dm6$Var1
  )
]


# Add OSC counts.

transposon_summary$OSCs <- TEs_OSC_full$Freq[
  match(
    transposon_summary$transposon,
    TEs_OSC_full$Var1
  )
]


# TEs with no detected copies receive a count of zero.

transposon_summary$dm6[
  is.na(transposon_summary$dm6)
] <- 0


transposon_summary$OSCs[
  is.na(transposon_summary$OSCs)
] <- 0


# Remove TEs absent from both genomes.

transposon_summary <- transposon_summary[
  transposon_summary$dm6 +
    transposon_summary$OSCs != 0,
  ,
  drop = FALSE
]


# ==============================================================================
# 15. CONVERT COPY-NUMBER TABLE TO LONG FORMAT
# ==============================================================================

n_TE <- nrow(
  transposon_summary
)


df <- data.frame(

  TE_name = rep(
    transposon_summary$transposon,
    2
  ),

  genome = c(
    rep("dm6", n_TE),
    rep("OSCs", n_TE)
  ),

  value = c(
    transposon_summary$dm6,
    transposon_summary$OSCs
  )
)


# ==============================================================================
# 16. PLOT dm6 VERSUS OSC TE COPY NUMBER
# ==============================================================================

ggplot(
  data = df,
  aes(
    x = TE_name,
    y = value,
    fill = genome
  )
) +

  geom_bar(
    stat = "identity",
    position = "dodge"
  ) +

  coord_flip() +

  theme_bw() +

  ggtitle(
    paste0(
      "RepeatMasker\nTEs with at least ",
      FULL_LENGTH_PERCENT,
      "% length"
    )
  ) +

  labs(
    fill = "Genome:"
  ) +

  xlab("") +

  ylab(
    "Number of full-length copies"
  ) +

  theme(
    text = element_text(
      size = 12
    )
  )


# ==============================================================================
# 17. ROO SEQUENCE COMPOSITION
# ==============================================================================

# Extract Roo consensus sequence from the TE library.

roo_TE <- TE_library[
  names(TE_library) == TARGET_TE
]


if (length(roo_TE) != 1) {

  stop(
    "Expected exactly one Roo sequence in the TE library."
  )
}


roo_sequence <- as.character(
  roo_TE
)


# Convert sequence into individual nucleotides.

roo_bases <- unlist(
  strsplit(
    roo_sequence,
    ""
  )
)


sequence_df <- data.frame(

  base = roo_bases,

  position = seq_along(
    roo_bases
  )
)


# Assign each nucleotide to a 50-bp window.
#
# ceiling(position / 50) gives:
#
#   1–50    -> window 1
#   51–100  -> window 2
#   etc.

sequence_df$window <- ceiling(
  sequence_df$position /
    SEQUENCE_WINDOW
)


# ==============================================================================
# 18. CALCULATE A/T/G/C CONTENT PER WINDOW
# ==============================================================================

ATGC <- sequence_df %>%

  count(
    window,
    base
  ) %>%

  complete(

    window = seq_len(
      max(sequence_df$window)
    ),

    base = c(
      "A",
      "T",
      "G",
      "C"
    ),

    fill = list(
      n = 0
    )
  ) %>%

  group_by(
    window
  ) %>%

  mutate(

    percent = (
      n /
      sum(n) *
      100
    )

  ) %>%

  ungroup()


# ==============================================================================
# 19. PLOT ROO NUCLEOTIDE COMPOSITION
# ==============================================================================

ggplot(
  data = ATGC,
  aes(
    x = window,
    y = percent,
    color = base
  )
) +

  geom_line() +

  theme_bw() +

  theme(
    text = element_text(
      size = 14
    )
  ) +

  ggtitle(
    "Roo element sequence composition"
  ) +

  xlab(
    paste0(
      SEQUENCE_WINDOW,
      "-bp window number"
    )
  ) +

  ylab(
    "Percent"
  ) +

  labs(
    color = "Base:"
  )


# ==============================================================================
# 20. EXTRACT ROO INSERTION SEQUENCES FROM OSC GENOME
# ==============================================================================

# Load OSC genome.

genome_OSC <- readDNAStringSet(
  OSC_GENOME_FILE
)


# Extract genomic sequences corresponding to full-length Roo insertions.

Roo_GR$sequence <- getSeq(
  genome_OSC,
  Roo_GR
)


# Optional inspection:
#
# View(
#   as.data.frame(Roo_GR)
# )


# ==============================================================================
# 21. WRITE SIMPLE THREE-COLUMN ROO BED FILE
# ==============================================================================

bed_roo_simple <- as.data.frame(
  Roo_GR
)


bed_roo_simple <- bed_roo_simple[
  ,
  1:3
]


write.table(

  bed_roo_simple,

  file = file.path(
    OUTPUT_DIR,
    "Roo_insertions_simple.bed"
  ),

  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,

  sep = "\t"
)


# ==============================================================================
# 22. CREATE STANDARD BED6 FILE
# ==============================================================================

# BED6 format:
#
# chromosome
# start
# end
# name
# score
# strand
#
# BED coordinates use a 0-based start coordinate.
# RepeatMasker coordinates are 1-based, therefore start - 1 is used here.

roo_start <- as.numeric(
  as.character(
    roo$V6
  )
)


roo_end <- as.numeric(
  as.character(
    roo$V7
  )
)


my_bed <- data.frame(

  chr = as.character(
    roo$V5
  ),

  start = roo_start - 1,

  end = roo_end,

  name = paste(
    roo$V5,
    roo$V6,
    roo$V7,
    sep = ":"
  ),

  score = 0,

  strand = ifelse(
    roo$V9 == "+",
    "+",
    "-"
  )
)


write.table(

  my_bed,

  file = file.path(
    OUTPUT_DIR,
    "Roo_insertions.bed"
  ),

  row.names = FALSE,

  col.names = FALSE,

  quote = FALSE,

  sep = "\t"
)


# ==============================================================================
# 23. ROO INSERTION IDENTIFICATION FROM ALIGNMENT FILES
# ==============================================================================

# ------------------------------------------------------------------------------
# Helper function:
# safely combine scanBam output into a data frame
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# Convert BAM alignment into data frame
# ------------------------------------------------------------------------------

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


  bam_df <- data.frame(
    bam_df
  )


  return(
    bam_df
  )
}


# ==============================================================================
# 24. ALIGNMENT FILE LOCATIONS
# ==============================================================================

# Change filenames here if required.
#
# NOTE:
# scanBam() expects BAM-format input.
# If these files are plain-text SAM files, convert them to BAM first.

GFP4_FILE <- file.path(
  INSERTION_ANALYSIS_DIR,
  "GFP4_align_to_163.sam.MAPPED"
)


MAEL4_FILE <- file.path(
  INSERTION_ANALYSIS_DIR,
  "Mael4_align_to_163.sam.MAPPED"
)


# ==============================================================================
# 25. GFP4
# ==============================================================================

GFP4 <- bam_to_df(
  GFP4_FILE
)


# Extract read count encoded after "=" in the read name.

GFP4$count <- vapply(

  strsplit(
    as.character(GFP4$qname),
    "=",
    fixed = TRUE
  ),

  function(x) {

    if (length(x) >= 2) {

      x[2]

    } else {

      NA_character_

    }

  },

  character(1)
)


GFP4$count <- as.numeric(
  GFP4$count
)


# Sum read counts for each Roo insertion.

GFP4_summary <- GFP4 %>%

  group_by(
    rname
  ) %>%

  summarise(

    Freq = sum(
      count,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ==============================================================================
# 26. MAEL4
# ==============================================================================

Mael4 <- bam_to_df(
  MAEL4_FILE
)


Mael4$count <- vapply(

  strsplit(
    as.character(Mael4$qname),
    "=",
    fixed = TRUE
  ),

  function(x) {

    if (length(x) >= 2) {

      x[2]

    } else {

      NA_character_

    }

  },

  character(1)
)


Mael4$count <- as.numeric(
  Mael4$count
)


Mael4_summary <- Mael4 %>%

  group_by(
    rname
  ) %>%

  summarise(

    Freq = sum(
      count,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ==============================================================================
# 27. IDENTIFY ROO INSERTIONS PRESENT IN MAEL4 BUT NOT GFP4
# ==============================================================================

Mael4_specific_insertions <- unique(

  Mael4_summary[
    !Mael4_summary$rname %in%
      GFP4_summary$rname,
  ]$rname

)


print(
  Mael4_specific_insertions
)


# ==============================================================================
# 28. SAVE SUMMARY TABLES
# ==============================================================================

write.table(

  transposon_summary,

  file = file.path(
    OUTPUT_DIR,
    "TE_copy_number_dm6_vs_OSC.tsv"
  ),

  sep = "\t",

  row.names = FALSE,

  quote = FALSE
)


write.table(

  GFP4_summary,

  file = file.path(
    OUTPUT_DIR,
    "GFP4_Roo_insertions.tsv"
  ),

  sep = "\t",

  row.names = FALSE,

  quote = FALSE
)


write.table(

  Mael4_summary,

  file = file.path(
    OUTPUT_DIR,
    "Mael4_Roo_insertions.tsv"
  ),

  sep = "\t",

  row.names = FALSE,

  quote = FALSE
)


# ==============================================================================
# End of analysis
# ==============================================================================
