# RepeatMasker / Roo analysis
#
# Main steps:
# - find Roo and other TE copies in the OSC genome
# - compare TE copy number between OSCs and dm6
# - look at TEs around Roo insertions
# - check Roo sequence composition
# - extract Roo coordinates/sequences
# - compare reads mapping to Roo insertion sites


library(Biostrings)
library(GenomicRanges)
library(Rsamtools)

library(ggplot2)
library(dplyr)
library(tidyr)


# ------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------

OSC_REPEATMASKER_FILE <- "/path/to/OSC_repeatmasker.out"
DM6_REPEATMASKER_FILE <- "/path/to/dm6_repeatmasker.out"

TE_LIBRARY_FILE <- "/path/to/TE_library.fasta"
OSC_GENOME_FILE <- "/path/to/OSC_genome.fasta"

INSERTION_ANALYSIS_DIR <- "/path/to/insertion_analysis"
OUTPUT_DIR <- "/path/to/output"


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# analysis settings
TARGET_TE <- "roo"

# keep OSC RepeatMasker hits covering at least 40% of the TE consensus
OSC_MIN_PERCENT <- 40

# use 90% as the cutoff for approximately full-length elements
FULL_LENGTH_PERCENT <- 90

# region around Roo insertions used for overlap analysis
FLANK_SIZE <- 5000

# window size for Roo sequence composition
SEQUENCE_WINDOW <- 50


# ------------------------------------------------------------------
# TE library
# ------------------------------------------------------------------

TE_library <- readDNAStringSet(
  TE_LIBRARY_FILE
)


if (!TARGET_TE %in% names(TE_library)) {

  stop(
    paste(
      TARGET_TE,
      "was not found in the TE library."
    )
  )
}


# ------------------------------------------------------------------
# Read RepeatMasker output
# ------------------------------------------------------------------

read_repeatmasker <- function(path) {

  RM <- read.table(
    path,
    fill = TRUE,
    stringsAsFactors = FALSE
  )


  # first two rows are RepeatMasker headers
  if (nrow(RM) >= 3) {

    RM <- RM[
      3:nrow(RM),
      ,
      drop = FALSE
    ]
  }


  RM <- na.omit(
    RM
  )


  return(
    RM
  )
}


# ------------------------------------------------------------------
# OSC RepeatMasker data
# ------------------------------------------------------------------

RM_OSC <- read_repeatmasker(
  OSC_REPEATMASKER_FILE
)


# V10 contains the TE name
# match each hit to the corresponding TE consensus length
RM_OSC$TE_length <- width(TE_library)[
  match(
    as.character(RM_OSC$V10),
    names(TE_library)
  )
]


# remove elements that are not present in the TE library
RM_OSC <- RM_OSC[
  !is.na(RM_OSC$TE_length),
  ,
  drop = FALSE
]


# length of each RepeatMasker hit
RM_OSC$range <- (
  as.numeric(as.character(RM_OSC$V7)) -
  as.numeric(as.character(RM_OSC$V6))
)


# fraction of the TE consensus represented by the hit
RM_OSC$percent <- (
  RM_OSC$range /
  RM_OSC$TE_length *
  100
)


# keep reasonably complete TE copies
RM_OSC <- RM_OSC[
  RM_OSC$percent >= OSC_MIN_PERCENT,
  ,
  drop = FALSE
]


# ------------------------------------------------------------------
# Roo insertions
# ------------------------------------------------------------------

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


write.table(
  roo,
  file = file.path(
    OUTPUT_DIR,
    "Roo_insertions_repeatmasker.out"
  ),
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------------
# TE abundance in OSC genome
# ------------------------------------------------------------------

TEs_OSC <- data.frame(
  table(
    RM_OSC$V10
  )
)


TEs_OSC <- TEs_OSC[
  order(
    TEs_OSC$Freq,
    decreasing = TRUE
  ),
]


ggplot(
  TEs_OSC,
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
  xlab(
    "TE name"
  ) +
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
    "OSC genome - TE statistics"
  )


# ------------------------------------------------------------------
# Roo GRanges
# ------------------------------------------------------------------

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


# keep approximately full-length Roo copies
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


# ------------------------------------------------------------------
# All other TEs
# ------------------------------------------------------------------

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


# remove Roo itself
All_TEs <- All_TEs[
  All_TEs$Element != TARGET_TE
]


cat(
  "Number of non-Roo TE hits:",
  length(All_TEs),
  "\n"
)


# ------------------------------------------------------------------
# Other TEs around Roo insertions
# ------------------------------------------------------------------

# make +/- 5 kb regions around Roo insertions
# pmax() prevents coordinates below 1

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


# Roo regions containing another TE
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


# View(
#   as.data.frame(Roo_overlap_regions)
# )


# ------------------------------------------------------------------
# dm6 RepeatMasker data
# ------------------------------------------------------------------

RM_dm6 <- read_repeatmasker(
  DM6_REPEATMASKER_FILE
)


RM_dm6$TE_length <- width(TE_library)[
  match(
    as.character(RM_dm6$V10),
    names(TE_library)
  )
]


RM_dm6 <- RM_dm6[
  !is.na(RM_dm6$TE_length),
  ,
  drop = FALSE
]


RM_dm6$range <- (
  as.numeric(as.character(RM_dm6$V7)) -
  as.numeric(as.character(RM_dm6$V6))
)


RM_dm6$percent <- (
  RM_dm6$range /
  RM_dm6$TE_length *
  100
)


# full-length dm6 TE copies
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


# ------------------------------------------------------------------
# Full-length OSC TE copies
# ------------------------------------------------------------------

# use the same 90% cutoff for OSCs and dm6 for this comparison

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


# ------------------------------------------------------------------
# Compare TE copy numbers between dm6 and OSCs
# ------------------------------------------------------------------

transposon_summary <- data.frame(

  transposon = names(
    TE_library
  ),

  dm6 = NA,

  OSCs = NA,

  stringsAsFactors = FALSE
)


transposon_summary$dm6 <- TEs_dm6$Freq[
  match(
    transposon_summary$transposon,
    TEs_dm6$Var1
  )
]


transposon_summary$OSCs <- TEs_OSC_full$Freq[
  match(
    transposon_summary$transposon,
    TEs_OSC_full$Var1
  )
]


# missing copy numbers = 0
transposon_summary$dm6[
  is.na(transposon_summary$dm6)
] <- 0


transposon_summary$OSCs[
  is.na(transposon_summary$OSCs)
] <- 0


# remove TEs absent from both genomes
transposon_summary <- transposon_summary[
  transposon_summary$dm6 +
    transposon_summary$OSCs != 0,
  ,
  drop = FALSE
]


# long format for plotting
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


ggplot(
  df,
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


# ------------------------------------------------------------------
# Roo sequence composition
# ------------------------------------------------------------------

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


# split Roo sequence into single bases
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


# 50-bp windows: 1-50, 51-100, ...
sequence_df$window <- ceiling(
  sequence_df$position /
    SEQUENCE_WINDOW
)


# calculate A/T/G/C content in each window
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


ggplot(
  ATGC,
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


# ------------------------------------------------------------------
# Extract Roo sequences from the OSC genome
# ------------------------------------------------------------------

genome_OSC <- readDNAStringSet(
  OSC_GENOME_FILE
)


Roo_GR$sequence <- getSeq(
  genome_OSC,
  Roo_GR
)


# View(
#   as.data.frame(Roo_GR)
# )


# ------------------------------------------------------------------
# Simple BED file
# ------------------------------------------------------------------

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


# ------------------------------------------------------------------
# BED6 file
# ------------------------------------------------------------------

# BED uses a 0-based start coordinate,
# so RepeatMasker start is converted with start - 1

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


# ------------------------------------------------------------------
# Roo insertion mapping
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


# convert BAM alignment to a data frame
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


# scanBam() expects BAM-format input
# change filenames here if needed

GFP4_FILE <- file.path(
  INSERTION_ANALYSIS_DIR,
  "GFP4_align_to_163.sam.MAPPED"
)


MAEL4_FILE <- file.path(
  INSERTION_ANALYSIS_DIR,
  "Mael4_align_to_163.sam.MAPPED"
)


# ------------------------------------------------------------------
# GFP4
# ------------------------------------------------------------------

GFP4 <- bam_to_df(
  GFP4_FILE
)


# read abundance is stored after "=" in the read name
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


# total reads for each Roo insertion
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


# ------------------------------------------------------------------
# Mael4
# ------------------------------------------------------------------

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


# ------------------------------------------------------------------
# Roo insertions found in Mael4 but not GFP4
# ------------------------------------------------------------------

Mael4_specific_insertions <- unique(

  Mael4_summary[
    !Mael4_summary$rname %in%
      GFP4_summary$rname,
  ]$rname
)


print(
  Mael4_specific_insertions
)


# ------------------------------------------------------------------
# Save summary tables
# ------------------------------------------------------------------

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
