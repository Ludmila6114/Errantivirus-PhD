# Rover sequence analysis
# Find long Rover copies, extract them with flanking sequence,
# and compare gag and LTR coordinates.

library(Biostrings)
library(GenomicRanges)
library(ggplot2)


# Input files
GENOME_FILE <- "/path/to/genome.fasta"
TE_LIBRARY_FILE <- "/path/to/TE_library.fasta"
REPEATMASKER_FILE <- "/path/to/repeatmasker.out"

# Files with gag and LTR coordinates
GAG_FILE <- "/path/to/gag_coordinates.txt"
LTR_FILE <- "/path/to/LTR_coordinates.txt"

# Output folder
OUTPUT_DIR <- "/path/to/output"


# Analysis settings
TARGET_TE <- "rover"
MIN_TE_PERCENT <- 80
FLANK_SIZE <- 500
LEFT_LTR_MAX_POSITION <- 1000


# Create output folder if it does not exist
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# Read genome and TE consensus sequences
genome <- readDNAStringSet(GENOME_FILE)
TE_library <- readDNAStringSet(TE_LIBRARY_FILE)


# Check that the TE is present in the library
if (!TARGET_TE %in% names(TE_library)) {
  stop(paste(TARGET_TE, "was not found in the TE library"))
}


# Read RepeatMasker output
RM <- read.table(
  REPEATMASKER_FILE,
  fill = TRUE,
  stringsAsFactors = FALSE
)


# Remove RepeatMasker header lines
RM <- RM[3:nrow(RM), , drop = FALSE]
RM <- na.omit(RM)


# Match each RepeatMasker hit to the corresponding TE consensus length
RM$TE_width <- width(TE_library)[
  match(
    as.character(RM$V10),
    names(TE_library)
  )
]


# Remove hits that are not present in the TE library
RM <- RM[
  !is.na(RM$TE_width),
  ,
  drop = FALSE
]


# Calculate how much of the TE consensus is covered by each hit
RM$hit_width <- (
  as.numeric(as.character(RM$V7)) -
  as.numeric(as.character(RM$V6)) +
  1
)

RM$percent <- (
  RM$hit_width /
  RM$TE_width *
  100
)


# Keep relatively complete TE copies
RM_filtered <- RM[
  RM$percent > MIN_TE_PERCENT,
  ,
  drop = FALSE
]


# Select Rover
target_hits <- RM_filtered[
  RM_filtered$V10 == TARGET_TE,
  ,
  drop = FALSE
]


cat(
  "Number of",
  TARGET_TE,
  "copies:",
  nrow(target_hits),
  "\n"
)


# Convert Rover coordinates to GRanges
TE_GR <- GRanges(

  seqnames = as.character(
    target_hits$V5
  ),

  ranges = IRanges(
    start = as.numeric(as.character(target_hits$V6)),
    end = as.numeric(as.character(target_hits$V7))
  ),

  strand = ifelse(
    as.character(target_hits$V9) == "+",
    "+",
    "-"
  ),

  Element = as.character(
    target_hits$V10
  ),

  percent = target_hits$percent
)


# Check that chromosome names match the genome FASTA
missing_contigs <- setdiff(
  unique(as.character(seqnames(TE_GR))),
  names(genome)
)

if (length(missing_contigs) > 0) {
  stop(
    paste(
      "Contigs missing from genome FASTA:",
      paste(missing_contigs, collapse = ", ")
    )
  )
}


# Add flanking sequence around each insertion
genome_lengths <- setNames(
  width(genome),
  names(genome)
)

new_start <- pmax(
  1,
  start(TE_GR) - FLANK_SIZE
)

new_end <- pmin(
  genome_lengths[
    as.character(seqnames(TE_GR))
  ],
  end(TE_GR) + FLANK_SIZE
)

ranges(TE_GR) <- IRanges(
  start = new_start,
  end = new_end
)


# Extract sequences from the genome
TE_GR$sequence <- getSeq(
  genome,
  TE_GR
)


# Convert extracted sequences to FASTA
TE_fasta <- DNAStringSet(
  TE_GR$sequence
)


# Add useful names to the FASTA entries
names(TE_fasta) <- paste0(
  as.character(seqnames(TE_GR)),
  ":",
  start(TE_GR),
  "-",
  end(TE_GR),
  ":",
  as.character(strand(TE_GR))
)


# Save FASTA
FASTA_OUTPUT <- file.path(
  OUTPUT_DIR,
  paste0(
    TARGET_TE,
    "_sequences_with_",
    FLANK_SIZE,
    "bp_flanks.fasta"
  )
)

writeXStringSet(
  TE_fasta,
  FASTA_OUTPUT
)


# ------------------------------------------------------------------
# gag and LTR analysis
# ------------------------------------------------------------------

# Read gag and LTR coordinates
gag <- read.table(
  GAG_FILE,
  stringsAsFactors = FALSE
)

LTR <- read.table(
  LTR_FILE,
  stringsAsFactors = FALSE
)


# If needed, select only specific gag entries
# gag <- gag[2:6, , drop = FALSE]


# Keep LTRs belonging to the same sequences as gag
LTR <- LTR[
  LTR$V1 %in% gag$V1,
  ,
  drop = FALSE
]


# Keep LTRs close to the beginning of the sequence
# These should correspond to the left / 5' LTR
LTR_left <- LTR[
  as.numeric(LTR$V3) < LEFT_LTR_MAX_POSITION,
  ,
  drop = FALSE
]


# Keep one LTR per sequence
LTR_left <- LTR_left[
  !duplicated(LTR_left$V1),
  ,
  drop = FALSE
]


# Match LTR end coordinates to gag
gag$LTR_end <- LTR_left$V4[
  match(
    gag$V1,
    LTR_left$V1
  )
]


# Remove sequences without a matching LTR
gag <- gag[
  !is.na(gag$LTR_end),
  ,
  drop = FALSE
]


# Calculate distance between gag and the left LTR
gag$gag_LTR_distance <- (
  as.numeric(gag$V4) -
  as.numeric(gag$LTR_end)
)


# Look at the result
print(
  gag[
    ,
    c(
      "V1",
      "V4",
      "LTR_end",
      "gag_LTR_distance"
    )
  ]
)


# Plot the distribution of gag-LTR distances
distance_plot <- ggplot(
  gag,
  aes(
    x = gag_LTR_distance
  )
) +
  geom_density(
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    title = paste(TARGET_TE, "gag-LTR distance"),
    x = "Distance between LTR end and gag coordinate (bp)",
    y = "Density"
  ) +
  theme(
    text = element_text(size = 14)
  )


print(distance_plot)


# Save plot
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    paste0(
      TARGET_TE,
      "_gag_LTR_distance.pdf"
    )
  ),
  plot = distance_plot,
  width = 7,
  height = 5,
  units = "in"
)
