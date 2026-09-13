# TED-seq follow-up: Tirant insertion classification
#
# Classify candidate Tirant/LTR hits, remove background hits,
# identify possible first/second LTR pairs, and count supporting
# reads from the iso1 and gag libraries.

library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(Biostrings)


# Input files

# Exact primer alignments
ISO1_PRIMER_ALIGNMENT <- "/path/to/primer_to_iso1.sam"
OSC_PRIMER_ALIGNMENT <- "/path/to/primer_to_OSC.sam"

# Reference genomes
ISO1_GENOME_FASTA <- "/path/to/iso1_genome.fa"
OSC_GENOME_FASTA <- "/path/to/OSC_genome.fa"

# TED-seq libraries
ISO1_LIBRARY_ALIGNMENT <- "/path/to/iso1_library.sam"
GAG_LIBRARY_ALIGNMENT <- "/path/to/gag_library.sam"

# Output folder
OUTPUT_DIR <- "/path/to/output"


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# Sequences used for classification

SECOND_LTR_UPSTREAM_SEQ <- "CTTCGAACTGGGGGGGGAGG"
EXPECTED_DOWNSTREAM_SEQ <- "CCCCACGCCTCTAA"

UPSTREAM_WIDTH <- 20
DOWNSTREAM_WIDTH <- 14

# Search region for matching first and second LTRs
LTR_PAIR_DISTANCE <- 8200

# Region upstream of each insertion used for read counting
READ_FLANK_SIZE <- 150


# Convert SAM to BAM if needed

prepare_bam <- function(path, output_name) {

  if (!file.exists(path)) {

    stop(
      paste(
        "Alignment file not found:",
        basename(path)
      )
    )
  }


  # Use BAM directly if the input is already BAM
  if (grepl(
    "\\.bam$",
    path,
    ignore.case = TRUE
  )) {

    return(path)
  }


  bam_prefix <- file.path(
    OUTPUT_DIR,
    output_name
  )


  asBam(
    path,
    destination = bam_prefix,
    overwrite = TRUE
  )


  paste0(
    bam_prefix,
    ".bam"
  )
}


# Read exact primer hits in the iso1 genome

iso1_primer_bam <- prepare_bam(
  ISO1_PRIMER_ALIGNMENT,
  "primer_to_iso1"
)


aln <- readGAlignments(
  iso1_primer_bam,
  use.names = TRUE
)


gr <- granges(
  aln
)


# Specific contigs can be excluded here if needed
#
# gr <- gr[
#   !as.character(seqnames(gr)) %in%
#     c("contig1", "contig2")
# ]


# Table used to classify the hits

hit_summary <- data.frame(

  name = paste0(
    seqnames(gr),
    ":",
    start(gr),
    "-",
    end(gr)
  ),

  index = seq_along(
    gr
  ),

  comment = NA_character_,

  stringsAsFactors = FALSE
)


# Extract 20 bp upstream of each primer hit

iso1_fasta <- FaFile(
  ISO1_GENOME_FASTA
)


gr_upstream20 <- promoters(
  gr,
  upstream = UPSTREAM_WIDTH,
  downstream = 0
)


seq_upstream20 <- getSeq(
  iso1_fasta,
  gr_upstream20
)


hit_summary$upstream20 <- as.character(
  seq_upstream20
)


# This upstream sequence is expected for a possible second Tirant LTR

hit_summary$comment <- ifelse(

  hit_summary$upstream20 ==
    SECOND_LTR_UPSTREAM_SEQ,

  "potential second LTR",

  NA_character_
)


table(
  hit_summary$comment,
  useNA = "ifany"
)


# Find background primer hits in the OSC genome

OSC_primer_bam <- prepare_bam(
  OSC_PRIMER_ALIGNMENT,
  "primer_to_OSC"
)


OSC_aln <- readGAlignments(
  OSC_primer_bam,
  use.names = TRUE
)


OSCgr <- granges(
  OSC_aln
)


OSC_fasta <- FaFile(
  OSC_GENOME_FASTA
)


OSC_upstream20 <- promoters(
  OSCgr,
  upstream = UPSTREAM_WIDTH,
  downstream = 0
)


OSC_seq_upstream20 <- getSeq(
  OSC_fasta,
  OSC_upstream20
)


background_seqs <- as.character(
  OSC_seq_upstream20
)


# Look at similarity between the OSC background sequences

# Hamming distance = number of different bases between two sequences

if (length(background_seqs) > 1) {

  mat <- do.call(
    rbind,
    strsplit(
      background_seqs,
      ""
    )
  )


  ham <- outer(

    seq_len(
      nrow(mat)
    ),

    seq_len(
      nrow(mat)
    ),

    Vectorize(
      function(i, j) {

        sum(
          mat[i, ] !=
            mat[j, ]
        )
      }
    )
  )


  rownames(ham) <- background_seqs
  colnames(ham) <- background_seqs


  d <- as.dist(
    ham
  )


  hc <- hclust(
    d,
    method = "average"
  )


  plot(
    hc,
    labels = background_seqs,
    main = "20-bp upstream sequences: OSC genome background hits",
    xlab = "",
    sub = ""
  )
}


# Label hits also found in the OSC genome as background

hit_summary$comment <- ifelse(

  hit_summary$upstream20 %in%
    background_seqs,

  "background",

  hit_summary$comment
)


# Everything else is unknown for now

hit_summary$comment[
  is.na(hit_summary$comment)
] <- "unknown"


table(
  hit_summary$comment
)


# Add the classification to the genomic ranges

gr$index <- hit_summary$index
gr$comment <- hit_summary$comment


# Separate possible second LTRs and unknown hits

gr_secondLTR <- gr[
  gr$comment == "potential second LTR"
]


gr_unknown <- gr[
  gr$comment == "unknown"
]


# Look 8.2 kb upstream from each possible second LTR

gr_secondLTR_search <- promoters(
  gr_secondLTR,
  upstream = LTR_PAIR_DISTANCE,
  downstream = 0
)


# Find unknown hits within this region on the same strand

hits <- findOverlaps(
  query = gr_secondLTR_search,
  subject = gr_unknown,
  ignore.strand = FALSE
)


if (length(hits) > 0) {

  second_ltr_hit <- queryHits(
    hits
  )

  unknown_hit <- subjectHits(
    hits
  )


  # Give possible first and second LTRs the same index

  gr_unknown$index[
    unknown_hit
  ] <- gr_secondLTR$index[
    second_ltr_hit
  ]
}


# Put the hits back together

gr_known <- gr[
  gr$comment != "unknown"
]


gr <- c(
  gr_known,
  gr_unknown
)


# Find indices occurring more than once

idx <- as.vector(
  gr$index
)


duplicated_index <- (
  duplicated(idx) |
    duplicated(
      idx,
      fromLast = TRUE
    )
)


# Unknown hits sharing an index with a second LTR
# are labelled as possible first LTRs

rows_to_change <- (
  duplicated_index &
    !is.na(gr$comment) &
    gr$comment == "unknown"
)


gr$comment[
  rows_to_change
] <- "potential first LTR"


table(
  gr$comment
)


# Save the initial classification

write.csv(
  data.frame(gr),
  file = file.path(
    OUTPUT_DIR,
    "classified_Tirant_hits.csv"
  ),
  quote = FALSE,
  row.names = FALSE
)


# Check the sequence downstream of each hit

# Extract 14 bp downstream in strand-aware orientation

gr_down14 <- flank(
  gr,
  width = DOWNSTREAM_WIDTH,
  start = FALSE,
  ignore.strand = FALSE
)


seq_down14 <- getSeq(
  iso1_fasta,
  gr_down14
)


gr$downstream14 <- as.character(
  seq_down14
)


# Check whether the expected sequence is present

gr$expect_to_find <- ifelse(

  gr$downstream14 ==
    EXPECTED_DOWNSTREAM_SEQ,

  "should find",

  "PBS is not correct"
)


table(
  gr$expect_to_find
)


# Read one TED-seq alignment library

# Properly paired primary alignments are kept.
#
# unique:
# MAPQ > 0 and no XA tag
#
# multi-hit:
# MAPQ = 0 or an XA alternative-alignment tag is present

read_library <- function(
    path,
    output_name
) {

  bam_file <- prepare_bam(
    path,
    output_name
  )


  param <- ScanBamParam(

    flag = scanBamFlag(
      isPaired = TRUE,
      isProperPair = TRUE,
      isUnmappedQuery = FALSE,
      hasUnmappedMate = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE
    ),

    what = c(
      "qname",
      "flag",
      "mapq"
    ),

    tag = "XA"
  )


  ga <- readGAlignments(
    bam_file,
    use.names = TRUE,
    param = param
  )


  gr_lib <- granges(
    ga,
    use.names = TRUE,
    use.mcols = TRUE
  )


  gr_lib$qname <- names(
    ga
  )


  # Check that unwanted alignment types were removed

  flag <- as.integer(
    gr_lib$flag
  )


  cat(
    output_name,
    "- unmapped:",
    sum(
      bitwAnd(flag, 4L) != 0L
    ),
    "\n"
  )


  cat(
    output_name,
    "- mate unmapped:",
    sum(
      bitwAnd(flag, 8L) != 0L
    ),
    "\n"
  )


  cat(
    output_name,
    "- secondary:",
    sum(
      bitwAnd(flag, 256L) != 0L
    ),
    "\n"
  )


  cat(
    output_name,
    "- supplementary:",
    sum(
      bitwAnd(flag, 2048L) != 0L
    ),
    "\n"
  )


  mapq <- as.integer(
    gr_lib$mapq
  )


  xa <- as.character(
    gr_lib$XA
  )


  has_XA <- (
    !is.na(xa) &
      nzchar(xa)
  )


  gr_lib$comment <- "unknown"


  gr_lib$comment[
    has_XA |
      (!is.na(mapq) & mapq == 0L)
  ] <- "multi-hit"


  gr_lib$comment[
    !has_XA &
      !is.na(mapq) &
      mapq > 0L
  ] <- "unique"


  cat(
    "\n",
    output_name,
    " alignment classes:\n",
    sep = ""
  )


  print(
    table(
      gr_lib$comment,
      useNA = "ifany"
    )
  )


  return(
    gr_lib
  )
}


# Make a 150-bp upstream region for each classified hit

gr_flank <- promoters(
  gr,
  upstream = READ_FLANK_SIZE,
  downstream = 0
)


# Count unique and multi-hit reads overlapping each region

count_library_overlaps <- function(
    target_gr,
    read_gr,
    prefix
) {

  is_unique <- (
    !is.na(read_gr$comment) &
      read_gr$comment == "unique"
  )


  is_multihit <- (
    !is.na(read_gr$comment) &
      read_gr$comment == "multi-hit"
  )


  # Unique alignments

  hits_unique <- findOverlaps(
    query = read_gr[is_unique],
    subject = target_gr,
    ignore.strand = TRUE
  )


  target_gr[
    [
      paste0(
        prefix,
        "_unique_intersections"
      )
    ]
  ] <- tabulate(
    subjectHits(
      hits_unique
    ),
    nbins = length(
      target_gr
    )
  )


  # Multi-hit alignments

  hits_multihit <- findOverlaps(
    query = read_gr[is_multihit],
    subject = target_gr,
    ignore.strand = TRUE
  )


  target_gr[
    [
      paste0(
        prefix,
        "_multihit_intersections"
      )
    ]
  ] <- tabulate(
    subjectHits(
      hits_multihit
    ),
    nbins = length(
      target_gr
    )
  )


  return(
    target_gr
  )
}


# iso1 library

gr_iso1 <- read_library(
  ISO1_LIBRARY_ALIGNMENT,
  "iso1_library"
)


gr_flank <- count_library_overlaps(
  target_gr = gr_flank,
  read_gr = gr_iso1,
  prefix = "iso1"
)


# gag library

gr_gag <- read_library(
  GAG_LIBRARY_ALIGNMENT,
  "gag_library"
)


gr_flank <- count_library_overlaps(
  target_gr = gr_flank,
  read_gr = gr_gag,
  prefix = "gag"
)


# Check total numbers of overlapping alignments

cat(
  "\niso1 multi-hit overlaps:",
  sum(
    gr_flank$iso1_multihit_intersections
  ),
  "\n"
)


cat(
  "iso1 unique overlaps:",
  sum(
    gr_flank$iso1_unique_intersections
  ),
  "\n"
)


cat(
  "gag multi-hit overlaps:",
  sum(
    gr_flank$gag_multihit_intersections
  ),
  "\n"
)


cat(
  "gag unique overlaps:",
  sum(
    gr_flank$gag_unique_intersections
  ),
  "\n"
)


# Save final statistics

write.csv(
  data.frame(gr_flank),
  file = file.path(
    OUTPUT_DIR,
    "final_Tirant_insertion_statistics.csv"
  ),
  quote = FALSE,
  row.names = FALSE
)
