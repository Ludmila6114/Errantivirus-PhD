# ==============================================================================
# Drosophila ovary protein candidate selection for protein-interaction screening
# ==============================================================================
#
# Purpose:
#
#   1. Identify genes expressed in Drosophila ovaries using RNA-seq TPM values.
#   2. Prioritize expressed proteins associated with:
#        - transmembrane domains
#        - membrane-related GO terms
#        - receptor-related GO terms
#   3. Convert selected genes to UniProt identifiers.
#   4. Generate a UniProt ID list for the complete Drosophila proteome.
#
# Inputs:
#
#   - Salmon gene-level quantification files from biological replicates
#   - BioMart annotation table containing:
#         Gene.name
#         GO.term.name
#         Transmembrane.helices
#   - Drosophila protein FASTA file containing FlyBase and UniProt annotations
#
# IMPORTANT:
#   No local/private paths are included in this script.
#   Replace the paths in the USER SETTINGS section with your own.
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

library(Biostrings)


# ==============================================================================
# 2. USER SETTINGS
# ==============================================================================

# RNA-seq quantification files.
QUANT_FILES <- c(
  "/path/to/replicate1.sf",
  "/path/to/replicate2.sf",
  "/path/to/replicate3.sf"
)

# BioMart annotation table.
GO_ANNOTATION_FILE <- "/path/to/mart_export.csv"

# Drosophila translated proteome FASTA.
PROTEOME_FASTA <- "/path/to/dmel_proteome.fasta"

# Directory for output files.
OUTPUT_DIR <- "/path/to/output"


# Expression threshold.
TPM_THRESHOLD <- 1


# Minimum number of replicates in which a gene must pass the TPM threshold.
#
# 1 = expressed in at least one replicate
# 2 = expressed in at least two replicates
# 3 = expressed in all three replicates
#
# Setting this to 1 reproduces the logic of the original analysis.

MIN_REPLICATES <- 1


# Create output directory if necessary.
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. READ RNA-SEQ QUANTIFICATION FILES
# ==============================================================================

# Salmon quantification files normally contain:
#
#   Name
#   Length
#   EffectiveLength
#   TPM
#   NumReads
#
# This helper extracts only the gene name and TPM columns.

read_quantification <- function(path) {

  if (!file.exists(path)) {
    stop(
      paste(
        "Quantification file not found:",
        path
      )
    )
  }


  x <- read.delim(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  if (!"TPM" %in% colnames(x)) {
    stop(
      paste(
        "No TPM column found in:",
        path
      )
    )
  }


  # The first column contains the gene identifier.
  gene_column <- colnames(x)[1]


  data.frame(
    gene = as.character(x[[gene_column]]),
    TPM = as.numeric(x$TPM),
    stringsAsFactors = FALSE
  )
}


# Read all biological replicates.
quant_data <- lapply(
  QUANT_FILES,
  read_quantification
)


# ==============================================================================
# 4. IDENTIFY GENES EXPRESSED IN THE OVARY
# ==============================================================================

# Determine which genes have TPM >= threshold in each replicate.

expressed_per_replicate <- lapply(

  quant_data,

  function(x) {

    unique(
      x$gene[
        !is.na(x$TPM) &
        x$TPM >= TPM_THRESHOLD
      ]
    )
  }
)


# Count in how many replicates each gene passes the TPM threshold.

expression_counts <- table(
  unlist(
    expressed_per_replicate
  )
)


# Keep genes passing the threshold in the requested number of replicates.

expressed_genes <- names(
  expression_counts[
    expression_counts >= MIN_REPLICATES
  ]
)


cat(
  "Genes passing expression filter:",
  length(expressed_genes),
  "\n"
)


# ==============================================================================
# 5. LOAD BIOMART / GO ANNOTATIONS
# ==============================================================================

go <- read.csv(
  GO_ANNOTATION_FILE,
  stringsAsFactors = FALSE,
  check.names = TRUE
)


# Columns expected from the original BioMart export.

required_columns <- c(
  "Gene.name",
  "GO.term.name",
  "Transmembrane.helices"
)


missing_columns <- setdiff(
  required_columns,
  colnames(go)
)


if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing annotation columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# Keep only genes detected as expressed in the ovary.

go_expressed <- go[
  go$Gene.name %in% expressed_genes,
  ,
  drop = FALSE
]


# ==============================================================================
# 6. IDENTIFY MEMBRANE / TRANSMEMBRANE PROTEINS
# ==============================================================================

# Genes with an annotated transmembrane helix.

genes_TM <- unique(

  go_expressed$Gene.name[
    !is.na(go_expressed$Transmembrane.helices) &
    go_expressed$Transmembrane.helices == "TMhelix"
  ]

)


# Genes containing "membrane" in their GO annotation.

genes_membrane <- unique(

  go_expressed$Gene.name[
    !is.na(go_expressed$GO.term.name) &
    grepl(
      "membrane",
      go_expressed$GO.term.name,
      ignore.case = TRUE
    )
  ]

)


# Main priority set:
#
# expressed in ovary AND
# either transmembrane or membrane-associated.

priority_genes <- unique(
  c(
    genes_TM,
    genes_membrane
  )
)


cat(
  "Expressed membrane/TM candidate genes:",
  length(priority_genes),
  "\n"
)


# ==============================================================================
# 7. IDENTIFY ADDITIONAL RECEPTOR GENES
# ==============================================================================

# Find expressed genes with receptor-related GO annotations.

genes_receptor <- unique(

  go_expressed$Gene.name[
    !is.na(go_expressed$GO.term.name) &
    grepl(
      "receptor",
      go_expressed$GO.term.name,
      ignore.case = TRUE
    )
  ]

)


# Keep receptors that were NOT already included in the membrane/TM list.

additional_receptor_genes <- setdiff(
  genes_receptor,
  priority_genes
)


cat(
  "Additional receptor genes:",
  length(additional_receptor_genes),
  "\n"
)


# ==============================================================================
# 8. LOAD DROSOPHILA PROTEOME
# ==============================================================================

# The FASTA contains protein sequences, therefore AAStringSet is used.

proteome <- readAAStringSet(
  PROTEOME_FASTA
)


protein_headers <- names(
  proteome
)


cat(
  "Protein sequences in proteome:",
  length(proteome),
  "\n"
)


# ==============================================================================
# 9. HELPER FUNCTIONS FOR FASTA HEADER PARSING
# ==============================================================================

# Extract everything after a particular annotation tag.

extract_after_tag <- function(x, tag) {

  vapply(

    strsplit(
      x,
      tag,
      fixed = TRUE
    ),

    function(parts) {

      if (length(parts) >= 2) {

        parts[2]

      } else {

        NA_character_
      }
    },

    character(1)
  )
}


# Extract the first identifier following a particular tag.
#
# Annotation entries in the FASTA headers may be separated by:
#
#   ;
#   ,
#   whitespace

extract_identifier <- function(x, tag) {

  value <- extract_after_tag(
    x,
    tag
  )


  value <- sub(
    "[;,[:space:]].*$",
    "",
    value
  )


  value[
    value == ""
  ] <- NA_character_


  return(
    value
  )
}


# ==============================================================================
# 10. PARSE FLYBASE AND UNIPROT ANNOTATIONS
# ==============================================================================

# Extract FlyBase annotation.

flybase_annotation <- extract_identifier(
  protein_headers,
  "FlyBase_Annotation_IDs:"
)


# Remove protein-isoform suffix.
#
# Example:
#
#   gene-PA -> gene
#   gene-PB -> gene

gene_name <- sub(
  "-P.*$",
  "",
  flybase_annotation
)


# Extract reviewed Swiss-Prot identifiers.

uniprot_swiss <- extract_identifier(
  protein_headers,
  "UniProt/Swiss-Prot:"
)


# Extract TrEMBL identifiers.

uniprot_trembl <- extract_identifier(
  protein_headers,
  "UniProt/TrEMBL:"
)


# Build protein annotation table.

protein_map <- data.frame(

  gene_name = gene_name,

  uniprot_swiss = uniprot_swiss,

  uniprot_trembl = uniprot_trembl,

  stringsAsFactors = FALSE
)


# ==============================================================================
# 11. CREATE A PREFERRED UNIPROT IDENTIFIER
# ==============================================================================

# Prefer Swiss-Prot when available.
#
# If no Swiss-Prot ID exists, use the corresponding TrEMBL ID.

protein_map$preferred_uniprot <- ifelse(

  !is.na(
    protein_map$uniprot_swiss
  ),

  protein_map$uniprot_swiss,

  protein_map$uniprot_trembl

)


# ==============================================================================
# 12. HELPER FUNCTION: GET ALL UNIPROT IDs FOR A GENE SET
# ==============================================================================

get_uniprot_ids <- function(genes) {

  x <- protein_map[
    protein_map$gene_name %in% genes,
    ,
    drop = FALSE
  ]


  # Include both Swiss-Prot and TrEMBL IDs where available.
  ids <- c(
    x$uniprot_swiss,
    x$uniprot_trembl
  )


  # Remove missing and empty identifiers.
  ids <- ids[
    !is.na(ids) &
    ids != ""
  ]


  unique(
    ids
  )
}


# ==============================================================================
# 13. UNIPROT IDS FOR MEMBRANE / TM CANDIDATES
# ==============================================================================

priority_uniprot <- get_uniprot_ids(
  priority_genes
)


priority_output <- data.frame(
  ID = priority_uniprot
)


cat(
  "UniProt IDs for membrane/TM candidates:",
  nrow(priority_output),
  "\n"
)


write.csv(

  priority_output,

  file = file.path(
    OUTPUT_DIR,
    "dmel_membrane_expressed_uniprots.csv"
  ),

  quote = FALSE,

  row.names = FALSE
)


# ==============================================================================
# 14. UNIPROT IDS FOR ADDITIONAL RECEPTORS
# ==============================================================================

receptor_uniprot <- get_uniprot_ids(
  additional_receptor_genes
)


receptor_output <- data.frame(
  ID = receptor_uniprot
)


cat(
  "UniProt IDs for additional receptors:",
  nrow(receptor_output),
  "\n"
)


write.csv(

  receptor_output,

  file = file.path(
    OUTPUT_DIR,
    "dmel_additional_receptors.csv"
  ),

  quote = FALSE,

  row.names = FALSE
)


# ==============================================================================
# 15. UNIPROT IDS FOR THE COMPLETE DROSOPHILA PROTEOME
# ==============================================================================

# For the complete proteome:
#
#   Swiss-Prot is preferred when available.
#   Otherwise, TrEMBL is used.

complete_proteome_ids <- unique(
  protein_map$preferred_uniprot[
    !is.na(
      protein_map$preferred_uniprot
    ) &
    protein_map$preferred_uniprot != ""
  ]
)


complete_proteome_output <- data.frame(
  ID = complete_proteome_ids
)


cat(
  "Proteins with UniProt identifiers:",
  nrow(complete_proteome_output),
  "\n"
)


write.csv(

  complete_proteome_output,

  file = file.path(
    OUTPUT_DIR,
    "filtered_proteome_dmel.csv"
  ),

  quote = FALSE,

  row.names = FALSE
)


# ==============================================================================
# 16. SUMMARY
# ==============================================================================

cat("\n")
cat("Analysis complete\n")
cat("-----------------\n")

cat(
  "Ovary-expressed genes:",
  length(expressed_genes),
  "\n"
)

cat(
  "Membrane/TM candidate genes:",
  length(priority_genes),
  "\n"
)

cat(
  "Additional receptor genes:",
  length(additional_receptor_genes),
  "\n"
)

cat(
  "Membrane/TM UniProt IDs:",
  length(priority_uniprot),
  "\n"
)

cat(
  "Additional receptor UniProt IDs:",
  length(receptor_uniprot),
  "\n"
)

cat(
  "Complete proteome UniProt IDs:",
  length(complete_proteome_ids),
  "\n"
)


# ==============================================================================
# End of script
# ==============================================================================
