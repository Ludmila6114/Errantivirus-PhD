# Drosophila ovary protein candidate selection
#
# Find genes expressed in ovary, prioritize membrane/TM/receptor proteins,
# and convert the selected genes to UniProt IDs.


library(Biostrings)

# Input files
QUANT_FILES <- c(
  "/path/to/replicate1.sf",
  "/path/to/replicate2.sf",
  "/path/to/replicate3.sf"
)

GO_ANNOTATION_FILE <- "/path/to/mart_export.csv"
PROTEOME_FASTA <- "/path/to/dmel_proteome.fasta"

OUTPUT_DIR <- "/path/to/output"


# expression cutoff
TPM_THRESHOLD <- 1

# number of replicates in which a gene must pass the TPM cutoff
# 1 = at least one replicate
# 2 = at least two replicates
# 3 = all three replicates
MIN_REPLICATES <- 1


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# RNA-seq data
# ------------------------------------------------------------------

# read Salmon quantification file and keep gene name + TPM
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


  gene_column <- colnames(x)[1]


  data.frame(
    gene = as.character(x[[gene_column]]),
    TPM = as.numeric(x$TPM),
    stringsAsFactors = FALSE
  )
}


# read all replicates
quant_data <- lapply(
  QUANT_FILES,
  read_quantification
)


# genes above the TPM cutoff in each replicate
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


# count in how many replicates each gene is expressed
expression_counts <- table(
  unlist(
    expressed_per_replicate
  )
)


# keep genes passing the cutoff in enough replicates
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


# GO / BioMart annotation
# ------------------------------------------------------------------

go <- read.csv(
  GO_ANNOTATION_FILE,
  stringsAsFactors = FALSE,
  check.names = TRUE
)


# columns expected in the BioMart table
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


# keep only genes expressed in ovary
go_expressed <- go[
  go$Gene.name %in% expressed_genes,
  ,
  drop = FALSE
]


# Membrane and transmembrane proteins
# ------------------------------------------------------------------

# genes with a predicted/annotated transmembrane helix
genes_TM <- unique(

  go_expressed$Gene.name[
    !is.na(go_expressed$Transmembrane.helices) &
      go_expressed$Transmembrane.helices == "TMhelix"
  ]
)


# genes with "membrane" in a GO term
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


# main set of membrane/TM candidates
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


# Additional receptor genes
# ------------------------------------------------------------------

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


# keep receptors not already in the membrane/TM list
additional_receptor_genes <- setdiff(
  genes_receptor,
  priority_genes
)


cat(
  "Additional receptor genes:",
  length(additional_receptor_genes),
  "\n"
)


# ------------------------------------------------------------------
# Drosophila proteome
# ------------------------------------------------------------------

# protein FASTA, so use AAStringSet
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


# Parse FASTA headers
# ------------------------------------------------------------------

# extract everything after a tag
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


# get the first identifier after a tag
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


# FlyBase annotation
flybase_annotation <- extract_identifier(
  protein_headers,
  "FlyBase_Annotation_IDs:"
)


# remove protein isoform suffix
# e.g. gene-PA -> gene
gene_name <- sub(
  "-P.*$",
  "",
  flybase_annotation
)


# UniProt IDs
uniprot_swiss <- extract_identifier(
  protein_headers,
  "UniProt/Swiss-Prot:"
)

uniprot_trembl <- extract_identifier(
  protein_headers,
  "UniProt/TrEMBL:"
)


protein_map <- data.frame(

  gene_name = gene_name,

  uniprot_swiss = uniprot_swiss,

  uniprot_trembl = uniprot_trembl,

  stringsAsFactors = FALSE
)


# prefer Swiss-Prot if available, otherwise use TrEMBL
protein_map$preferred_uniprot <- ifelse(

  !is.na(
    protein_map$uniprot_swiss
  ),

  protein_map$uniprot_swiss,

  protein_map$uniprot_trembl
)

# Get UniProt IDs for a set of genes
# ------------------------------------------------------------------

get_uniprot_ids <- function(genes) {

  x <- protein_map[
    protein_map$gene_name %in% genes,
    ,
    drop = FALSE
  ]


  # keep both Swiss-Prot and TrEMBL IDs if available
  ids <- c(
    x$uniprot_swiss,
    x$uniprot_trembl
  )


  ids <- ids[
    !is.na(ids) &
      ids != ""
  ]


  unique(
    ids
  )
}


# Membrane / TM candidates
# ------------------------------------------------------------------

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

# Additional receptor candidates
# ------------------------------------------------------------------

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



# one preferred UniProt ID per protein:
# Swiss-Prot first, otherwise TrEMBL

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



cat("\n")
cat("Analysis complete\n")

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
