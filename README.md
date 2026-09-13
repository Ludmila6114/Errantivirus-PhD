# PhD-project related scripts that I used to analyse individual project-related experiments.

The code covers TED-seq, Tirant insertion analysis, piRNA and small-RNA targeting, RepeatMasker analysis, TE copy number, Roo and Rover sequence analysis, qPCR plotting and protein candidate selection.

## 1. TED-seq processing

`ted_seq.sh` processes paired-end TED-seq data.

The main steps are:

1. check the expected LTR motif and its orientation
2. keep read pairs in which R1 contains the reverse-complement LTR sequence
3. remove the LTR sequence and keep the genomic flanking sequence
4. filter low-quality and low-complexity reads with `fastp`
5. remove PCR duplicates
6. align the remaining reads to the reference genome with `bwa mem`
7. keep high-confidence properly paired alignments with `samtools`
8. generate mapping statistics

Example:

```bash
R1="/path/to/read1.fastq.gz" \
R2="/path/to/read2.fastq.gz" \
GENOME="/path/to/reference.fa" \
SAMPLE="sample_name" \
OUTDIR="./ted_seq_output" \
./ted_seq.sh
```

The pipeline uses:

```text
seqkit
fastp
bwa
samtools
gzip
awk
sort
```

Input paths are provided when the script is run and are not stored in the public version.

## 2. Tirant insertion classification after TED-seq

This R analysis is used after the TED-seq processing step.

It starts from primer/alignment hits and looks at the genomic sequence around each candidate site.

Hits are classified as:

```text
background
potential first LTR
potential second LTR
unknown
```

The analysis:

- extracts 20 bp upstream of each candidate hit
- checks for the sequence expected next to a possible second Tirant LTR
- uses hits in the OSC genome to define background sequences
- searches upstream of possible second LTRs for matching first-LTR candidates
- checks the expected 14-bp sequence downstream of each hit
- separates unique and multi-hit alignments
- counts supporting reads around each candidate site

Possible first and second LTRs are linked using genomic position and strand.

Alignments with MAPQ 0 or an `XA` alternative-alignment tag are treated as multi-hit. Alignments with MAPQ above 0 and no `XA` tag are treated as unique.

## 3. Tirant 25-mer analysis

This analysis splits the Tirant consensus into overlapping 25-bp sequences:

```text
1-25
2-26
3-27
4-28
...
```

The 25-mers are written to FASTA and can then be aligned to a TE library or genome.

The downstream R analysis is used to:

- find Tirant regions that are also similar to other TE families
- mark potentially non-specific 25-mers
- remove these positions from the genome search
- identify Tirant-like regions in the OSC genome
- merge nearby mappings into larger regions
- add genomic flanking sequence
- extract candidate Tirant regions
- save the regions as FASTA
- calculate coverage along the Tirant consensus

The alignment itself is done separately. The R script prepares sequences and works with the alignment results.

## 4. piRNA-seq TE analysis

This script compares antisense piRNA abundance across TE families and Drosophila libraries.

TEs are grouped into related LTR-retrotransposon groups, including:

```text
gypsy
Tabor
mdg3
Osvaldo
bel
copia
```

Current and older datasets are matched by TE name rather than by row order.

The analysis was used to compare piRNA levels across different strains and libraries, with particular attention to **Tirant** and **ZAM**.

The script also contains additional classification of other elements into LINE and DNA transposon groups.

## 5. Tirant piRNA coverage across Drosophila strains

This script compares strand-specific small-RNA/piRNA coverage along the Tirant consensus between different Drosophila strains.

The analysis includes libraries such as:

```text
attp2
attp40
Iso1
ZH11
Chr2 balancer
Chr3 balancer
double balancer
```

For each library, the script:

- selects reads mapping to Tirant
- separates sense and antisense alignments
- recovers read abundance from collapsed read names
- calculates coverage along the consensus
- normalizes using a sample-specific coefficient
- summarizes coverage in 100-bp windows

Sense coverage is plotted above zero and antisense coverage below zero.

The same analysis can also be used for another TE by changing the TE name.

## 6. Small-RNA targeting along TE consensus sequences

This is a more general version of the TE coverage analysis.

It calculates small-RNA targeting profiles for many different TEs and compares several libraries.

Examples of libraries used here are:

```text
Control KD
Panx KD
Piwi KD
Nxf2 KD
```

For every TE, the script:

- selects reads mapping to that TE
- separates sense and antisense reads
- uses collapsed read counts when available
- calculates strand-specific coverage
- normalizes each library
- averages coverage in fixed-size windows
- plots the libraries together

It can make plots for a single TE or generate a multi-page PDF with many TE profiles.

Examples include:

```text
Stalker
ZAM
gypsy
mdg1
roo
Tirant
```

## 7. RepeatMasker and Roo analysis

This script works with RepeatMasker annotations from the OSC genome and dm6.

The analysis includes:

- calculating how much of each TE consensus is covered by a RepeatMasker hit
- identifying Roo insertions
- counting TE copies in the OSC genome
- selecting approximately full-length TE copies
- checking for other TEs around Roo insertions
- comparing full-length TE copy numbers between OSCs and dm6
- examining Roo sequence composition
- extracting Roo sequences from the OSC genome
- creating BED files
- comparing reads mapping to Roo insertion loci

For the initial OSC analysis, a more relaxed cutoff can be used to keep partial TE copies.

For the dm6 versus OSC comparison, the same full-length cutoff is used for both genomes so that the copy numbers are comparable.

The script also divides the Roo consensus into 50-bp windows and calculates local A, T, G and C content.

## 8. Rover sequence analysis

This script uses RepeatMasker annotations to find relatively complete Rover copies.

It:

- matches RepeatMasker hits to the Rover consensus length
- filters for relatively complete Rover copies
- converts the coordinates to `GRanges`
- adds genomic sequence on both sides of each insertion
- extracts the sequences from the genome
- saves the regions as FASTA
- compares gag and LTR coordinates

The TE completeness cutoff and flank size can be changed at the beginning of the script.

## 9. qPCR plotting

This script contains a reusable function for plotting qPCR experiments.

The main input columns are:

```text
TE
treatment
value
```

An optional `experiment` column can be used for faceting.

The script contains examples from experiments involving:

- Mael
- Gtsf1
- Piwi
- TE expression
- rescue constructs

The plotting function can also be reused for other qPCR datasets with the same simple input format.

## 10. Ovary protein candidate selection

This script combines ovary RNA-seq data, BioMart annotations and the Drosophila proteome to select possible candidates for protein-interaction experiments.

The workflow:

1. filters genes using RNA-seq TPM values
2. keeps genes expressed in the requested number of biological replicates
3. finds expressed genes with transmembrane helices
4. finds expressed genes with membrane-related GO terms
5. finds additional receptor-related genes
6. connects FlyBase protein annotations to UniProt IDs
7. exports candidate UniProt lists

For the complete proteome list, Swiss-Prot is preferred when available and TrEMBL is used otherwise.

The number of replicates required to pass the expression filter can be changed in the script.

For example:

```r
MIN_REPLICATES <- 1
```

means that a gene needs to pass the TPM threshold in at least one replicate.

## Input data

Large input files are not included in this repository.

Depending on the analysis, the scripts may require:

```text
FASTQ files
SAM or BAM files
genome FASTA files
TE consensus FASTA files
RepeatMasker output
Salmon quantification files
BioMart annotation tables
small-RNA TE count tables
```

Paths are set at the beginning of the R scripts, for example:

```r
TE_LIBRARY_FILE <- "/path/to/TE_library.fasta"
GENOME_FILE <- "/path/to/genome.fasta"
OUTPUT_DIR <- "/path/to/output"
```

Private and local filesystem paths have been removed from the public versions.

## Main R packages

The scripts use different combinations of:

```text
Biostrings
GenomicRanges
GenomicAlignments
Rsamtools
IRanges
ggplot2
ggbeeswarm
dplyr
tidyr
```

Not every script requires all of them.

Bioconductor packages can be installed with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "Biostrings",
  "GenomicRanges",
  "GenomicAlignments",
  "Rsamtools",
  "IRanges"
))
```

Other R packages can be installed with:

```r
install.packages(c(
  "ggplot2",
  "ggbeeswarm",
  "dplyr",
  "tidyr"
))
```

## Notes

These scripts were developed during an experimental PhD project, so some parts are specific to the original datasets.

Before using them with another dataset, it is worth checking:

- sample names
- TE names
- normalization coefficients
- thresholds
- genome and contig names
- input table structure
- filenames

Some older alignment files have extensions such as `.sam.out` or `.sam.MAPPED`. `Rsamtools::scanBam()` and `GenomicAlignments::readGAlignments()` require BAM-compatible input, so plain-text SAM files need to be converted to BAM first.

Several small-RNA datasets contain collapsed reads where the read count is stored in the read name. The relevant scripts recover this value and use it when calculating coverage.

This repository contains the analysis code, but not all raw sequencing data, genome assemblies or other large reference files used during the project.

## Author

**Liudmila Protsenko**

PhD researcher working on transposable elements and errantivirus biology in *Drosophila melanogaster*.

My background is in Applied Mathematics and Physics, and my PhD work combines molecular biology, genomics and bioinformatics.
