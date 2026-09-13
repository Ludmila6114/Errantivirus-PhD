#!/usr/bin/env bash

# TE/LTR paired-end processing pipeline
#
# Input paths are provided at runtime so they do not need to be stored here.
#
# Example:
#
# R1="/path/to/read1.fastq.gz" \
# R2="/path/to/read2.fastq.gz" \
# GENOME="/path/to/reference.fa" \
# SAMPLE="sample_name" \
# OUTDIR="./ted_seq_output" \
# ./ted_seq_pipeline.sh
#
# Required programs:
# seqkit, fastp, bwa, samtools, gzip, awk, sort
#
# Optional variables:
# SAMPLE            output prefix (default: sample)
# OUTDIR            output directory (default: ./ted_seq_output)
# MOTIF             LTR sequence
# BWA_THREADS       bwa threads (default: 12)
# SAMTOOLS_THREADS  samtools threads (default: 8)
# FASTP_THREADS     fastp threads (default: 8)
#
# SEQKIT, FASTP, BWA and SAMTOOLS can also be changed if the programs
# are installed under different executable names.
#
# The pipeline assumes that R1 and R2 have matching read IDs.
# The LTR sequence below is the sequence used in the original analysis.

set -euo pipefail

# Do not use shell tracing here because it would print input paths.
set +x


# Input files
R1="${R1:-}"
R2="${R2:-}"
GENOME="${GENOME:-}"

SAMPLE="${SAMPLE:-sample}"
OUTDIR="${OUTDIR:-./ted_seq_output}"


# LTR motif
MOTIF="${MOTIF:-AGTTACCACCCCACCCCCTAAACCCCCACGCCTCTAAACAAAT}"


# Threads
BWA_THREADS="${BWA_THREADS:-12}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-8}"
FASTP_THREADS="${FASTP_THREADS:-8}"


# Program names
SEQKIT="${SEQKIT:-seqkit}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"


# Check that required input variables were provided
if [[ -z "$R1" || -z "$R2" || -z "$GENOME" ]]; then

    cat >&2 <<'USAGE'

ERROR: R1, R2 and GENOME must be provided.

Example:

  R1="/path/to/read1.fastq.gz" \
  R2="/path/to/read2.fastq.gz" \
  GENOME="/path/to/reference.fa" \
  SAMPLE="sample_name" \
  OUTDIR="./ted_seq_output" \
  ./ted_seq_pipeline.sh

USAGE

    exit 1
fi


# Check that input files exist
for input_file in "$R1" "$R2" "$GENOME"; do

    if [[ ! -f "$input_file" ]]; then

        echo "ERROR: Required input file was not found." >&2
        echo "Check the path supplied through R1, R2 or GENOME." >&2

        exit 1
    fi

done


# Check required programs
for cmd in "$SEQKIT" "$FASTP" "$BWA" "$SAMTOOLS" gzip awk sort; do

    if ! command -v "$cmd" >/dev/null 2>&1; then

        echo "ERROR: Required command is not available: $cmd" >&2

        exit 1
    fi

done


# Output folders
FILTER_DIR="$OUTDIR/01_ltr_filter"
TRIM_DIR="$OUTDIR/02_trimmed_flanks"
QC_DIR="$OUTDIR/03_fastp"
ALIGN_DIR="$OUTDIR/04_alignment"
STATS_DIR="$OUTDIR/05_stats"

mkdir -p \
    "$FILTER_DIR" \
    "$TRIM_DIR" \
    "$QC_DIR" \
    "$ALIGN_DIR" \
    "$STATS_DIR"


# Output files
KEPT_R1="$FILTER_DIR/${SAMPLE}.kept_R1.fastq.gz"
KEPT_R2="$FILTER_DIR/${SAMPLE}.kept_R2.fastq.gz"
KEPT_IDS="$FILTER_DIR/${SAMPLE}.kept_ids.txt"

TRIMMED_R1="$TRIM_DIR/${SAMPLE}.trimmed_R1.fastq.gz"
TRIMMED_R2="$TRIM_DIR/${SAMPLE}.trimmed_R2.fastq.gz"

CLEAN_R1="$QC_DIR/${SAMPLE}.clean_R1.fastq.gz"
CLEAN_R2="$QC_DIR/${SAMPLE}.clean_R2.fastq.gz"

DEDUP_R1="$QC_DIR/${SAMPLE}.dedup_R1.fastq.gz"
DEDUP_R2="$QC_DIR/${SAMPLE}.dedup_R2.fastq.gz"

FASTP_FILTER_HTML="$QC_DIR/${SAMPLE}.fastp_filter.html"
FASTP_FILTER_JSON="$QC_DIR/${SAMPLE}.fastp_filter.json"

FASTP_DEDUP_HTML="$QC_DIR/${SAMPLE}.fastp_dedup.html"
FASTP_DEDUP_JSON="$QC_DIR/${SAMPLE}.fastp_dedup.json"

SAM_FILE="$ALIGN_DIR/${SAMPLE}.sam"
BAM_FILE="$ALIGN_DIR/${SAMPLE}.bam"
STRICT_BAM="$ALIGN_DIR/${SAMPLE}.strict.bam"
SORTED_BAM="$ALIGN_DIR/${SAMPLE}.strict.sorted.bam"

ORIENTATION_REPORT="$STATS_DIR/${SAMPLE}.ltr_orientation.tsv"
FLAGSTAT_REPORT="$STATS_DIR/${SAMPLE}.flagstat.txt"
CHROMOSOME_COUNTS="$STATS_DIR/${SAMPLE}.chromosome_counts.tsv"


# Reverse complement of the LTR motif
RC="$({
    printf '>motif\n%s\n' "$MOTIF"
} | "$SEQKIT" seq -r -p 2>/dev/null \
  | awk '!/^>/{printf "%s",$0} END{print ""}')"


if [[ -z "$RC" ]]; then

    echo "ERROR: Could not calculate reverse complement of the LTR motif." >&2

    exit 1
fi


echo "Starting TE/LTR pipeline for sample: $SAMPLE"
echo "Output directory: $OUTDIR"



# ------------------------------------------------------------------------------
# Step 1: check LTR orientation
# ------------------------------------------------------------------------------

# In the original analysis the useful reads contained the LTR in reverse
# complement orientation in R1.
#
# Count both orientations in R1 and R2 first as a quick check.

echo "[Step 1/6] Checking exact LTR motif orientation..."


count_exact_hits() {

    local fastq="$1"
    local pattern="$2"

    "$SEQKIT" grep \
        -s \
        -P \
        -p "$pattern" \
        "$fastq" \
        2>/dev/null \
        | "$SEQKIT" stats -T 2>/dev/null \
        | awk 'NR==2{print $4; found=1} END{if(!found) print 0}'
}


R1_FORWARD_COUNT="$(count_exact_hits "$R1" "$MOTIF")"
R1_RC_COUNT="$(count_exact_hits "$R1" "$RC")"

R2_FORWARD_COUNT="$(count_exact_hits "$R2" "$MOTIF")"
R2_RC_COUNT="$(count_exact_hits "$R2" "$RC")"


{
    printf 'read\torientation\texact_match_count\n'

    printf 'R1\tforward\t%s\n' "$R1_FORWARD_COUNT"
    printf 'R1\treverse_complement\t%s\n' "$R1_RC_COUNT"

    printf 'R2\tforward\t%s\n' "$R2_FORWARD_COUNT"
    printf 'R2\treverse_complement\t%s\n' "$R2_RC_COUNT"

} > "$ORIENTATION_REPORT"


if [[ "$R1_RC_COUNT" =~ ^[0-9]+$ ]] && (( R1_RC_COUNT == 0 )); then

    echo "WARNING: No exact reverse-complement motif matches were found in R1." >&2
    echo "Check the motif and orientation report before interpreting the result." >&2

fi



# ------------------------------------------------------------------------------
# Step 2: select read pairs containing the LTR
# ------------------------------------------------------------------------------

# Keep R1 reads containing the exact reverse-complement LTR.
# Then use their read IDs to retrieve the matching R2 reads.

echo "[Step 2/6] Selecting R1 reads containing the exact RC LTR motif..."


"$SEQKIT" grep \
    -s \
    -P \
    -p "$RC" \
    "$R1" \
    -o "$KEPT_R1"


"$SEQKIT" seq \
    -n \
    -i \
    "$KEPT_R1" \
    > "$KEPT_IDS"


"$SEQKIT" grep \
    -f "$KEPT_IDS" \
    "$R2" \
    -o "$KEPT_R2"


KEPT_COUNT="$(
    "$SEQKIT" stats -T "$KEPT_R1" 2>/dev/null \
    | awk 'NR==2{print $4}'
)"


if [[ -z "$KEPT_COUNT" || "$KEPT_COUNT" == "0" ]]; then

    echo "ERROR: No read pairs remained after exact LTR filtering." >&2

    exit 1
fi


echo "  Retained read pairs: $KEPT_COUNT"



# ------------------------------------------------------------------------------
# Step 3: remove the LTR and keep genomic flanking sequence
# ------------------------------------------------------------------------------

# R2:
# keep the sequence before the forward LTR.
#
# R1:
# remove everything through the reverse-complement LTR and keep the sequence
# after it.
#
# Sequence and quality strings are trimmed together.

echo "[Step 3/6] Trimming LTR sequence and retaining genomic flanks..."


# Trim R2
gzip -cd "$KEPT_R2" \
    | awk -v m="$MOTIF" '

        NR%4==1 {h=$0}
        NR%4==2 {s=$0}
        NR%4==3 {p=$0}

        NR%4==0 {

            q=$0
            i=index(s,m)

            if (i>0) {

                s=substr(s,1,i-1)
                q=substr(q,1,i-1)

            }

            print h
            print s
            print p
            print q
        }

    ' \
    | gzip -c \
    > "$TRIMMED_R2"


# Trim R1
gzip -cd "$KEPT_R1" \
    | awk -v m="$RC" '

        NR%4==1 {h=$0}
        NR%4==2 {s=$0}
        NR%4==3 {p=$0}

        NR%4==0 {

            q=$0
            i=index(s,m)

            if (i>0) {

                cut=i+length(m)-1

                s=substr(s,cut+1)
                q=substr(q,cut+1)

            }

            print h
            print s
            print p
            print q
        }

    ' \
    | gzip -c \
    > "$TRIMMED_R1"



# ------------------------------------------------------------------------------
# Step 4: fastp filtering
# ------------------------------------------------------------------------------

# Remove low-complexity, low-quality and very short reads.
#
# Parameters from the original analysis:
#   complexity threshold = 30
#   minimum base quality = Q20
#   maximum unqualified bases = 40%
#   maximum N bases = 5
#   minimum read length = 10
#
# Adapter trimming is disabled.

echo "[Step 4/6] Running fastp quality and low-complexity filtering..."


"$FASTP" \
    -i "$TRIMMED_R1" \
    -I "$TRIMMED_R2" \
    -o "$CLEAN_R1" \
    -O "$CLEAN_R2" \
    --thread "$FASTP_THREADS" \
    --disable_adapter_trimming \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --length_required 10 \
    --html "$FASTP_FILTER_HTML" \
    --json "$FASTP_FILTER_JSON"



# ------------------------------------------------------------------------------
# Step 5: remove PCR duplicates
# ------------------------------------------------------------------------------

# Run fastp again for deduplication only.
# Quality and length filtering were already done in the previous step.

echo "[Step 5/6] Removing PCR duplicates with fastp..."


"$FASTP" \
    -i "$CLEAN_R1" \
    -I "$CLEAN_R2" \
    -o "$DEDUP_R1" \
    -O "$DEDUP_R2" \
    --thread "$FASTP_THREADS" \
    --dedup \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    --html "$FASTP_DEDUP_HTML" \
    --json "$FASTP_DEDUP_JSON"



# ------------------------------------------------------------------------------
# Step 6: genome alignment
# ------------------------------------------------------------------------------

# Align genomic flanks with bwa mem.
#
# bwa:
#   -T 30     minimum alignment score
#
# samtools filtering:
#   -f 2      keep properly paired reads
#   -F 2308   remove unmapped, secondary and supplementary alignments
#   -q 30     MAPQ >= 30

echo "[Step 6/6] Aligning retained genomic flanks to the reference genome..."


# Make BWA index if it is not already present.
# Index files are written next to the genome FASTA.
if [[ ! -f "${GENOME}.bwt" || ! -f "${GENOME}.sa" ]]; then

    echo "  BWA index not detected; creating it now..."

    "$BWA" index "$GENOME"

fi


# Alignment
"$BWA" mem \
    -t "$BWA_THREADS" \
    -T 30 \
    "$GENOME" \
    "$DEDUP_R1" \
    "$DEDUP_R2" \
    > "$SAM_FILE"


# SAM -> BAM
"$SAMTOOLS" view \
    -@ "$SAMTOOLS_THREADS" \
    -b \
    -o "$BAM_FILE" \
    "$SAM_FILE"


# Keep confident properly paired primary alignments
"$SAMTOOLS" view \
    -@ "$SAMTOOLS_THREADS" \
    -b \
    -f 2 \
    -F 2308 \
    -q 30 \
    -o "$STRICT_BAM" \
    "$BAM_FILE"


# Sort and index
"$SAMTOOLS" sort \
    -@ "$SAMTOOLS_THREADS" \
    -o "$SORTED_BAM" \
    "$STRICT_BAM"


"$SAMTOOLS" index \
    "$SORTED_BAM"


# Alignment statistics
"$SAMTOOLS" flagstat \
    "$SORTED_BAM" \
    > "$FLAGSTAT_REPORT"


# Number of mapped reads per chromosome/reference sequence
"$SAMTOOLS" idxstats \
    "$SORTED_BAM" \
    | awk '$1!="*"{print $1"\t"$3}' \
    | sort -k2,2nr \
    > "$CHROMOSOME_COUNTS"



# Done
echo
echo "Pipeline completed successfully."
echo
echo "Key outputs:"
echo "  LTR orientation report: $ORIENTATION_REPORT"
echo "  Deduplicated reads:      $DEDUP_R1 and $DEDUP_R2"
echo "  Strict sorted BAM:       $SORTED_BAM"
echo "  Alignment summary:       $FLAGSTAT_REPORT"
echo "  Reads per chromosome:    $CHROMOSOME_COUNTS"
echo
echo "The script itself does not store private input paths."
echo "Generated reports should still be checked before publishing in case"
echo "external programs included filenames or metadata."
