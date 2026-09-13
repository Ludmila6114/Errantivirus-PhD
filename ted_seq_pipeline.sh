#!/usr/bin/env bash

# TE/LTR paired-end processing pipeline
# ------------------------------------
# This script turns the workflow from the lab journal into a reusable pipeline.
# It intentionally contains NO private filesystem paths.
#
# Supply input paths at runtime as environment variables instead of editing them
# into this file. Example:
#
#   R1="/path/to/read1.fastq.gz" \
#   R2="/path/to/read2.fastq.gz" \
#   GENOME="/path/to/reference.fa" \
#   SAMPLE="sample_name" \
#   OUTDIR="./ted_seq_output" \
#   ./ted_seq_pipeline.sh
#
# Required software on PATH:
#   seqkit, fastp, bwa, samtools, gzip, awk, sort
#
# Optional environment variables:
#   SAMPLE            Output prefix. Default: sample
#   OUTDIR            Output directory. Default: ./ted_seq_output
#   MOTIF             LTR motif. Default: sequence recorded in the lab journal
#   BWA_THREADS       Threads for bwa mem. Default: 12
#   SAMTOOLS_THREADS  Threads for samtools. Default: 8
#   FASTP_THREADS     Threads for fastp. Default: 8
#   SEQKIT            seqkit executable. Default: seqkit
#   FASTP             fastp executable. Default: fastp
#   BWA               bwa executable. Default: bwa
#   SAMTOOLS          samtools executable. Default: samtools
#
# Important:
# - Do not commit private input paths to this script.
# - The workflow assumes R1/R2 identifiers have matching core read IDs.
# - The lab journal describes the LTR by a stated length and also provides a
#   literal sequence. This script uses the literal sequence exactly as recorded.

set -euo pipefail

# Never enable shell tracing here: `set -x` would print user-supplied paths.
set +x

###############################################################################
# 1. User-provided inputs and configurable parameters
###############################################################################

R1="${R1:-}"
R2="${R2:-}"
GENOME="${GENOME:-}"

SAMPLE="${SAMPLE:-sample}"
OUTDIR="${OUTDIR:-./ted_seq_output}"

# LTR motif recorded in the lab journal.
MOTIF="${MOTIF:-AGTTACCACCCCACCCCCTAAACCCCCACGCCTCTAAACAAAT}"

BWA_THREADS="${BWA_THREADS:-12}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-8}"
FASTP_THREADS="${FASTP_THREADS:-8}"

# Executable names can also be overridden at runtime if needed.
SEQKIT="${SEQKIT:-seqkit}"
FASTP="${FASTP:-fastp}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"

if [[ -z "$R1" || -z "$R2" || -z "$GENOME" ]]; then
    cat >&2 <<'USAGE'
ERROR: R1, R2, and GENOME must be provided as environment variables.

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

for input_file in "$R1" "$R2" "$GENOME"; do
    if [[ ! -f "$input_file" ]]; then
        echo "ERROR: Required input file was not found." >&2
        echo "Check the path supplied through R1, R2, or GENOME." >&2
        exit 1
    fi
done

# Check required programs without printing private paths.
for cmd in "$SEQKIT" "$FASTP" "$BWA" "$SAMTOOLS" gzip awk sort; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required command is not available: $cmd" >&2
        exit 1
    fi
done

###############################################################################
# 2. Output structure
###############################################################################

FILTER_DIR="$OUTDIR/01_ltr_filter"
TRIM_DIR="$OUTDIR/02_trimmed_flanks"
QC_DIR="$OUTDIR/03_fastp"
ALIGN_DIR="$OUTDIR/04_alignment"
STATS_DIR="$OUTDIR/05_stats"

mkdir -p "$FILTER_DIR" "$TRIM_DIR" "$QC_DIR" "$ALIGN_DIR" "$STATS_DIR"

# Output file names are generic and do not contain the private input paths.
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

###############################################################################
# 3. Calculate the reverse complement of the LTR motif
###############################################################################

RC="$({
    printf '>motif\n%s\n' "$MOTIF"
} | "$SEQKIT" seq -r -p 2>/dev/null | awk '!/^>/{printf "%s",$0} END{print ""}')"

if [[ -z "$RC" ]]; then
    echo "ERROR: Could not calculate the reverse complement of MOTIF." >&2
    exit 1
fi

echo "Starting TE/LTR pipeline for sample: $SAMPLE"
echo "Output directory: $OUTDIR"

###############################################################################
# Step 1: Check LTR orientation in the input reads
#
# The original workflow checked the forward motif and its reverse complement.
# The journal concluded that R1 contains the LTR in reverse-complement
# orientation, without mismatches. We report counts for both R1 and R2 as a
# sanity check, but filtering in Step 2 follows the journal and uses RC in R1.
###############################################################################

echo "[Step 1/6] Checking exact LTR motif orientation..."

count_exact_hits() {
    local fastq="$1"
    local pattern="$2"

    "$SEQKIT" grep -s -P -p "$pattern" "$fastq" 2>/dev/null \
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
    echo "Review the motif/orientation report before continuing with biological interpretation." >&2
fi

###############################################################################
# Step 2: Keep read pairs for which R1 contains the exact RC LTR motif
#
# 1. Select R1 reads containing the exact reverse-complement LTR sequence.
# 2. Extract their read IDs.
# 3. Retrieve the corresponding R2 reads using those IDs.
###############################################################################

echo "[Step 2/6] Selecting R1 reads containing the exact RC LTR motif..."

"$SEQKIT" grep -s -P -p "$RC" "$R1" -o "$KEPT_R1"
"$SEQKIT" seq -n -i "$KEPT_R1" > "$KEPT_IDS"
"$SEQKIT" grep -f "$KEPT_IDS" "$R2" -o "$KEPT_R2"

KEPT_COUNT="$($SEQKIT stats -T "$KEPT_R1" 2>/dev/null | awk 'NR==2{print $4}')"
if [[ -z "$KEPT_COUNT" || "$KEPT_COUNT" == "0" ]]; then
    echo "ERROR: No read pairs remained after exact LTR filtering." >&2
    exit 1
fi

echo "  Retained read pairs: $KEPT_COUNT"

###############################################################################
# Step 3: Trim reads so that only genomic flanking sequence is retained
#
# R2:
#   If the forward LTR motif is present, keep only sequence BEFORE the motif.
#   The motif itself and everything after it are removed.
#
# R1:
#   The RC motif is expected from Step 2. Remove everything through the end of
#   that RC motif and retain only sequence AFTER it.
#
# Sequence and quality strings are trimmed by the same number of bases.
###############################################################################

echo "[Step 3/6] Trimming LTR sequence and retaining genomic flanks..."

# Trim R2 at the first forward-motif occurrence.
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
    | gzip -c > "$TRIMMED_R2"

# Trim R1 at the first reverse-complement-motif occurrence.
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
    | gzip -c > "$TRIMMED_R1"

###############################################################################
# Step 4: Remove low-complexity / low-quality / too-short reads with fastp
#
# Parameters are preserved from the lab journal:
#   - adapter trimming disabled
#   - low-complexity filtering enabled, threshold 30
#   - qualified base threshold Q20
#   - at most 40% unqualified bases
#   - at most 5 N bases
#   - minimum length 10 nt
###############################################################################

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

###############################################################################
# Step 5: Remove PCR duplicates with fastp
#
# The second fastp pass performs deduplication only. Adapter trimming, quality
# filtering, and length filtering are disabled here because they were handled
# in Step 4.
###############################################################################

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

###############################################################################
# Step 6: Align to the reference genome with bwa mem and post-process with
# samtools
#
# Alignment threshold from the journal:
#   bwa mem -T 30
#
# Strict BAM filtering from the journal:
#   -f 2     keep properly paired alignments
#   -F 2308  remove unmapped, secondary, and supplementary alignments
#   -q 30    require MAPQ >= 30
###############################################################################

echo "[Step 6/6] Aligning retained genomic flanks to the reference genome..."

# Build the BWA index only if the standard BWA index files are not present.
# This writes index files alongside GENOME, so the reference location must be
# writable the first time the index is created.
if [[ ! -f "${GENOME}.bwt" || ! -f "${GENOME}.sa" ]]; then
    echo "  BWA index not detected; creating it now..."
    "$BWA" index "$GENOME"
fi

# BWA itself may print its command line to stderr. The script does not enable
# shell tracing and does not save that stderr output into repository files.
"$BWA" mem \
    -t "$BWA_THREADS" \
    -T 30 \
    "$GENOME" \
    "$DEDUP_R1" \
    "$DEDUP_R2" \
    > "$SAM_FILE"

# Convert SAM to BAM.
"$SAMTOOLS" view \
    -@ "$SAMTOOLS_THREADS" \
    -b \
    -o "$BAM_FILE" \
    "$SAM_FILE"

# Keep only confident, properly paired primary alignments.
"$SAMTOOLS" view \
    -@ "$SAMTOOLS_THREADS" \
    -b \
    -f 2 \
    -F 2308 \
    -q 30 \
    -o "$STRICT_BAM" \
    "$BAM_FILE"

# Coordinate-sort and index the strict BAM.
"$SAMTOOLS" sort \
    -@ "$SAMTOOLS_THREADS" \
    -o "$SORTED_BAM" \
    "$STRICT_BAM"

"$SAMTOOLS" index "$SORTED_BAM"

# Mapping summary.
"$SAMTOOLS" flagstat "$SORTED_BAM" > "$FLAGSTAT_REPORT"

# Count mapped reads per reference sequence/chromosome.
# Column 3 of `samtools idxstats` is the number of mapped reads.
"$SAMTOOLS" idxstats "$SORTED_BAM" \
    | awk '$1!="*"{print $1"\t"$3}' \
    | sort -k2,2nr \
    > "$CHROMOSOME_COUNTS"

###############################################################################
# Finished
###############################################################################

echo
echo "Pipeline completed successfully."
echo "Key outputs:"
echo "  LTR orientation report: $ORIENTATION_REPORT"
echo "  Deduplicated reads:      $DEDUP_R1 and $DEDUP_R2"
echo "  Strict sorted BAM:       $SORTED_BAM"
echo "  Alignment summary:       $FLAGSTAT_REPORT"
echo "  Reads per chromosome:    $CHROMOSOME_COUNTS"
echo
echo "Privacy note: the script itself contains no private input paths."
echo "Before publishing generated reports, review them for filenames or metadata"
echo "that external tools may have embedded."
