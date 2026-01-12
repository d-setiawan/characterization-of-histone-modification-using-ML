#!/usr/bin/env bash
set -euo pipefail 

# How to use:
usage() {
  echo "USAGE:"
  echo "  bash chipseq_simple.sh --config config/chipseq.env --sra SRR123 --outdir results [--sample-id NAME]"
  echo "  bash chipseq_simple.sh --config config/chipseq.env --fastq1 R1.fq.gz [--fastq2 R2.fq.gz] --outdir results --sample-id NAME"
  echo ""
  echo "NOTES:"
  echo "  - Config provides paths like GENOME_INDEX, CHROMSIZES, REGIONS_BED, THREADS"
  exit 1
}

# 1) Parse through args
CONFIG=""
SRA=""
FASTQ1=""
FASTQ2=""
OUTDIR="results"
SAMPLE_ID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
    --sra) SRA="$2"; shift 2;;
    --fastq1) FASTQ1="$2"; shift 2;;
    --fastq2) FASTQ2="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --sample-id) SAMPLE_ID="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown argument: $1"; usage;;
  esac
done

[[ -n "$CONFIG" ]] || usage
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG"; exit 1; }

# 2) Source config
# shellcheck source=/dev/null
source "$CONFIG"

# Required config variables
: "${THREADS:?Missing THREADS in config}"
: "${GENOME_INDEX:?Missing GENOME_INDEX in config}"
: "${REGIONS_BED:?Missing REGIONS_BED in config}"

# Optional config with defaults
MAPQ="${MAPQ:-30}"
BOWTIE2_EXTRA="${BOWTIE2_EXTRA:-}"
FASTP_EXTRA="${FASTP_EXTRA:-}"
DEDUP="${DEDUP:-1}"
PICARD_CMD="${PICARD_CMD:-picard}"   # requires "picard" in PATH
BAMCOVERAGE_NORM="${BAMCOVERAGE_NORM:-RPKM}"
BAMCOVERAGE_BIN_SIZE="${BAMCOVERAGE_BIN_SIZE:-10}"
BLACKLIST_BED="${BLACKLIST_BED:-}"
MATRIX_BIN_SIZE="${MATRIX_BIN_SIZE:-50}"
MATRIX_UPSTREAM="${MATRIX_UPSTREAM:-2000}"
MATRIX_DOWNSTREAM="${MATRIX_DOWNSTREAM:-2000}"
MATRIX_EXTRA="${MATRIX_EXTRA:-}"

# 3) Check input mode
if [[ -n "$SRA" ]]; then
  [[ -z "$FASTQ1" ]] || { echo "Use either --sra OR --fastq1, not both."; exit 1; }
  [[ -z "$FASTQ2" ]] || { echo "Use either --sra OR --fastq2, not both."; exit 1; }
else
  [[ -n "$FASTQ1" ]] || { echo "Provide --sra OR --fastq1"; exit 1; }
fi

# Decide sample id
if [[ -z "$SAMPLE_ID" ]]; then
  if [[ -n "$SRA" ]]; then
    SAMPLE_ID="$SRA"
  else
    # derive from fastq filename
    base="$(basename "$FASTQ1")"
    base="${base%.fastq.gz}"
    base="${base%.fq.gz}"
    base="${base%.fastq}"
    base="${base%.fq}"
    base="${base%_R1}"
    base="${base%_1}"
    SAMPLE_ID="$base"
  fi
fi

echo "[INFO] Sample: $SAMPLE_ID"

# 4) Create output directories for easier cleanup
SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
QC_DIR="${SAMPLE_DIR}/qc"
ALIGN_DIR="${SAMPLE_DIR}/align"
TRACK_DIR="${SAMPLE_DIR}/tracks"
MATRIX_DIR="${SAMPLE_DIR}/matrix"
TMP_DIR="${SAMPLE_DIR}/tmp"

mkdir -p "$QC_DIR" "$ALIGN_DIR" "$TRACK_DIR" "$MATRIX_DIR" "$TMP_DIR"

# Step 5.1: Download SRA if needed
if [[ -n "$SRA" ]]; then
  echo "[INFO] Downloading SRA: $SRA"
  fasterq-dump "$SRA" --threads "$THREADS" --outdir "$TMP_DIR"

  # detect paired-end vs single-end output
  if [[ -f "${TMP_DIR}/${SRA}_1.fastq" && -f "${TMP_DIR}/${SRA}_2.fastq" ]]; then
    gzip -f "${TMP_DIR}/${SRA}_1.fastq" "${TMP_DIR}/${SRA}_2.fastq"
    FASTQ1="${TMP_DIR}/${SRA}_1.fastq.gz"
    FASTQ2="${TMP_DIR}/${SRA}_2.fastq.gz"
  elif [[ -f "${TMP_DIR}/${SRA}.fastq" ]]; then
    gzip -f "${TMP_DIR}/${SRA}.fastq"
    FASTQ1="${TMP_DIR}/${SRA}.fastq.gz"
    FASTQ2=""
  else
    echo "[ERROR] Could not find FASTQ output from fasterq-dump for $SRA"
    exit 1
  fi
fi

# Step 5.2: Check FASTQ file exists
[[ -f "$FASTQ1" ]] || { echo "[ERROR] FASTQ1 not found: $FASTQ1"; exit 1; }
if [[ -n "$FASTQ2" ]]; then
  [[ -f "$FASTQ2" ]] || { echo "[ERROR] FASTQ2 not found: $FASTQ2"; exit 1; }
fi

# Step 6: Perform quality control and cleanup using fastp
echo "[INFO] fastp..."
CLEAN1="${TMP_DIR}/${SAMPLE_ID}.clean_R1.fastq.gz"
CLEAN2="${TMP_DIR}/${SAMPLE_ID}.clean_R2.fastq.gz"
FASTP_JSON="${QC_DIR}/fastp.${SAMPLE_ID}.json"
FASTP_HTML="${QC_DIR}/fastp.${SAMPLE_ID}.html"

if [[ -n "$FASTQ2" ]]; then
  fastp -w "$THREADS" -i "$FASTQ1" -I "$FASTQ2" -o "$CLEAN1" -O "$CLEAN2" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" $FASTP_EXTRA
else
  fastp -w "$THREADS" -i "$FASTQ1" -o "$CLEAN1" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" $FASTP_EXTRA
  CLEAN2=""
fi

# Step 7: Align and sort
echo "[INFO] bowtie2 alignment..."
SAM="${TMP_DIR}/${SAMPLE_ID}.sam"
BAM="${TMP_DIR}/${SAMPLE_ID}.bam"
SORTED_BAM="${ALIGN_DIR}/${SAMPLE_ID}.sorted.bam"

if [[ -n "$CLEAN2" ]]; then
  bowtie2 -x "$GENOME_INDEX" -1 "$CLEAN1" -2 "$CLEAN2" -p "$THREADS" $BOWTIE2_EXTRA > "$SAM"
else
  bowtie2 -x "$GENOME_INDEX" -U "$CLEAN1" -p "$THREADS" $BOWTIE2_EXTRA > "$SAM"
fi

samtools view -b -q "$MAPQ" "$SAM" > "$BAM"
samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM"
samtools index "$SORTED_BAM"

# Step 7.1: Mark and Remove Duplicates Using Picard (only if DEDUP=1)
FINAL_BAM="$SORTED_BAM"
if [[ "$DEDUP" == "1" ]]; then
  echo "[INFO] MarkDuplicates..."
  DEDUP_BAM="${ALIGN_DIR}/${SAMPLE_ID}.sorted.dedup.bam"
  METRICS="${ALIGN_DIR}/${SAMPLE_ID}.dedup.metrics.txt"
  picard MarkDuplicates I="$SORTED_BAM" O="$DEDUP_BAM" M="$METRICS" REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
  samtools index "$DEDUP_BAM"
  FINAL_BAM="$DEDUP_BAM"
else
  echo "[INFO] Skipping MarkDuplicates (DEDUP=0)"
fi

# Step 8: Generate BigWig file
echo "[INFO] bamCoverage..."
BW="${TRACK_DIR}/${SAMPLE_ID}.${BAMCOVERAGE_NORM}.bw"

BL_ARG=""
if [[ -n "$BLACKLIST_BED" ]]; then
  BL_ARG="--blackListFileName $BLACKLIST_BED"
fi

NORM_ARG=""
if [[ "$BAMCOVERAGE_NORM" != "None" ]]; then
  NORM_ARG="--normalizeUsing $BAMCOVERAGE_NORM"
fi

bamCoverage -b "$FINAL_BAM" -o "$BW" -p "$THREADS" --binSize "$BAMCOVERAGE_BIN_SIZE" $NORM_ARG $BL_ARG

# Step 9: Compute matrix of peaks using bigwig file
echo "[INFO] computeMatrix..."
MAT="${MATRIX_DIR}/${SAMPLE_ID}.matrix.gz"
TAB="${MATRIX_DIR}/${SAMPLE_ID}.matrix.tab"

computeMatrix reference-point \
  -S "$BW" \
  -R "$REGIONS_BED" \
  --referencePoint TSS \
  -b "$MATRIX_UPSTREAM" -a "$MATRIX_DOWNSTREAM" \
  --binSize "$MATRIX_BIN_SIZE" \
  -p "$THREADS" \
  -o "$MAT" \
  --outFileNameMatrix "$TAB" \
  $MATRIX_EXTRA

echo "[INFO] DONE"
echo "[INFO] Outputs in: $SAMPLE_DIR"
echo "  BAM:    $FINAL_BAM"
echo "  bigWig: $BW"
echo "  Matrix: $MAT"

# Optional: cleanup temp files 
# rm -rf "$TMP_DIR"
