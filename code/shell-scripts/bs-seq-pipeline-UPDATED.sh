#!/usr/bin/env bash
set -euo pipefail 

# How to use:
usage() {
  echo "USAGE:"
  echo "  # Local FASTQ (single-end):"
  echo "  bash pipelines/bsseq_simple_matrix.sh --config config/bsseq.env --fastq1 sample.fastq.gz --outdir results [--sample-id NAME]"
  echo ""
  echo "  # Local FASTQ (paired-end):"
  echo "  bash pipelines/bsseq_simple_matrix.sh --config config/bsseq.env --fastq1 sample_R1.fastq.gz --fastq2 sample_R2.fastq.gz --outdir results [--sample-id NAME]"
  echo ""
  echo "NOTES:"
  echo "  - This script matches your original flow:"
  echo "    fastp -> bismark -> deduplicate -> methylation_extractor -> bismark2bedGraph -> bedGraphToBigWig -> computeMatrix"
  echo "  - Config provides: BISMARK_GENOME_DIR, CHROMSIZES, GENE_REGIONS_BED, THREADS, and computeMatrix params"
  exit 1
}

# 1) Parse through args
CONFIG=""
FASTQ1=""
FASTQ2=""
OUTDIR="results"
SAMPLE_ID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
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
: "${BISMARK_GENOME_DIR:?Missing BISMARK_GENOME_DIR in config}"
: "${CHROMSIZES:?Missing CHROMSIZES in config}"
: "${GENE_REGIONS_BED:?Missing GENE_REGIONS_BED in config}"

# Optional config with defaults 
FASTP_EXTRA="${FASTP_EXTRA:-}"
BISMARK_EXTRA="${BISMARK_EXTRA:-}"
METH_EXTRACTOR_PARALLEL="${METH_EXTRACTOR_PARALLEL:-3}"

# computeMatrix params 
MATRIX_MODE="${MATRIX_MODE:-scale-regions}"
MATRIX_REGION_BODY_LENGTH="${MATRIX_REGION_BODY_LENGTH:-1}"
MATRIX_BIN_SIZE="${MATRIX_BIN_SIZE:-1}"
MATRIX_THREADS="${MATRIX_THREADS:-max/2}"   # deepTools accepts "max/2" string

# 3) Check input mode (this cleaned version is local FASTQ only)
[[ -n "$FASTQ1" ]] || { echo "Provide --fastq1"; exit 1; }

# Check FASTQ file exists
[[ -f "$FASTQ1" ]] || { echo "[ERROR] FASTQ1 not found: $FASTQ1"; exit 1; }
if [[ -n "$FASTQ2" ]]; then
  [[ -f "$FASTQ2" ]] || { echo "[ERROR] FASTQ2 not found: $FASTQ2"; exit 1; }
fi

# Decide sample id
if [[ -z "$SAMPLE_ID" ]]; then
  # derive from fastq filename
  base="$(basename "$FASTQ1")"
  base="${base%.fastq.gz}"
  base="${base%.fq.gz}"
  base="${base%.fastq}"
  base="${base%.fq}"
  base="${base%_R1}"
  base="${base%_1}"

  SAMPLE_ID="$(echo "$base" | cut -d'_' -f1)"
fi

echo "[INFO] Sample: $SAMPLE_ID"

# 4) Create output directories for easier cleanup
SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
QC_DIR="${SAMPLE_DIR}/qc"
ALIGN_DIR="${SAMPLE_DIR}/align"
METH_DIR="${SAMPLE_DIR}/methylation"
TRACK_DIR="${SAMPLE_DIR}/tracks"
MATRIX_DIR="${SAMPLE_DIR}/matrix"
TMP_DIR="${SAMPLE_DIR}/tmp"

mkdir -p "$QC_DIR" "$ALIGN_DIR" "$METH_DIR" "$TRACK_DIR" "$MATRIX_DIR" "$TMP_DIR"

# Step 5: Perform QC and trimming using fastp
# Output trimmed fastqs are placed in tmp/ and named consistently
echo "[INFO] fastp..."
TRIM1="${TMP_DIR}/${SAMPLE_ID}_1-trimmed.fastq.gz"
TRIM2="${TMP_DIR}/${SAMPLE_ID}_2-trimmed.fastq.gz"
FASTP_JSON="${QC_DIR}/fastp.${SAMPLE_ID}.json"
FASTP_HTML="${QC_DIR}/fastp.${SAMPLE_ID}.html"

if [[ -n "$FASTQ2" ]]; then
  fastp -w "$THREADS" -i "$FASTQ1" -I "$FASTQ2" -o "$TRIM1" -O "$TRIM2" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" $FASTP_EXTRA
else

  TRIM1="${TMP_DIR}/${SAMPLE_ID}-trimmed.fastq.gz"
  fastp -w "$THREADS" -i "$FASTQ1" -o "$TRIM1" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" $FASTP_EXTRA
  TRIM2=""
fi

# Step 6: Align using Bismark (writes BAM and reports into align/)
echo "[INFO] bismark alignment..."
if [[ -n "$TRIM2" ]]; then
  bismark --genome "$BISMARK_GENOME_DIR" \
    -1 "$TRIM1" -2 "$TRIM2" \
    -p "$THREADS" \
    -o "$ALIGN_DIR" \
    $BISMARK_EXTRA
else
  bismark --genome "$BISMARK_GENOME_DIR" \
    "$TRIM1" \
    -p "$THREADS" \
    -o "$ALIGN_DIR" \
    $BISMARK_EXTRA
fi

# Step 6.1: Identify the expected Bismark BAM name (deterministic based on our trimmed filenames)
if [[ -n "$TRIM2" ]]; then
  # paired-end: <trim1_base>_bismark_bt2_pe.bam
  TRIM1_BASE="$(basename "$TRIM1")"
  TRIM1_BASE="${TRIM1_BASE%.fastq.gz}"
  TRIM1_BASE="${TRIM1_BASE%.fastq}"
  TRIM1_BASE="${TRIM1_BASE%.fq.gz}"
  TRIM1_BASE="${TRIM1_BASE%.fq}"
  BISMARK_BAM="${ALIGN_DIR}/${TRIM1_BASE}_bismark_bt2_pe.bam"
else
  # single-end: <trim_base>_bismark_bt2.bam
  TRIM_BASE="$(basename "$TRIM1")"
  TRIM_BASE="${TRIM_BASE%.fastq.gz}"
  TRIM_BASE="${TRIM_BASE%.fastq}"
  TRIM_BASE="${TRIM_BASE%.fq.gz}"
  TRIM_BASE="${TRIM_BASE%.fq}"
  BISMARK_BAM="${ALIGN_DIR}/${TRIM_BASE}_bismark_bt2.bam"
fi

[[ -f "$BISMARK_BAM" ]] || { echo "[ERROR] Expected Bismark BAM not found: $BISMARK_BAM"; exit 1; }
echo "[INFO] Bismark BAM: $BISMARK_BAM"

# Step 7: Deduplicate the samples (produces *.deduplicated.bam next to input)
echo "[INFO] deduplicate_bismark..."
if [[ -n "$TRIM2" ]]; then
  deduplicate_bismark --bam --paired "$BISMARK_BAM"
  DEDUP_BAM="${BISMARK_BAM%.bam}.deduplicated.bam"
else
  deduplicate_bismark --bam "$BISMARK_BAM"
  DEDUP_BAM="${BISMARK_BAM%.bam}.deduplicated.bam"
fi

[[ -f "$DEDUP_BAM" ]] || { echo "[ERROR] Deduplicated BAM not found: $DEDUP_BAM"; exit 1; }
echo "[INFO] Deduplicated BAM: $DEDUP_BAM"

# Step 8: Extract methylation values 

echo "[INFO] bismark_methylation_extractor..."
bismark_methylation_extractor \
  --CX \
  --bedGraph \
  --parallel "$METH_EXTRACTOR_PARALLEL" \
  --gzip \
  -o "$METH_DIR" \
  "$DEDUP_BAM"

# Step 9: Generate bedGraphs for cytosine contexts (CpG, CHG, CHH) using bismark2bedGraph
# Run from within METH_DIR so the CpG_OT*/CHG_OT* globs resolve correctly.
echo "[INFO] bismark2bedGraph (CpG/CHG/CHH)..."
pushd "$METH_DIR" >/dev/null

bismark2bedGraph -o "${SAMPLE_ID}.cpg.bedgraph" CpG_OT* CpG_OB*
bismark2bedGraph --CX -o "${SAMPLE_ID}.chg.bedgraph" CHG_OT* CHG_OB*
bismark2bedGraph --CX -o "${SAMPLE_ID}.chh.bedgraph" CHH_OT* CHH_OB*

# Step 9.1: Unzip bedGraphs 
gunzip -f "${SAMPLE_ID}.cpg.bedgraph.gz"
gunzip -f "${SAMPLE_ID}.chg.bedgraph.gz"
gunzip -f "${SAMPLE_ID}.chh.bedgraph.gz"

# Step 10: Remove header line and sort bedGraphs (required for bedGraphToBigWig)
echo "[INFO] sort bedGraphs..."
sed '1d' "${SAMPLE_ID}.cpg.bedgraph" | sort -k1,1 -k2,2n > "${TRACK_DIR}/${SAMPLE_ID}-sorted.cpg.bedgraph"
sed '1d' "${SAMPLE_ID}.chg.bedgraph" | sort -k1,1 -k2,2n > "${TRACK_DIR}/${SAMPLE_ID}-sorted.chg.bedgraph"
sed '1d' "${SAMPLE_ID}.chh.bedgraph" | sort -k1,1 -k2,2n > "${TRACK_DIR}/${SAMPLE_ID}-sorted.chh.bedgraph"

popd >/dev/null

# Step 11: Convert bedGraph to BigWig
echo "[INFO] bedGraphToBigWig..."
CPG_BW="${TRACK_DIR}/${SAMPLE_ID}-sorted.cpg.bw"
CHG_BW="${TRACK_DIR}/${SAMPLE_ID}-sorted.chg.bw"
CHH_BW="${TRACK_DIR}/${SAMPLE_ID}-sorted.chh.bw"

bedGraphToBigWig "${TRACK_DIR}/${SAMPLE_ID}-sorted.cpg.bedgraph" "$CHROMSIZES" "$CPG_BW"
bedGraphToBigWig "${TRACK_DIR}/${SAMPLE_ID}-sorted.chg.bedgraph" "$CHROMSIZES" "$CHG_BW"
bedGraphToBigWig "${TRACK_DIR}/${SAMPLE_ID}-sorted.chh.bedgraph" "$CHROMSIZES" "$CHH_BW"

# Step 12: Compute Matrix for all cytosine contexts 
echo "[INFO] computeMatrix (CpG/CHG/CHH)..."
CPG_MAT="${MATRIX_DIR}/${SAMPLE_ID}-sorted.cpg-matrix.gz"
CHG_MAT="${MATRIX_DIR}/${SAMPLE_ID}-sorted.chg-matrix.gz"
CHH_MAT="${MATRIX_DIR}/${SAMPLE_ID}-sorted.chh-matrix.gz"

computeMatrix "$MATRIX_MODE" -S "$CPG_BW" -R "$GENE_REGIONS_BED" -m "$MATRIX_REGION_BODY_LENGTH" --binSize "$MATRIX_BIN_SIZE" -p "$MATRIX_THREADS" -o "$CPG_MAT"
computeMatrix "$MATRIX_MODE" -S "$CHG_BW" -R "$GENE_REGIONS_BED" -m "$MATRIX_REGION_BODY_LENGTH" --binSize "$MATRIX_BIN_SIZE" -p "$MATRIX_THREADS" -o "$CHG_MAT"
computeMatrix "$MATRIX_MODE" -S "$CHH_BW" -R "$GENE_REGIONS_BED" -m "$MATRIX_REGION_BODY_LENGTH" --binSize "$MATRIX_BIN_SIZE" -p "$MATRIX_THREADS" -o "$CHH_MAT"

echo "[INFO] Pipeline completed for ${SAMPLE_ID}."
echo "[INFO] Outputs in: $SAMPLE_DIR"
echo "  QC:      $QC_DIR"
echo "  Align:   $ALIGN_DIR"
echo "  Meth:    $METH_DIR"
echo "  Tracks:  $TRACK_DIR"
echo "  Matrix:  $MATRIX_DIR"
echo "    - $CPG_MAT"
echo "    - $CHG_MAT"
echo "    - $CHH_MAT"

# Optional: cleanup temp files
# rm -rf "$TMP_DIR"
