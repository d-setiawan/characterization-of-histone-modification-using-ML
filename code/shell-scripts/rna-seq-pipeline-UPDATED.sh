#!/usr/bin/env bash
set -euo pipefail

# How to use:
usage() {
  echo "USAGE:"
  echo "  # Single-end SRA:"
  echo "  bash pipelines/rnaseq_rsem_simple.sh --config config/rnaseq.env --sra SRR123 --outdir results --matrix-name sampleA"
  echo ""
  echo "  # Paired-end SRA:"
  echo "  bash pipelines/rnaseq_rsem_simple.sh --config config/rnaseq.env --sra SRR123 --paired-end --outdir results --matrix-name sampleA"
  echo ""
  echo "NOTES:"
  echo "  - This matches your original flow:"
  echo "    fasterq-dump -> fastp -> rsem-prepare-reference -> rsem-calculate-expression"
  echo "  - Config provides: REF_FASTA, GFF3, RSEM_REF_PREFIX, THREADS"
  echo "  - If BUILD_RSEM_REF=0 in config, it will NOT rebuild the reference each run (recommended)."
  exit 1
}

# 1) Parse through args
CONFIG=""
SRA=""
OUTDIR="results"
MATRIX_NAME=""
PAIRED_END="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
    --sra) SRA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --matrix-name) MATRIX_NAME="$2"; shift 2;;
    --paired-end) PAIRED_END="1"; shift;;
    -h|--help) usage;;
    *) echo "Unknown argument: $1"; usage;;
  esac
done

[[ -n "$CONFIG" ]] || usage
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG"; exit 1; }
[[ -n "$SRA" ]] || { echo "Missing --sra"; usage; }
[[ -n "$MATRIX_NAME" ]] || { echo "Missing --matrix-name"; usage; }

# 2) Source config
# shellcheck source=/dev/null
source "$CONFIG"

# Required config variables
: "${THREADS:?Missing THREADS in config}"
: "${REF_FASTA:?Missing REF_FASTA in config}"
: "${GFF3:?Missing GFF3 in config}"
: "${RSEM_REF_PREFIX:?Missing RSEM_REF_PREFIX in config}"

# Optional config with defaults
FASTP_EXTRA="${FASTP_EXTRA:-}"
BUILD_RSEM_REF="${BUILD_RSEM_REF:-0}"       # recommended default: 0 (build once, reuse rsem reference)
RSEM_PREP_EXTRA="${RSEM_PREP_EXTRA:-}"
RSEM_CALC_EXTRA="${RSEM_CALC_EXTRA:-}"

# 3) Create output directories for easier cleanup
SAMPLE_ID="$SRA"
SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
QC_DIR="${SAMPLE_DIR}/qc"
FASTQ_DIR="${SAMPLE_DIR}/fastq"
RSEM_DIR="${SAMPLE_DIR}/rsem"
TMP_DIR="${SAMPLE_DIR}/tmp"

mkdir -p "$QC_DIR" "$FASTQ_DIR" "$RSEM_DIR" "$TMP_DIR"

echo "[INFO] Sample: $SAMPLE_ID"
echo "[INFO] Output dir: $SAMPLE_DIR"

# 4) Download SRA (fasterq-dump)
echo "[INFO] Downloading and converting SRA to FASTQ..."
fasterq-dump "$SRA" --progress --outdir "$FASTQ_DIR"

# Step 4.1: Determine input FASTQs based on paired-end flag
if [[ "$PAIRED_END" == "1" ]]; then
  FASTQ1="${FASTQ_DIR}/${SRA}_1.fastq"
  FASTQ2="${FASTQ_DIR}/${SRA}_2.fastq"
  [[ -f "$FASTQ1" ]] || { echo "[ERROR] Missing FASTQ1: $FASTQ1"; exit 1; }
  [[ -f "$FASTQ2" ]] || { echo "[ERROR] Missing FASTQ2: $FASTQ2"; exit 1; }
else
  FASTQ1="${FASTQ_DIR}/${SRA}.fastq"
  [[ -f "$FASTQ1" ]] || { echo "[ERROR] Missing FASTQ: $FASTQ1"; exit 1; }
fi

# 5) Perform quality control and trimming using fastp
# IMPORTANT: RSEM should use the TRIMMED outputs, not the raw FASTQs.
echo "[INFO] Trimming and quality control (fastp)..."
FASTP_JSON="${QC_DIR}/fastp.${SAMPLE_ID}.json"
FASTP_HTML="${QC_DIR}/fastp.${SAMPLE_ID}.html"

if [[ "$PAIRED_END" == "1" ]]; then
  TRIM1="${FASTQ_DIR}/${SRA}_1_trimmed.fastq"
  TRIM2="${FASTQ_DIR}/${SRA}_2_trimmed.fastq"
  fastp -w "$THREADS" \
    -i "$FASTQ1" -I "$FASTQ2" \
    -o "$TRIM1" -O "$TRIM2" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" \
    $FASTP_EXTRA
else
  TRIM1="${FASTQ_DIR}/${SRA}_trimmed.fastq"
  fastp -w "$THREADS" \
    -i "$FASTQ1" \
    -o "$TRIM1" \
    --json "$FASTP_JSON" --html "$FASTP_HTML" \
    $FASTP_EXTRA
fi

# 6) Create RSEM reference (optional per run)
# Set BUILD_RSEM_REF=1 if you want to rebuild each run.
if [[ "$BUILD_RSEM_REF" == "1" ]]; then
  echo "[INFO] Creating RSEM reference..."
  rsem-prepare-reference -p "$THREADS" --bowtie2 --gff3 "$GFF3" "$REF_FASTA" "$RSEM_REF_PREFIX" $RSEM_PREP_EXTRA
else
  echo "[INFO] Skipping RSEM reference build (BUILD_RSEM_REF=0). Assuming reference already exists at prefix: $RSEM_REF_PREFIX"
fi

# 7) Calculate expression (RSEM) to create output matrix/results
# Outputs will be written with prefix: <RSEM_DIR>/<MATRIX_NAME>.*
echo "[INFO] Calculating expression (RSEM)..."
OUT_PREFIX="${RSEM_DIR}/${MATRIX_NAME}"

if [[ "$PAIRED_END" == "1" ]]; then
  rsem-calculate-expression \
    --bowtie2 --estimate-rspd -p "$THREADS" --no-bam-output --paired-end \
    "$TRIM1" "$TRIM2" \
    "$RSEM_REF_PREFIX" \
    "$OUT_PREFIX" \
    $RSEM_CALC_EXTRA
else
  rsem-calculate-expression \
    --bowtie2 --estimate-rspd -p "$THREADS" --no-bam-output \
    "$TRIM1" \
    "$RSEM_REF_PREFIX" \
    "$OUT_PREFIX" \
    $RSEM_CALC_EXTRA
fi

echo "[INFO] DONE"
echo "[INFO] Outputs in: $SAMPLE_DIR"
echo "  QC:    $QC_DIR"
echo "  FASTQ: $FASTQ_DIR"
echo "  RSEM:  $RSEM_DIR"
echo "  Prefix: $OUT_PREFIX"
