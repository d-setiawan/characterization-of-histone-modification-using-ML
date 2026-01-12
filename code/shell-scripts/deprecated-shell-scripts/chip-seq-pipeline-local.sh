#!/bin/bash

# Prompt user for FASTQ files and other inputs
FASTQ1=$1
read -p "Is the data paired-end? (y/n): " PAIRED_END
if [ "$PAIRED_END" == "y" ]; then
    read -p "Enter the path to the second FASTQ file: " FASTQ2
fi
read -p "Enter name of the output matrix: " MATRIX_NAME

REF_INDEX="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/tair10-index"
REGIONS_BED="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/TAIR10_genes_only.bed"
GREENSCREEN_BED="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/chip-seq-data/arabidopsis_greenscreen_20inputs.bed"

# Extract base name for output files
BASE_NAME=$(basename "$FASTQ1" | sed 's/_1.fastq$//; s/.fastq$//')

# Step 1: QC & Trim
if [ "$PAIRED_END" == "y" ]; then
    fastp -i ${FASTQ1} -I ${FASTQ2} -o ${BASE_NAME}_1_trimmed.fastq -O ${BASE_NAME}_2_trimmed.fastq -h fastp_report.html
else
    fastp -i ${FASTQ1} -o ${BASE_NAME}_trimmed.fastq -h fastp_report.html
fi

# Step 2: Align to reference genome

if [ "$PAIRED_END" == "y" ]; then
bowtie2 -x ${REF_INDEX} -1 ${BASE_NAME}_1_trimmed.fastq -2 ${BASE_NAME}_2_trimmed.fastq -S ${BASE_NAME}.sam --local -p 4 --no-unal \
 --rg-id ${BASE_NAME} --rg "SM:${BASE_NAME}" --rg "LB:${BASE_NAME}" --rg "PL:ILLUMINA" --rg "PU:${BASE_NAME}"
else
bowtie2 -x ${REF_INDEX} -U ${BASE_NAME}_trimmed.fastq -S ${BASE_NAME}.sam --local -p 4 --no-unal \
 --rg-id ${BASE_NAME} --rg "SM:${BASE_NAME}" --rg "LB:${BASE_NAME}" --rg "PL:ILLUMINA" --rg "PU:${BASE_NAME}"
fi

# Step 3: Convert SAM to BAM and sort
samtools view -S -b ${BASE_NAME}.sam > ${BASE_NAME}.bam
samtools sort ${BASE_NAME}.bam -o ${BASE_NAME}_sorted.bam

# Step 4: Mark and remove duplicates using Picard (via Conda)
picard MarkDuplicates \
    -I ${BASE_NAME}_sorted.bam \
    -O ${BASE_NAME}_deduplicated.bam \
    -M ${BASE_NAME}_duplication_metrics.txt \
    -REMOVE_DUPLICATES true

# Step 5: Index the deduplicated BAM file
samtools index ${BASE_NAME}_deduplicated.bam

# Step 6: Generate BigWig file using deepTools
bamCoverage -b ${BASE_NAME}_deduplicated.bam -o ${BASE_NAME}.bw --normalizeUsing RPKM -bl ${GREENSCREEN_BED}

# Step 7: Compute matrix using BigWig file
computeMatrix scale-regions -S ${BASE_NAME}.bw -R ${REGIONS_BED} -m 1000 --binSize 1000 -p max/2 -o ${MATRIX_NAME}.gz

echo "Process completed for ${BASE_NAME}."
