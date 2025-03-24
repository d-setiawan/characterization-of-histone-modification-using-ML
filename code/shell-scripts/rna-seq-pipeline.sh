#!/bin/bash

# Step 1: Prompt for SRA accession number, paired end, and matrix name
read -p "Enter the SRA accession number: " SRA_ACCESSION
read -p "Is the data paired-end (y/n)? " PAIRED_END
read -p "Enter the name for the end matrix for RSEM calculate expression: " MATRIX_NAME

REF_GENOME_PATH="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/TAIR10_chr_all.fas"
GENE_ANNOTATION_PATH="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/TAIR10_GFF3_genes.gff"

# Step 2: Download data from SRA and convert to FASTQ
echo "Downloading and converting SRA to FASTQ..."
fasterq-dump $SRA_ACCESSION --progress

# Step 3: Define input and output file names based on paired-end status
if [ "$PAIRED_END" == "y" ]; then
    FASTQ1="${SRA_ACCESSION}_1.fastq"
    FASTQ2="${SRA_ACCESSION}_2.fastq"
else
    FASTQ1="${SRA_ACCESSION}.fastq"
fi

# Step 4: Trim and perform quality control
echo "Trimming and quality control..."
if [ "$PAIRED_END" == "y" ]; then
    fastp -i ${FASTQ1} -I ${FASTQ2} -o ${SRA_ACCESSION}_1_trimmed.fastq -O ${SRA_ACCESSION}_2_trimmed.fastq -h fastp_report.html
else
    fastp -i ${FASTQ1} -o ${SRA_ACCESSION}_trimmed.fastq -h fastp_report.html
fi

# Step 5: Create RSEM reference
echo "Creating RSEM reference..."
rsem-prepare-reference -p 10 --bowtie2 --gff3 ${GENE_ANNOTATION_PATH} ${REF_GENOME_PATH} arabidopsis-ref

# Step 6: Calculate expression to create matrix
echo "Calculating expression..."
if [ "$PAIRED_END" == "y" ]; then
    rsem-calculate-expression --bowtie2 --estimate-rspd -p 10 --no-bam-output --paired-end ${FASTQ1} ${FASTQ2} arabidopsis-ref $MATRIX_NAME
else
    rsem-calculate-expression --bowtie2 --estimate-rspd -p 10 --no-bam-output ${FASTQ1} arabidopsis-ref $MATRIX_NAME
fi

echo "Process completed."
