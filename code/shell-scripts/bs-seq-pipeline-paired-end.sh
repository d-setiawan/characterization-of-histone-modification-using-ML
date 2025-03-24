#!/bin/bash

set -x

# Initialize Variables
read -p "Enter the input FASTQ file (for paired-end, provide the first FASTQ file): " INPUT_FASTQ
read -p "Is the data paired-end? (y/n): " PAIRED_END

BASE_NAME=$(basename "$INPUT_FASTQ" .fastq | cut -d'_' -f1)
REF_GENOME="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/bs-seq-ref"
CHROM_SIZE="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/TAIR10.chrom.sizes"
GENE_REGIONS="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/TAIR10_genes_only.bed"

# Perform QC and Trimming
if [ "$PAIRED_END" == "y" ]; then
    read -p "Enter the second paired FASTQ file: " INPUT_FASTQ2
    fastp -i "$INPUT_FASTQ" -I "$INPUT_FASTQ2" -o "${BASE_NAME}_1-trimmed.fastq" -O "${BASE_NAME}_2-trimmed.fastq"
else
    fastp -i "$INPUT_FASTQ" -o "${BASE_NAME}-trimmed.fastq"
fi

# Align and Quantify
if [ "$PAIRED_END" == "y" ]; then
    bismark --genome $REF_GENOME -1 "${BASE_NAME}_1-trimmed.fastq" -2 "${BASE_NAME}_2-trimmed.fastq" -p 8
else
    bismark --genome $REF_GENOME "${BASE_NAME}-trimmed.fastq" -p 8
fi

# Deduplicate the Samples
if [ "$PAIRED_END" == "y" ]; then
    deduplicate_bismark --bam "${BASE_NAME}_1-trimmed_bismark_bt2_pe.bam"
    DEDUPLICATED_BAM="${BASE_NAME}_1-trimmed_bismark_bt2_pe.deduplicated.bam"
else
    deduplicate_bismark --bam "${BASE_NAME}-trimmed_bismark_bt2.bam"
    DEDUPLICATED_BAM="${BASE_NAME}-trimmed_bismark_bt2.deduplicated.bam"
fi

# Extract the Methylation Values
bismark_methylation_extractor --CX --bedgraph "$DEDUPLICATED_BAM" --parallel 3

# Generate Bedgraphs for Cytosine Contexts
bismark2bedGraph -o "${BASE_NAME}.cpg.bedgraph" CpG_OT* CpG_OB*
bismark2bedGraph --CX -o "${BASE_NAME}.chg.bedgraph" CHG_OT* CHG_OB*
bismark2bedGraph --CX -o "${BASE_NAME}.chh.bedgraph" CHH_OT* CHH_OB*

gunzip "${BASE_NAME}.cpg.bedgraph.gz"
gunzip "${BASE_NAME}.chg.bedgraph.gz"
gunzip "${BASE_NAME}.chh.bedgraph.gz"

# Remove and Sort Bedgraph Files
sed '1d' "${BASE_NAME}.cpg.bedgraph" | sort -k1,1 -k2,2n > "${BASE_NAME}-sorted.cpg.bedgraph"
sed '1d' "${BASE_NAME}.chg.bedgraph" | sort -k1,1 -k2,2n > "${BASE_NAME}-sorted.chg.bedgraph"
sed '1d' "${BASE_NAME}.chh.bedgraph" | sort -k1,1 -k2,2n > "${BASE_NAME}-sorted.chh.bedgraph"

# Convert Bedgraph to BigWig
bedGraphToBigWig "${BASE_NAME}-sorted.cpg.bedgraph" $CHROM_SIZE "${BASE_NAME}-sorted.cpg.bw"
bedGraphToBigWig "${BASE_NAME}-sorted.chg.bedgraph" $CHROM_SIZE "${BASE_NAME}-sorted.chg.bw"
bedGraphToBigWig "${BASE_NAME}-sorted.chh.bedgraph" $CHROM_SIZE "${BASE_NAME}-sorted.chh.bw"

# Compute Matrix for All Cytosine Contexts
computeMatrix scale-regions -S "${BASE_NAME}-sorted.cpg.bw" -R $GENE_REGIONS -m 1 --binSize 1 -p max/2 -o "${BASE_NAME}-sorted.cpg-matrix.gz"
computeMatrix scale-regions -S "${BASE_NAME}-sorted.chg.bw" -R $GENE_REGIONS -m 1 --binSize 1 -p max/2 -o "${BASE_NAME}-sorted.chg-matrix.gz"
computeMatrix scale-regions -S "${BASE_NAME}-sorted.chh.bw" -R $GENE_REGIONS -m 1 --binSize 1 -p max/2 -o "${BASE_NAME}-sorted.chh-matrix.gz"

echo "Pipeline completed for ${BASE_NAME}."
