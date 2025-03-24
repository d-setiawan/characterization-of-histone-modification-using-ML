#!/bin/bash

# Prompt for necessary inputs
read -p "Is the data paired-end? (y/n): " PAIRED_END

SRA_LIST=$1

REF_INDEX="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/tair10-index"
REGIONS_BED="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/TAIR10-ref-genome/TAIR10_genes_only.bed"
GREENSCREEN_BED="/home/dgsetiawan/MachineLearning/ZhongLab/ProofOfConcept/Data/chip-seq-data/arabidopsis_greenscreen_20inputs.bed"

# Read the SRA accession numbers from the list
while IFS= read -r SRA_ACCESSION; do
    # Create a folder for the current SRA accession number
    mkdir -p ${SRA_ACCESSION}
    cd ${SRA_ACCESSION}

    echo "Processing ${SRA_ACCESSION}..."

    # Step 1: Download the file
    fasterq-dump ${SRA_ACCESSION} --progress

    # Step 2: Define FASTQ filenames
    if [ "$PAIRED_END" == "y" ]; then
        FASTQ1="${SRA_ACCESSION}_1.fastq"
        FASTQ2="${SRA_ACCESSION}_2.fastq"
    else
        FASTQ1="${SRA_ACCESSION}.fastq"
    fi

    # Step 3: QC & Trim
    if [ "$PAIRED_END" == "y" ]; then
        fastp -i ${FASTQ1} -I ${FASTQ2} -o ${SRA_ACCESSION}_1_trimmed.fastq -O ${SRA_ACCESSION}_2_trimmed.fastq -h fastp_report.html
    else
        fastp -i ${FASTQ1} -o ${SRA_ACCESSION}_trimmed.fastq -h fastp_report.html
    fi

    # Step 4: Align to reference genome
    if [ "$PAIRED_END" == "y" ]; then
    bowtie2 -x ${REF_INDEX} -1 ${SRA_ACCESSION}_1_trimmed.fastq -2 ${SRA_ACCESSION}_2_trimmed.fastq -S ${SRA_ACCESSION}.sam --local -p 4 --no-unal \
    --rg-id ${SRA_ACCESSION} --rg "SM:${SRA_ACCESSION}" --rg "LB:${SRA_ACCESSION}" --rg "PL:ILLUMINA" --rg "PU:${SRA_ACCESSION}"
    else
    bowtie2 -x ${REF_INDEX} -U ${SRA_ACCESSION}_trimmed.fastq -S ${SRA_ACCESSION}.sam --local -p 4 --no-unal \
    --rg-id ${SRA_ACCESSION} --rg "SM:${SRA_ACCESSION}" --rg "LB:${SRA_ACCESSION}" --rg "PL:ILLUMINA" --rg "PU:${SRA_ACCESSION}"
    fi

    # Step 5: Convert SAM to BAM and sort
    samtools view -S -b ${SRA_ACCESSION}.sam > ${SRA_ACCESSION}.bam
    samtools sort ${SRA_ACCESSION}.bam -o ${SRA_ACCESSION}_sorted.bam

    # Step 5.1: Mark and Remove Duplicates Using Picard
    picard MarkDuplicates \
        -I ${SRA_ACCESSION}_sorted.bam \
        -O ${SRA_ACCESSION}_deduplicated.bam \
        -M ${SRA_ACCESSION}_duplication_metrics.txt \
        -REMOVE_DUPLICATES true

    # Step 5.2: Index the Deduplicated BAM File
    samtools index ${SRA_ACCESSION}_deduplicated.bam

    # Step 6: Generate BigWig file using deepTools
    bamCoverage -b ${SRA_ACCESSION}_deduplicated.bam -o ${SRA_ACCESSION}.bw --normalizeUsing RPKM -bl ${GREENSCREEN_BED}

    # Step 7: Compute matrix using BigWig file
    computeMatrix scale-regions -S ${SRA_ACCESSION}.bw -R ${REGIONS_BED} -m 5 --binSize 5 -p max/2 -o ${SRA_ACCESSION}.mat.gz

    echo "Processing of ${SRA_ACCESSION} completed."

    # Move back to the parent directory
    cd ..
done < ${SRA_LIST}

echo "All processes completed."
