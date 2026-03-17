#!/bin/bash

#$ -N PRJNA122271-STAR-hg38
#$ -cwd
#$ -pe smp 192
#$ -q all.q
#$ -wd /master/amendez/qsubs


# Author: Armando Mendez <amendez@txbiomed.org>
# Created on: 2025-11-13
# Modified on: 
#==========#
# Comments #
#==========#

#==========#
# Versions #
#==========#

# v1.0 - 2024-06-05: script created

#=================#
# Initialize Conda #
#=================#
source /master/amendez/miniconda3/etc/profile.d/conda.sh || {
    echo "ERROR: Failed to initialize Conda"
    exit 1
}
conda activate rnaseq || {
    echo "ERROR: Failed to activate rnaseq environment"
    exit 1
}

# fastqc and multiqc have been done

# Directory where the data are
DATA_DIR="/master/amendez/nhpdata/rna-seq/homo-sapiens/PRJNA122271"
cd "$DATA_DIR"

# Directory to store the output files
OUTPUT_DIR="$DATA_DIR/trimmed-reads"
mkdir -p "$OUTPUT_DIR"

# Trimmomatic adapter file
ADAPTER_FILE="/master/amendez/nhpdata/truseq-adapters/TruSeq-adapter.fasta"

#################### Run trimmomatic
for infile in *.fastq.gz; do
    base=$(basename "$infile" .fastq.gz)
    output_file="$OUTPUT_DIR/${base}_trimmed.fastq.gz"

    # Skip if output file already exists
    if [ -f "$output_file" ]; then
        echo "Skipping $infile: $output_file already exists."
        continue
    fi

    # Run trimmomatic
    trimmomatic SE -phred33 -threads 4 "$infile" \
        "$output_file" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
        HEADCROP:13 TRAILING:28 SLIDINGWINDOW:4:20 MINLEN:20
done

##################### Run STAR aligner

# Directory containing your STAR index
GENOME="/master/amendez/nhpdata/genomes/all-genomes/Homo-sapiens-hg38/ensembl"

# Output directory for STAR
STAR_OUT_DIR="$OUTPUT_DIR/star"
mkdir -p "$STAR_OUT_DIR"

# Iterate over all trimmed FASTQ files (assumes single-end data)
for FASTQ_FILE in "$OUTPUT_DIR"/*_trimmed.fastq.gz; do
    BASE_NAME=$(basename "$FASTQ_FILE" _trimmed.fastq.gz)
    OUTPUT_PREFIX="$STAR_OUT_DIR/${BASE_NAME}_hg38_"

    /master/amendez/app/STAR/STAR --runThreadN 12 \
        --genomeDir "$GENOME" \
        --readFilesIn "$FASTQ_FILE" \
        --outFileNamePrefix "$OUTPUT_PREFIX" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMprimaryFlag AllBestScore \
        --outSAMmultNmax 1

done

################### convert, merge and index files

BAM_FILES_PATH="$STAR_OUT_DIR"

# Index individual BAM files if not already indexed
for bam_file in "$BAM_FILES_PATH"/*hg38_Aligned.sortedByCoord.out.bam; do
    if [ ! -f "$bam_file.bai" ]; then
        samtools index "$bam_file"
    fi
done

cd "$STAR_OUT_DIR"
if [ ! -f "PRJNA122271-merged-hg38.bam" ]; then
    samtools merge -o "PRJNA214703-merged-hg38.bam" *hg38_Aligned.sortedByCoord.out.bam
    samtools index "PRJNA214703-merged-hg38.bam"
fi
