#!/bin/bash

#$ -N PRJNA292112-STAR-hg38
#$ -cwd
#$ -pe smp 192
#$ -q all.q
#$ -wd /master/amendez/qsubs


# Author: Armando Mendez <amendez@txbiomed.org>
# Created on: 2025-07-07
# Modified on: 
#==========#
# Comments #
#==========#
- rerunning with ensembl version of hg38
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
DATA_DIR="/master/amendez/nhpdata/ribo-seq/homo-sapiens/PRJNA292112"
cd "$DATA_DIR"

# Directory to store the output files
OUTPUT_DIR="$DATA_DIR/trimmed-reads"
mkdir -p "$OUTPUT_DIR"

# Trimmomatic adapter file
ADAPTER_FILE="/master/amendez/nhpdata/truseq-adapters/TruSeq-adapter.fasta"

#################### Run Trimmomatic
for infile in *.fastq.gz; do
    base=$(basename "$infile" .fastq.gz)
    trimmomatic_out="$OUTPUT_DIR/${base}_trimmed.fastq.gz"

    if [ -f "$trimmomatic_out" ]; then
        echo "Skipping $infile: $trimmomatic_out already exists."
        continue
    fi

    trimmomatic SE -phred33 -threads 4 "$infile" \
        "$trimmomatic_out" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
        HEADCROP:1 SLIDINGWINDOW:4:20 MINLEN:15
done

#################### Run BBDuk to remove rRNA
# Directory for BBDuk output
BBDUK_OUT_DIR="$OUTPUT_DIR/bbduk"
mkdir -p "$BBDUK_OUT_DIR"

for infile in "$OUTPUT_DIR"/*_trimmed.fastq.gz; do
    base=$(basename "$infile" _trimmed.fastq.gz)
    bbduk_out="$BBDUK_OUT_DIR/${base}_filtered.fastq.gz"
    bbduk_rrna="$BBDUK_OUT_DIR/${base}_rrna.fastq.gz"

    if [ -f "$bbduk_out" ]; then
        echo "Skipping $infile: $bbduk_out already exists."
        continue
    fi

    # Note: BBDuk can read and write gzipped files directly
    bbduk.sh in="$infile" out="$bbduk_out" outm="$bbduk_rrna" \
        ref=/master/amendez/nhpdata/ribo-seq/General_EUK_longread_v1.9.4.fasta k=31
done

##################### Run STAR aligner
GENOME="/master/amendez/nhpdata/genomes/all-genomes/Homo-sapiens-hg38/ensembl"
STAR_OUT_DIR="$OUTPUT_DIR/star"
mkdir -p "$STAR_OUT_DIR"

for FASTQ_FILE in "$BBDUK_OUT_DIR"/*_filtered.fastq.gz; do
    BASE_NAME=$(basename "$FASTQ_FILE" _filtered.fastq.gz)
    OUTPUT_PREFIX="$STAR_OUT_DIR/${BASE_NAME}_hg38_"

    /master/amendez/app/STAR/STAR --runThreadN 12 \
        --genomeDir "$GENOME" \
        --readFilesIn "$FASTQ_FILE" \
        --outFileNamePrefix "$OUTPUT_PREFIX" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
	--outSAMattributes All
done

################### convert, merge and index files
BAM_FILES_PATH="$STAR_OUT_DIR"

for bam_file in "$BAM_FILES_PATH"/*hg38_Aligned.sortedByCoord.out.bam; do
    if [ ! -f "$bam_file.bai" ]; then
        samtools index "$bam_file"
    fi
done

cd "$STAR_OUT_DIR"
if [ ! -f "PRJNA292112-merged-hg38.bam" ]; then
    samtools merge -o "PRJNA292112-merged-hg38.bam" *hg38_Aligned.sortedByCoord.out.bam
    samtools index "PRJNA292112-merged-hg38.bam"
fi
