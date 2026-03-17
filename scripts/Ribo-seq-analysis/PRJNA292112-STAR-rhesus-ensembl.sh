#!/bin/bash

#$ -N PRJNA292112-STAR-rhesus-ensembl
#$ -cwd
#$ -pe smp 192
#$ -q all.q
#$ -wd /master/amendez/qsubs
#$ -l h_vmem=980G
#$ -l mem_free=980G

# Author: Armando Mendez
# Created on: 2025-07-08

#==========#
# Comments #
#==========#
# genome already indexed

#==========#
# Versions #
#==========#

#===========#
# Variables #
#===========#

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

# Directory containing raw FASTQ files
RAW_FASTQ_DIR="/master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112"

# Directory for trimmed FASTQ files
FASTQ_DIR="/master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112/trimmed-reads"
mkdir -p "${FASTQ_DIR}"

# Directory containing your STAR index
genome="/master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl"

# Output directory for BBDuk
BBDUK_DIR="${FASTQ_DIR}/bbduk"
mkdir -p "${BBDUK_DIR}"

# Output directory for STAR
OUT_DIR="${FASTQ_DIR}/bbduk/star"
mkdir -p "${OUT_DIR}"

# BBDuk rRNA reference
BBDUK_REF="/master/amendez/nhpdata/ribo-seq/General_EUK_longread_v1.9.4.fasta"

# Trimmomatic adapter file
ADAPTER_FILE="/master/amendez/nhpdata/truseq-adapters/TruSeq-adapter.fasta"

#========================#
# Trimmomatic Trimming   #
#========================#
for infile in ${RAW_FASTQ_DIR}/*.fastq.gz; do
    base=$(basename "$infile" .fastq.gz)
    trimmomatic_out="${FASTQ_DIR}/${base}_trimmed.fastq.gz"

    if [ -f "$trimmomatic_out" ]; then
        echo "Skipping $infile: $trimmomatic_out already exists."
        continue
    fi

    echo "Running Trimmomatic for $infile..."
    trimmomatic SE -phred33 -threads 4 "$infile" \
        "$trimmomatic_out" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
        HEADCROP:1 SLIDINGWINDOW:4:20 MINLEN:15
done

#================#
# BBDuk Filtering
#================#
for TRIMMED_FASTQ in ${FASTQ_DIR}/*_trimmed.fastq.gz; do
    BASE_NAME=$(basename "${TRIMMED_FASTQ}" _trimmed.fastq.gz)
    FILTERED_FASTQ="${BBDUK_DIR}/${BASE_NAME}_filtered.fastq.gz"
    RRNA_FASTQ="${BBDUK_DIR}/${BASE_NAME}_rrna.fastq.gz"

    # Skip if already filtered
    if [ -f "${FILTERED_FASTQ}" ]; then
        echo "Skipping BBDuk for ${TRIMMED_FASTQ}: output exists."
        continue
    fi

    echo "Running BBDuk for ${TRIMMED_FASTQ}..."
    bbduk.sh in="${TRIMMED_FASTQ}" out="${FILTERED_FASTQ}" outm="${RRNA_FASTQ}" ref="${BBDUK_REF}" k=31
done

#================#
# STAR Alignment
#================#
for FILTERED_FASTQ in ${BBDUK_DIR}/*_filtered.fastq.gz; do
    BASE_NAME=$(basename "${FILTERED_FASTQ}" _filtered.fastq.gz)
    STAR_PREFIX="${OUT_DIR}/${BASE_NAME}_mmul_10_"

    /master/amendez/app/STAR/STAR \
        --runThreadN 12 \
        --genomeDir "${genome}" \
        --readFilesIn "${FILTERED_FASTQ}" \
        --outFileNamePrefix "${STAR_PREFIX}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
	--outSAMattributes All
done

#=============================#
# BAM Indexing and Merging    #
#=============================#
bam_files_path="${OUT_DIR}"

for bam_file in ${bam_files_path}/*mmul_10_Aligned.sortedByCoord.out.bam; do
    samtools index "${bam_file}"
done

cd "${bam_files_path}"
samtools merge -o PRJNA292112-merged-mmul_10.bam *mmul_10_Aligned.sortedByCoord.out.bam
samtools index PRJNA292112-merged-mmul_10.bam
