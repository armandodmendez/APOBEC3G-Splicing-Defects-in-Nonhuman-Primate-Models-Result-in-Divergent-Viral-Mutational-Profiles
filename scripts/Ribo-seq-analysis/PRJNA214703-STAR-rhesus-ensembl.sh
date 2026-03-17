#!/bin/bash

#$ -N PRJNA214703-STAR-rhesus-ensembl
#$ -cwd
#$ -pe smp 192
#$ -q all.q
#$ -wd /master/amendez/qsubs


# Author: Armando Mendez <amendez@txbiomed.org>
# Created on: 2025-07-08
# Modified on: 
#code to submit job: qsub -V -b n PRJNA214703-STAR.sh
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

# Define output directory
OUTPUT_DIR="/master/amendez/nhpdata/rna-seq/macaca-mulatta/PRJNA214703"
# Directory where the trimmed data are
DATA_DIR="/master/amendez/nhpdata/rna-seq/macaca-mulatta/PRJNA214703/trimmed-reads"
cd "$DATA_DIR"

##################### Run STAR aligner

# Directory containing your STAR index
GENOME="/master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl"

# Output directory for STAR
STAR_OUT_DIR="$OUTPUT_DIR/star"
mkdir -p "$STAR_OUT_DIR"

# Iterate over all trimmed FASTQ files (assumes single-end data)
for FASTQ_FILE in "$DATA_DIR"/*_trimmed.fastq.gz; do
    BASE_NAME=$(basename "$FASTQ_FILE" _trimmed.fastq.gz)
    OUTPUT_PREFIX="$STAR_OUT_DIR/${BASE_NAME}_mmul_10_"

    /master/amendez/app/STAR/STAR --runThreadN 12 \
        --genomeDir "$GENOME" \
        --readFilesIn "$FASTQ_FILE" \
        --outFileNamePrefix "$OUTPUT_PREFIX" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate
done

################### convert, merge and index files

BAM_FILES_PATH="$STAR_OUT_DIR"

# Index individual BAM files if not already indexed
for bam_file in "$BAM_FILES_PATH"/*mmul_10_Aligned.sortedByCoord.out.bam; do
    if [ ! -f "$bam_file.bai" ]; then
        samtools index "$bam_file"
    fi
done

# Merge and index BAM files if not already merged
cd "$STAR_OUT_DIR"
if [ ! -f "PRJNA214703-merged-rhemac10-2.bam" ]; then
    samtools merge -o "PRJNA214703-merged-mmul_10.bam" *mmul_10_Aligned.sortedByCoord.out.bam
    samtools index "PRJNA214703-merged-mmul_10.bam"
fi
