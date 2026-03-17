#!/bin/bash

#$ -N hg38-Ribo-TISH-predict-PRJNA292112
#$ -cwd
#$ -pe smp 192
#$ -q all.q
#$ -wd /master/amendez/qsubs


# Author: Armando Mendez <amendez@txbiomed.org>
# Created on: 2025-07-10
# Modified on: 
#==========#
# Comments #
#==========#
# rerunning with ensembl version of hg38
#==========#
# Versions #
#==========#

#=================#
# Initialize Conda #
#=================#
source /master/amendez/miniconda3/etc/profile.d/conda.sh || {
    echo "ERROR: Failed to initialize Conda"
    exit 1
}
conda activate ribo-tish || {
    echo "ERROR: Failed to activate ribo-tish environment"
    exit 1
}


# move into the directory where the data will be outputted

cd /master/amendez/nhpdata/ribo-seq/homo-sapiens/PRJNA292112/trimmed-reads/star/ribo-tish


ribotish predict -b /master/amendez/nhpdata/ribo-seq/homo-sapiens/PRJNA292112/trimmed-reads/star/PRJNA292112-merged-hg38.bam \
-g /master/amendez/nhpdata/genomes/all-genomes/Homo-sapiens-hg38/ensembl/Homo_sapiens.GRCh38.114.gtf \
-f /master/amendez/nhpdata/genomes/all-genomes/Homo-sapiens-hg38/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--framebest \
-o PRJNA292112-merged-hg38-ribo-tish-framebest.txt \
--seq \
--aaseq \
--blocks \
--alt \
-p 192 \

ribotish predict -b /master/amendez/nhpdata/ribo-seq/homo-sapiens/PRJNA292112/trimmed-reads/star/PRJNA292112-merged-hg38.bam \
-g /master/amendez/nhpdata/genomes/all-genomes/Homo-sapiens-hg38/ensembl/Homo_sapiens.GRCh38.114.gtf \
-f /master/amendez/nhpdata/genomes/all-genomes/Homo-sapiens-hg38/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--longest \
-o PRJNA292112-merged-hg38-ribo-tish-longest.txt \
--seq \
--aaseq \
--blocks \
--alt \
-p 192 \








