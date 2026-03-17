#!/bin/bash

#$ -N mmul10-Ribo-TISH-predict-PRJNA292112
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
# removed all problematic genes/transcripts with "Wrong CDS annotation", potentially due to overlapping genes
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

cd /master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112/trimmed-reads/bbduk/star/ribo-tish

#framebest Ribo-TISH
ribotish predict -b /master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112/trimmed-reads/bbduk/star/PRJNA292112-merged-mmul_10.bam \
-g /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.114.filtered2.gtf \
-f /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
--framebest \
-o PRJNA292112-merged-mmul_10_framebest.txt \
--seq \
--aaseq \
--blocks \
--alt \
-p 192

#framebest Ribo-TIH
ribotish predict -b /master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112/trimmed-reads/bbduk/star/PRJNA292112-merged-mmul_10.bam \
-g /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.114.filtered2.gtf \
-f /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
--framebest \
-o PRJNA292112-merged-mmul_10_framebest.txt \
--seq \
--aaseq \
--blocks \
--alt \
--genefilter ENSMMUG00000039346


#framebest Ribo-TIH
ribotish predict -b /master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112/trimmed-reads/bbduk/star/PRJNA292112-merged-mmul_10.bam \
-g /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.114.filtered2.gtf \
-f /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
--framebest \
-o PRJNA292112-merged-mmul_10_framebest.txt \
--seq \
--aaseq \
--blocks \
--genefilter ENSMMUG00000039346


ribotish predict -b /master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112/trimmed-reads/bbduk/star/PRJNA292112-merged-mmul_10.bam \
-g /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.114.filtered2.gtf \
-f /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
--longest \
-o PRJNA292112-merged-mmul_10_longest.txt \
--seq \
--aaseq \
--blocks \
--genefilter ENSMMUG00000039346

