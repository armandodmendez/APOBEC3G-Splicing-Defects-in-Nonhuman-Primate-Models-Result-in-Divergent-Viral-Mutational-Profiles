#!/bin/bash

#$ -N mmul10-Ribo-TISH-QC-PRJNA292112
#$ -cwd
#$ -pe smp 1
#$ -q all.q
#$ -wd /master/amendez/qsubs


# Author: Armando Mendez <amendez@txbiomed.org>
# Created on: 2025-07-08
# Modified on: 
#==========#
# Comments #
#==========#
- rerunning with ensembl version of hg38
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

ribotish quality -b /master/amendez/nhpdata/ribo-seq/macaca-mulatta/PRJNA292112/trimmed-reads/bbduk/star/PRJNA292112-merged-mmul_10.bam -g /master/amendez/nhpdata/genomes/all-genomes/Macaca-mulatta/ensembl/Macaca_mulatta.Mmul_10.114.gtf --nom0









