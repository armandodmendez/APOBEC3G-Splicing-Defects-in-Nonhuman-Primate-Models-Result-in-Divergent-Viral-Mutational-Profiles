#!/bin/bash
#$ -N pacbio-lima-HIV-SIV.sh
#$ -cwd
#$ -pe smp 192
#$ -q all.q
#$ -wd /master/amendez/qsubs

# Version: 1.2
# Author: Armando Mendez <amendez@txbiomed.org>
# Updated on: 2025-11-21
#==========#
# Changes: #
#==========#

#=================#
# Initialize Conda #
#=================#
source /master/amendez/miniconda3/etc/profile.d/conda.sh || {
    echo "ERROR: Failed to initialize Conda"
    exit 1
}
conda activate pacbio-isoseq || {
    echo "ERROR: Failed to activate pacbio-isoseq environment"
    exit 1
}

#===============================#
# Configuration and Directories #
#===============================#
THREADS_PER_JOB=24  # Based on 192 total cores / 8 parallel jobs
NUM_PARALLEL_JOBS=$((192 / THREADS_PER_JOB))

cd /master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio

lima RNA.m84193_251117_232709_s3.hifi_reads.bc2015.bam barcodes.fasta /master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/split-barcode/RNA/RNA.demux.bam --hifi-preset ASYMMETRIC --split-named

lima DNA.m84193_251117_232709_s3.hifi_reads.bc2014.bam barcodes.fasta /master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/split-barcode/DNA/DNA.demux.bam --hifi-preset ASYMMETRIC --split-named

