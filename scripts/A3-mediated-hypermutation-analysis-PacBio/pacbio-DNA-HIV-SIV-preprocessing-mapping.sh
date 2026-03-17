#!/bin/bash
#$ -N pacbio-DNA-HIV-SIV-preprocessing-mapping.sh
#$ -cwd
#$ -pe smp 24
#$ -q all.q
#$ -wd /master/amendez/qsubs

# Version: 1.2
# Author: Armando Mendez <amendez@txbiomed.org>
# Updated on: 2025-11-24
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
cd /master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio

REF_HIV=/master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/reference-genomes/pJH048_i3_HIV-NFL.fasta     
REF_SIV=/master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/reference-genomes/SIV239-SpX_i22-NFL.fasta

THREADS=24

############################
# DNA: with HIV/SIV subdirectories, mapped output, and samtools -F 2308
############################

cd /master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/split-barcode/DNA/trimmed || exit 1

for virus in HIV SIV; do
  cd "$virus" || { echo "Failed to enter directory $virus"; exit 1; }
  mkdir -p mapped

  for bam in *.trimmed.bam; do
    [ -e "$bam" ] || { echo "No BAM files found for $virus"; break; }
    sample="${bam%.trimmed.bam}"

    if [ ! -f "${bam}.pbi" ]; then
      pbindex "$bam"
    fi

    bam2fastq -o "${sample}" "$bam"

    if [ "$virus" = "HIV" ]; then
      REF="$REF_HIV"
    else
      REF="$REF_SIV"
    fi

    minimap2 -t "${THREADS}" -ax map-hifi "$REF" "${sample}.fastq.gz" > "mapped/${sample}.sam"

    samtools view -b -F 2308 "mapped/${sample}.sam" | samtools sort -o "mapped/${sample}.aligned.bam" -
    samtools index "mapped/${sample}.aligned.bam"

    samtools fasta "mapped/${sample}.aligned.bam" > "${sample}_queries.fasta"

    rm "mapped/${sample}.sam"
  done

  cd ..
done

