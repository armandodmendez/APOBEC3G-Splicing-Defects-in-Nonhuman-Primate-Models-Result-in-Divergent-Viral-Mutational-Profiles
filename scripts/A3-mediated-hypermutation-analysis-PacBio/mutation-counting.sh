#!/bin/bash
#SBATCH --job-name=mutation-counting
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-20%10
#SBATCH --chdir=/master/amendez/slurm-scripts
#SBATCH --output=logs/mutation-counting-%A_%a.out
#SBATCH --error=logs/mutation-counting-%A_%a.err

set -euo pipefail
echo "=== Task $SLURM_ARRAY_TASK_ID ($(date)) ==="

source /master/amendez/miniconda3/etc/profile.d/conda.sh
conda activate pacbio-isoseq

# Get sorted BAM lists (alphabetical)
mapfile -t HIV_BAMS < <(ls -1 /master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/mutational-calling/HIV/*.bam 2>/dev/null)
mapfile -t SIV_BAMS < <(ls -1 /master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/mutational-calling/SIV/*.bam 2>/dev/null)

TOTAL_BAMS=${#HIV_BAMS[@]}
TOTAL_BAMS=$((TOTAL_BAMS + ${#SIV_BAMS[@]}))

echo "Found ${#HIV_BAMS[@]} HIV BAMs, ${#SIV_BAMS[@]} SIV BAMs (total $TOTAL_BAMS)"

if [ $SLURM_ARRAY_TASK_ID -gt $TOTAL_BAMS ]; then
    echo "No more BAMs for task $SLURM_ARRAY_TASK_ID"
    exit 0
fi

if [ $SLURM_ARRAY_TASK_ID -le ${#HIV_BAMS[@]} ]; then
    # HIV
    BAM_PATH="${HIV_BAMS[$((SLURM_ARRAY_TASK_ID-1))]}"
    REF_PATH="/master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/reference-genomes/pJH048_i3_HIV-NFL.fasta"
    RESULT_DIR="/master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/mutational-calling/HIV/res/"
else
    # SIV  
    IDX=$((SLURM_ARRAY_TASK_ID - ${#HIV_BAMS[@]} - 1))
    BAM_PATH="${SIV_BAMS[$IDX]}"
    REF_PATH="/master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/reference-genomes/SIV239-SpX_i22-NFL.fasta"
    RESULT_DIR="/master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio/mutational-calling/SIV/res/"
fi

echo "Task $SLURM_ARRAY_TASK_ID -> $BAM_PATH"
ls -lh "$REF_PATH" "$BAM_PATH"
mkdir -p "$RESULT_DIR"
export REF_PATH BAM_PATH RESULT_DIR
python /master/amendez/slurm-scripts/python-scripts/mutation_counting2.py
