#!/bin/bash

#$ -N combine_hypermut3_final.sh
#$ -cwd
#$ -pe smp 8
#$ -q all.q
#$ -wd /master/amendez/qsubs

# Version: 1.4
# Author: Armando Mendez <amendez@txbiomed.org>
# Updated on: 2025-12-08
# combine_hypermut3_final.sh

#=================#
# Initialize Conda #
#=================#
source /master/amendez/miniconda3/etc/profile.d/conda.sh || {
    echo "ERROR: Failed to initialize Conda"
    exit 1
}
conda activate hypermut3 || {
    echo "ERROR: Failed to activate hypermut3 environment"
    exit 1
}



set -euo pipefail

ROOT="/master/amendez/Experimental-ADM/HIV-SIV-Hypermutation-Assay-PacBio"
OUT_DIR="${ROOT}/hypermut3-combined"

# 1) Make sure the list of all *_summary.csv files exists
find "${ROOT}" -type f -name "*_hypermut3_*_summary.csv" > "${OUT_DIR}/hypermut3_csv_list.txt"

# 2) Run the final combiner with BOTH required arguments
python3 "${OUT_DIR}/combine_hypermut3_final.py" \
  --list "${OUT_DIR}/hypermut3_csv_list.txt" \
  --out  "${OUT_DIR}/master_hypermut3_final.csv"

