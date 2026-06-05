#!/bin/bash

#PBS -N OT_GCST_download
#PBS -l nodes=1:ppn=1,mem=16gb,walltime=6:00:00
#PBS -j oe
#PBS -o /working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/API_opentargets/log/OT_GCST_download.log

set -euo pipefail

module load R/4.5.0
source ~/.proxy

LOG_DIR="/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/API_opentargets/log"
R_SCRIPT="/working/lab_jonathb/lefeiW/projects/CLEAR_working/CLEAR_working/OpenTargets_API/download_GCST_all.r"

mkdir -p "${LOG_DIR}"

echo "Host: $(hostname)"
echo "Started: $(date)"
echo "Script: ${R_SCRIPT}"
Rscript "${R_SCRIPT}"
echo "Finished: $(date)"
