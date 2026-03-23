#!/bin/bash

# Submit CLEAR trait downsampling sensitivity jobs (one qsub per schedule).
# Schedules: mild, medium, strong
# Repeats per schedule: 3

set -euo pipefail

# Path
WD=/working/lab_jonathb/lefeiW/projects/CLEAR_data
RSCRIPT=${WD}/traits/downsample_traits_signal.R
INPUT_DIR=/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait
OUT_ROOT=/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/downsampling_test
QSUB_DIR=${OUT_ROOT}/qsub_jobs

# ---------- Parameters ----------
FILE_PATTERN='_CLEAR_credible_set_variants\.gr\.rds$'
TARGET_FRACTION=0.5
MIN_PER_SIGNAL=20 # The Min and Max are override schedule. The hard cap
MAX_PER_SIGNAL=200 # if fractioned over 200, randomly sample 200 SNPs
BASE_SEED=1
N_REPEATS=3
SCHEDULES=(mild medium strong)

# ---------- Resources ----------
NCPUS=2
MEM=16G
WALLTIME=06:00:00

mkdir -p "${QSUB_DIR}" "${OUT_ROOT}"

for SCHED in "${SCHEDULES[@]}"; do
  OUT_DIR="${OUT_ROOT}/${SCHED}"
  mkdir -p "${OUT_DIR}"

  RUNSCRIPT="${QSUB_DIR}/run_downsample_${SCHED}.sh"

  cat > "${RUNSCRIPT}" << EOF
#!/bin/bash
set -euo pipefail
module load rstudio/R-4.3.1
cd "${OUT_DIR}"
Rscript "${RSCRIPT}" \
  "${INPUT_DIR}" \
  "${OUT_DIR}" \
  "${FILE_PATTERN}" \
  "${TARGET_FRACTION}" \
  "${MIN_PER_SIGNAL}" \
  "${MAX_PER_SIGNAL}" \
  "${BASE_SEED}" \
  "${SCHED}" \
  "${N_REPEATS}"
EOF

  chmod +x "${RUNSCRIPT}"

  qsub \
    -l ncpus=${NCPUS},mem=${MEM},walltime=${WALLTIME} \
    -N "ds_${SCHED}" \
    -o "${OUT_DIR}/downsample_${SCHED}.OU" \
    -e "${OUT_DIR}/downsample_${SCHED}.ER" \
    "${RUNSCRIPT}"

  echo "Submitted schedule: ${SCHED}"
  echo "  Output dir: ${OUT_DIR}"
done

echo "All sensitivity jobs submitted."
