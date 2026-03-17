#!/bin/bash

# Submit ArchR batch runs per tissue using a TSV input file.

set -euo pipefail

WD=/working/lab_jonathb/lefeiW/projects/CLEAR_data/ArchR_batch_processed
SD=${WD}/jobs
OD=/working/lab_jonathb/lefeiW/projects/CLEAR_data/ArchR_batch_processed
RSCRIPT=${WD}/R/ArchR_streamline_batch.R
INPUT_TSV=${WD}/R/archr_batch_inputs.tsv
GENOME=hg38
THREADS=24
MINTSS=4
MINFRAGS=3200 # log10(3200) = 3.5 

mkdir -p ${SD}/qsub ${OD}

# Extract unique tissues (skip header); works with tabs or spaces.
TISSUES=$(awk 'NR>1 {print $1}' ${INPUT_TSV} | sort -u)

for TISSUE in ${TISSUES}; do
  OUTDIR=${OD}/${TISSUE}
  mkdir -p ${OUTDIR}

  RUNSCRIPT=${SD}/qsub/archr_${TISSUE}.sh

  cat > ${RUNSCRIPT} << EOF
#!/bin/bash
module load rstudio/R-4.3.1
cd ${OUTDIR}
Rscript ${RSCRIPT} ${INPUT_TSV} ${OD} ${GENOME} ${THREADS} ${MINTSS} ${MINFRAGS} ${TISSUE}
EOF

  qsub -l ncpus=24,mem=96G,walltime=8:00:00 \
       -e ${OUTDIR}/archr_${TISSUE}.ER \
       -o ${OUTDIR}/archr_${TISSUE}.OU \
       -N archr_${TISSUE} \
       ${RUNSCRIPT}

  echo "Submitted: ${TISSUE}"

done
