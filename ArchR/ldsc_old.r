#!/bin/bash

# Note to check genome build of SumStats/annotations/reference 

WD=/working/lab_jonathb/lefeiW/projects/sc_eQTL

# hg38 liftover to hg19 peaks

# FACS sorted breast cell linage ATAC bed fiels
PEAK_DIR=${WD}/data/bc_annotation/cts_specific_ATAC_HiChIP

module load ldsc/20190815
module load bedtools/2.31.1

# Loop through each BED file in the peak directory
for bed_file in ${PEAK_DIR}/*_hg19.bed; do
  # Extract the tissue type from the file name (first part before the dot) NEED to check file names
  tissue_type=$(basename "${bed_file}" | cut -d '.' -f 1)

  # Set annotation name
  annotation_name="bc_annotation_CTS_ATAC_HiChIP"
  OUTDIR=${WD}/results/cts_specific_ATAC_HiChIP/${annotation_name}_${tissue_type}_partitioned

  # Set other variables
  bed=${WD}/data/s_LDSC_annot/5kb_HiChIP_FACS_${tissue_type}/bed
  annot=${WD}/data/s_LDSC_annot/5kb_HiChIP_FACS_${tissue_type}/annot
  
  KG_DIR=${WD}/data/1KG_LD_Phase3
  GWAS_DIR=${WD}/results/Cleaned_BC_subtype_GWAS                       # hg19
  GWAS_FILE=( ${GWAS_DIR}/*_GWAS_summary.txt  )
  SCRIPTDIR=${WD}/script/${annotation_name}_${tissue_type}_partitioned
  

  # Create directories 
  mkdir -p ${bed} ${annot} ${SCRIPTDIR}/{log,qsub} ${OUTDIR}

  # Copy the bed file to the working directory
  cp ${bed_file} ${bed}/${annotation_name}.sorted.bed

  # Split BED annotation into chr1..22 and remove 'chr' prefix # Check original bed file format
  for chr in {1..22}; do
    awk -v chr="chr$chr" '$1 == chr' ${bed}/${annotation_name}.sorted.bed > ${bed}/${annotation_name}.chr$chr.bed
  done

  # Making the annotation matrix file required by LDSC
  for chr in {1..22}; do
    python /software/ldsc/ldsc-20190815/make_annot.py \
      --bed-file ${bed}/${annotation_name}.chr$chr.bed \
      --bimfile ${KG_DIR}/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim \
      --annot-file ${annot}/${annotation_name}.chr$chr.annot.gz
  done

  # Recalculate the LD score for the annotation using qsub scripts and capture job IDs
  ldscore_job_ids=()
  for chr in {1..22}; do
    RUNSCRIPT=${SCRIPTDIR}/qsub/ldsc_recalc_ldscore.chr$chr.sh

    cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load ldsc/20190815

python /software/ldsc/ldsc-20190815/ldsc.py --l2 \
  --bfile ${KG_DIR}/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr \
  --ld-wind-cm 1 \
  --annot ${annot}/${annotation_name}.chr$chr.annot.gz \
  --thin-annot \
  --out ${annot}/${annotation_name}.chr$chr \
  --print-snps ${KG_DIR}/hm3_no_MHC.list.txt

EOF

    job_id=$(qsub -l ncpus=4,mem=16G,walltime=01:00:00 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log -N ldscore_chr$chr ${RUNSCRIPT})
    ldscore_job_ids+=($job_id)
  done

  # Convert job IDs array to colon-separated list
  ldscore_job_ids_str=$(IFS=:; echo "${ldscore_job_ids[*]}")

  # Final step in partitioned heritability
  for f in "${GWAS_FILE[@]}" ; do
    FILE=$( basename "${f%_GWAS_summary.txt}" )
    RUNSCRIPT=${SCRIPTDIR}/qsub/${annotation_name}_${tissue_type}_partitioned.${FILE}.sh

    cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load ldsc/20190815

python /software/ldsc/ldsc-20190815/ldsc.py \
  --h2 ${f} \
  --ref-ld-chr ${KG_DIR}/baselineLDV2.2/baselineLD.,${annot}/${annotation_name}.chr \
  --w-ld-chr ${KG_DIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --overlap-annot \
  --frqfile-chr ${KG_DIR}/1000G_Phase3_frq/1000G.EUR.QC. \
  --out ${OUTDIR}/${FILE}_${annotation_name} \
  --print-coefficients

EOF

    qsub -W depend=afterok:$ldscore_job_ids_str -l ncpus=2,mem=8G,walltime=01:00:00 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log -N parti_h2_${FILE} ${RUNSCRIPT}
  done

done

# Run the script
# bash make_annot.sh > make_annot.log 2>&1 &
