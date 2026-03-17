#!/bin/bash

# launch run_clear_pair_locus.R on all 10 snATAC x 15 trait combinations (150 jobs)

WD=/working/lab_jonathb/lefeiW/projects/CLEAR_2026
SD=${WD}/jobs
OD=${WD}/output
RSCRIPT=${WD}/R/run_clear_pair_locus.R
GTF=/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz
NPERM=100
KLINEAGE=4

# make dirs
mkdir -p ${SD}/qsub ${OD}

# ---- 10 snATAC SE inputs (name|path) ----
SE_INPUTS=(
  "breast_union_peak_330017|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_330017peak.rds"
  "breast_final_union_531279|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/SE_breast_final_union_531279_PeakMatrix.hg38.rds"
  "pbmc|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/10X_pbmc_q0.01_243562peaks.rds"
  "hten_dcis|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/dcis_hten_100653peaks.rds"
  "zhang2021_mammary_ArchR|/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/output/snATAC_SE_library/zhang2021_mammary/zhang2021_mammary_ArchR_se.rds"
  "regner_mammary_cellline|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/regner_mammary_celllines_3000950peaks.rds"
  "kanemaru_heart|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/SE_kanemaru2023_heart_ArchR_PeakMatrix.hg38.rds"
  "kidney_25Cluster|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/kidney_25clusters_364335peaks.rds"
  "kidney_combinedClusters|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/kidney_combinedClusters_354332peaks.rds"
  "pancreas_140514|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/pancreas_zhang2021_7CellTypes_140514peaks.rds"
)

# ---- 15 GWAS trait inputs (name|path) ----
TRAIT_INPUTS=(
  "bcac_overall_opentarget|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/breast_carcinoma_2017_michailidou_cs210__GCST004988_CLEAR_credible_set_variants.gr.rds"
  "bcac_overall_FM_Fachal_2020|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/BCAC_FM_GR.rds"
  "t1d_santiago|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/t1d_santiago.rds"
  "t2d|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/t2d_2024_suzuki_cs687__GCST90492734_CLEAR_credible_set_variants.gr.rds"
  "heart_failure|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/heart_failure_2025_lee_cs165__GCST90455657_CLEAR_credible_set_variants.gr.rds"
  "JT_interval|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/JT_interval_CLEAR.gr.rds"
  "height|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/height_finngen_cs1903__FINNGEN_R12_HEIGHT_IRN_CLEAR_credible_set_variants.gr.rds"
  "mammographic_density_nondensearea|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/mammographic_density_nondensearea_2020_sieh_cs11__GCST90011732_CLEAR_credible_set_variants.gr.rds"
  "mammographic_density_densearea|/mnt/lustre/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/data/opentarget_trait/mammographic_density_densearea_2020_sieh_cs10__GCST90011731/mammographic_density_densearea_2020_sieh_cs10__GCST90011731_CLEAR_credible_set_variants.gr.rds"
  "monocyte_counts|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/monocyte_count_2020_chen_cs374__GCST90002344_CLEAR_credible_set_variants.gr.rds"
  "platelet_counts|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/platelet_count_2020_chen_cs1514__GCST90002357_CLEAR_credible_set_variants.gr.rds"
  "ovarian_cancer|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/ovarian_cancer_2022_dareng_cs27__GCST90016665_CLEAR_credible_set_variants.gr.rds"
  "kidney_cancer|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/kidney_cancer_2024_purdue_2024_cs89__GCST90320054_CLEAR_credible_set_variants.gr.rds"
  "colorectal_cancer|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/colorectal_cancer_2024_tian_2024_cs70__GCST90435268_CLEAR_credible_set_variants.gr.rds"
  "brain_cancer|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/brain_cancer_2024_verma_cs84__GCST90479812_CLEAR_credible_set_variants.gr.rds"
  "alzheimers|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/alzheimer_2022_bellenguez_cs83__GCST90027158_CLEAR_credible_set_variants.gr.rds"
  "DCIS|/mnt/lustre/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/data/opentarget_trait/dcis_2023_kurki_cs8__FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC/dcis_2023_kurki_cs8__FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC_CLEAR_credible_set_variants.gr.rds"
  "BCAC_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/BCAC_FineMapping.rds" 
  "LumA_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/LumA_FineMapping.rds"
  "LumB_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/LumB_FineMapping.rds"
  "LumB_HER2Neg_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/LumB_HER2Neg_FineMapping.rds"
  "CIMBA_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/CIMBA_FineMapping.rds"
  "HER2_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/HER2_FineMapping.rds"
  "TripleNeg_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/TripleNeg_FineMapping.rds"

)

# ---- Loop over all SE x trait pairs ----
for SE_ENTRY in "${SE_INPUTS[@]}"; do

  SE_NAME=${SE_ENTRY%%|*}
  SE_PATH=${SE_ENTRY#*|}

  for TRAIT_ENTRY in "${TRAIT_INPUTS[@]}"; do

    TRAIT_NAME=${TRAIT_ENTRY%%|*}
    TRAIT_PATH=${TRAIT_ENTRY#*|}

    PAIR_NAME="${SE_NAME}_x_${TRAIT_NAME}"
    OUTDIR=${OD}/${PAIR_NAME}
    mkdir -p ${OUTDIR}

    RUNSCRIPT=${SD}/qsub/clear_${PAIR_NAME}.sh

    cat > ${RUNSCRIPT} << EOF
#!/bin/bash
module load rstudio/R-4.3.1
cd ${OUTDIR}
Rscript ${RSCRIPT} ${SE_NAME} ${SE_PATH} ${TRAIT_NAME} ${TRAIT_PATH} ${GTF} ${NPERM} ${KLINEAGE}
EOF

    qsub -l ncpus=2,mem=15G,walltime=1:00:00 \
         -e ${OUTDIR}/clear_${PAIR_NAME}.ER \
         -o ${OUTDIR}/clear_${PAIR_NAME}.OU \
         -N clear_${PAIR_NAME} \
         ${RUNSCRIPT}

    echo "Submitted: ${PAIR_NAME}"

  done
done
