#!/bin/bash

# launch run_clear_pair_locus.R on all 10 snATAC x 15 trait combinations (150 jobs)

WD=/working/lab_jonathb/lefeiW/projects/CLEAR_2026
SD=${WD}/jobs
<<<<<<< HEAD
OD=${WD}/output_19-03
PD=${WD}/processed_data
RSCRIPT=${WD}/R/run_clear_pair_locus_18-03.R
=======
OD=${WD}/output
PD=${WD}/processed_data
RSCRIPT=${WD}/R/run_clear_pair_locus.R
>>>>>>> 0466a99 (update)
GTF=/working/lab_jonathb/lefeiW/projects/CLEAR_data/gencode.v46.chr_patch_hapl_scaff.basic.annotation.gtf.gz
NPERM=1000
KLINEAGE=4

# make dirs
mkdir -p ${SD}/qsub ${OD} ${PD}

# ---- 12 snATAC SE inputs (name|path) ----
SE_INPUTS=(
  "breast_union_peak_330017|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_330017peak.rds"
  # "breast_final_union_531279|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/SE_breast_final_union_531279_PeakMatrix.hg38.rds"
  "pbmc|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/10X_pbmc_q0.01_243562peaks.rds"
  "hten_dcis|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/dcis_hten_100653peaks.rds"
  "zhang2021_mammary_ArchR|/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/output/snATAC_SE_library/zhang2021_mammary/zhang2021_mammary_ArchR_se.rds"
  "regner_mammary_cellline|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/regner_mammary_celllines_3000950peaks.rds"
  "kanemaru_heart|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/SE_kanemaru2023_heart_ArchR_PeakMatrix.hg38.rds"
  "kidney_25Cluster|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/kidney_25clusters_364335peaks.rds"
  "kidney_combinedClusters|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/kidney_combinedClusters_354332peaks.rds"
  "pancreas_140514|/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/pancreas_zhang2021_7CellTypes_140514peaks.rds"
  "breast_BRCA1_210591|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_BRCA1_210591peaks.rds"
  "breast_BRCA2_291074|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_BRCA2_291074peaks.rds"
  "breast_nonBRCA_264423|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/snATAC_ArchR_PeakMatrix/breast_nonBRCA_264423peaks.rds"
)

# ---- 15 GWAS trait inputs (name|path) ----
TRAIT_INPUTS=(
  # "bcac_overall_FM_Fachal_2020|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/BCAC_FM_GR.rds"
  # "t1d_santiago|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/t1d_santiago.rds"
  # "t2d|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/t2d_2024_suzuki_cs687__GCST90492734_CLEAR_credible_set_variants.gr.rds"
  # "JT_interval|/working/joint_projects/bc_risk_locus_multiomics/bc_risk_locus_funcannotation/openTargets_traits/JT_interval_CLEAR.gr.rds"
  # "DCIS|/mnt/lustre/working/lab_jonathb/lefeiW/projects/ATAC_BCAC/data/opentarget_trait/dcis_2023_kurki_cs8__FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC/dcis_2023_kurki_cs8__FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC_CLEAR_credible_set_variants.gr.rds"
  # "BCAC_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/BCAC_FineMapping.rds"
  # "LumA_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/LumA_FineMapping.rds"
  # "LumB_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/LumB_FineMapping.rds"
  # "LumB_HER2Neg_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/LumB_HER2Neg_FineMapping.rds"
  # "CIMBA_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/CIMBA_FineMapping.rds"
  # "HER2_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/HER2_FineMapping.rds"
  # "TripleNeg_susie|/mnt/lustre/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/TripleNeg_FineMapping.rds"
  # "alzheimer_2022_bellenguez_cs83__GCST90027158|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/alzheimer_2022_bellenguez_cs83__GCST90027158/alzheimer_2022_bellenguez_cs83__GCST90027158_CLEAR_credible_set_variants.gr.rds"
  # "bc_ovarian_prostate_pleiotropy_2016_kar_cs18__GCST010797|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/bc_ovarian_prostate_pleiotropy_2016_kar_cs18__GCST010797/bc_ovarian_prostate_pleiotropy_2016_kar_cs18__GCST010797_CLEAR_credible_set_variants.gr.rds"
  # "brain_and_nervous_cancer_2024_verma_cs52__GCST90479813|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/brain_and_nervous_cancer_2024_verma_cs52__GCST90479813/brain_and_nervous_cancer_2024_verma_cs52__GCST90479813_CLEAR_credible_set_variants.gr.rds"
  # "brain_cancer_2024_verma_cs84__GCST90479812|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/brain_cancer_2024_verma_cs84__GCST90479812/brain_cancer_2024_verma_cs84__GCST90479812_CLEAR_credible_set_variants.gr.rds"
  # "breast_carcinoma_2017_michailidou_cs210__GCST004988|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/breast_carcinoma_2017_michailidou_cs210__GCST004988/breast_carcinoma_2017_michailidou_cs210__GCST004988_CLEAR_credible_set_variants.gr.rds"
  # "colorectal_cancer_2024_tian_2024_cs70__GCST90435268|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/colorectal_cancer_2024_tian_2024_cs70__GCST90435268/colorectal_cancer_2024_tian_2024_cs70__GCST90435268_CLEAR_credible_set_variants.gr.rds"
  # "coronary_artery_disease_2022_aragam_cs250__GCST90132314|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/coronary_artery_disease_2022_aragam_cs250__GCST90132314/coronary_artery_disease_2022_aragam_cs250__GCST90132314_CLEAR_credible_set_variants.gr.rds"
  # "dcis_2023_kurki_cs8__FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/dcis_2023_kurki_cs8__FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC/dcis_2023_kurki_cs8__FINNGEN_R12_CD2_INSITU_BREAST_INTRADUCTAL_EXALLC_CLEAR_credible_set_variants.gr.rds"
  # "kidney_cancer_2024_purdue_2024_cs89__GCST90320054|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/kidney_cancer_2024_purdue_2024_cs89__GCST90320054/kidney_cancer_2024_purdue_2024_cs89__GCST90320054_CLEAR_credible_set_variants.gr.rds"
  # "mammographic_density_densearea_2020_sieh_cs10__GCST90011731|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/mammographic_density_densearea_2020_sieh_cs10__GCST90011731/mammographic_density_densearea_2020_sieh_cs10__GCST90011731_CLEAR_credible_set_variants.gr.rds"
  # "mammographic_density_nondensearea_2020_sieh_cs11__GCST90011732|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/mammographic_density_nondensearea_2020_sieh_cs11__GCST90011732/mammographic_density_nondensearea_2020_sieh_cs11__GCST90011732_CLEAR_credible_set_variants.gr.rds"
  # "monocyte_count_2020_chen_cs374__GCST90002344|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/monocyte_count_2020_chen_cs374__GCST90002344/monocyte_count_2020_chen_cs374__GCST90002344_CLEAR_credible_set_variants.gr.rds"
  # "ovarian_cancer_2022_dareng_cs27__GCST90016665|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/ovarian_cancer_2022_dareng_cs27__GCST90016665/ovarian_cancer_2022_dareng_cs27__GCST90016665_CLEAR_credible_set_variants.gr.rds"
  # "t2d_2024_suzuki_cs687__GCST90492734|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/t2d_2024_suzuki_cs687__GCST90492734/t2d_2024_suzuki_cs687__GCST90492734_CLEAR_credible_set_variants.gr.rds"
  # "tnbc_brca1_2020_zhang_cs18__GCST010100|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/tnbc_brca1_2020_zhang_cs18__GCST010100/tnbc_brca1_2020_zhang_cs18__GCST010100_CLEAR_credible_set_variants.gr.rds"
  # "heart_failure_2025_lee_cs165__GCST90455657|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/heart_failure_2025_lee_cs165__GCST90455657/heart_failure_2025_lee_cs165__GCST90455657_CLEAR_credible_set_variants.gr.rds"

  "platelet_count_2020_chen_cs1514__GCST90002357|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/platelet_count_2020_chen_cs1514__GCST90002357/platelet_count_2020_chen_cs1514__GCST90002357_CLEAR_credible_set_variants.gr.rds"
  "schizophrenia_2025_dang_cs247__GCST90503210|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/schizophrenia_2025_dang_cs247__GCST90503210/schizophrenia_2025_dang_cs247__GCST90503210_CLEAR_credible_set_variants.gr.rds"
  "sitting_height_2025_hu_cs686__GCST90310178|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/sitting_height_2025_hu_cs686__GCST90310178/sitting_height_2025_hu_cs686__GCST90310178_CLEAR_credible_set_variants.gr.rds"
  "systolic_blood_pressure_2024_keaton_cs1162__GCST90310294|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/systolic_blood_pressure_2024_keaton_cs1162__GCST90310294/systolic_blood_pressure_2024_keaton_cs1162__GCST90310294_CLEAR_credible_set_variants.gr.rds"
  "bmi_2022_huang_cs1071__GCST90255621|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/bmi_2022_huang_cs1071__GCST90255621/bmi_2022_huang_cs1071__GCST90255621_CLEAR_credible_set_variants.gr.rds"
  "educational_attainment_2022_okbay_cs1599__GCST90105038|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/educational_attainment_2022_okbay_cs1599__GCST90105038/educational_attainment_2022_okbay_cs1599__GCST90105038_CLEAR_credible_set_variants.gr.rds"
  "height_finngen_cs1903__FINNGEN_R12_HEIGHT_IRN|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/height_finngen_cs1903__FINNGEN_R12_HEIGHT_IRN/height_finngen_cs1903__FINNGEN_R12_HEIGHT_IRN_CLEAR_credible_set_variants.gr.rds"
  "ibd_2024_liu_cs293__GCST90292538|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/ibd_2024_liu_cs293__GCST90292538/ibd_2024_liu_cs293__GCST90292538_CLEAR_credible_set_variants.gr.rds"
  "diastolic_blood_pressure_2024_keaton_cs1108__GCST90310295|/working/lab_jonathb/lefeiW/projects/CLEAR_data/traits/opentarget_trait/diastolic_blood_pressure_2024_keaton_cs1108__GCST90310295/diastolic_blood_pressure_2024_keaton_cs1108__GCST90310295_CLEAR_credible_set_variants.gr.rds"

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
Rscript ${RSCRIPT} ${SE_NAME} ${SE_PATH} ${TRAIT_NAME} ${TRAIT_PATH} ${GTF} ${NPERM} ${KLINEAGE} ${PD}
EOF

    qsub -l ncpus=2,mem=15G,walltime=4:00:00 \
         -e ${OUTDIR}/clear_${PAIR_NAME}.ER \
         -o ${OUTDIR}/clear_${PAIR_NAME}.OU \
         -N clear_${PAIR_NAME} \
         ${RUNSCRIPT}

    echo "Submitted: ${PAIR_NAME}"

  done
done
