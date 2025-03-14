#!/usr/bin/zsh

# Example Process session with fast HR FLAIR (C), and sT1 MRI
# Diana Giraldo, Nov 2022
# Last update: March 2025

###################################################
# Requirements:

# ANTs directory
ANTS_DIR=/opt/ANTs/bin

# SPM directory
SPM_DIR=/home/vlab/spm12/

# Processing repo directory
SCR_DIR=/home/vlab/MS_MRI_pipelines/scripts
# LST auxiliary functions
SEGF_DIR=${SCR_DIR}/LSTfunctions/

###################################################
# Image directory containing raw T1 and FLAIR images
CASE=0229128
DATE=20160818
RAW_IM_DIR=/home/vlab/MS_proj/MS_MRI/sub-${CASE}/ses-${DATE}/anat

# Directory for processed images
IM_DIR=/home/vlab/testpipeline

###################################################################################
# Pre-process T1 and FLAIR images in IM_DIR
for SEQ in sT1 Flair_fast; 
do
    RAW_IM=$(ls ${RAW_IM_DIR}/*${SEQ}*.nii* | head -n 1 )
    zsh ${SCR_DIR}/preprocess.sh ${RAW_IM} ${IM_DIR} ${ANTS_DIR}
done

###################################################################################
# Align sT1 to fast FLAIR with rigid registration
REF_SEQ=Flair_fast
MOV_SEQ=sT1

REF_IM=$(ls ${IM_DIR}/*${REF_SEQ}*_preproc.nii* | head -n 1 )
REF_MASK=$(ls ${IM_DIR}/*${REF_SEQ}*_brainmask.nii* | head -n 1 )
MOV_IM=$(ls ${IM_DIR}/*${MOV_SEQ}*_preproc.nii* | head -n 1 )
MOV_MASK=$(ls ${IM_DIR}/*${MOV_SEQ}*_brainmask.nii* | head -n 1 )
OUT_PRE=${IM_DIR}/rigreg_${MOV_SEQ}_to_${REF_SEQ}_

antsRegistration --dimensionality 3 --output \[ ${OUT_PRE} \] \
    --collapse-output-transforms 1 \
    --interpolation BSpline \
    --initial-moving-transform \[ ${REF_IM},${MOV_IM},1 \] \
    --metric MI\[ ${REF_IM},${MOV_IM},1,32,Regular,0.25 \] \
    --transform Rigid\[ 0.1 \] \
    --convergence \[ 1000x500x250x0,1e-6,10 \] \
    --smoothing-sigmas 3x2x1x0vox \
    --shrink-factors 8x4x2x1 \
    --use-histogram-matching 0 \
    --winsorize-image-intensities \[ 0.005,0.995 \] \
    --masks \[ ${REF_MASK},${MOV_MASK} \] \
    --float 0 \
    --verbose 0

antsApplyTransforms --dimensionality 3 \
    --input ${MOV_IM} \
    --reference-image ${REF_IM} \
    --output ${OUT_PRE}transformed.nii.gz \
    --interpolation BSpline \
    --transform ${OUT_PRE}0GenericAffine.mat
mrcalc ${OUT_PRE}transformed.nii.gz 0 -lt 0 ${OUT_PRE}transformed.nii.gz -if ${OUT_PRE}transformed.nii.gz -force -quiet

###################################################################################
# Segmentation

FLAIR_IM=${IM_DIR}/HR_FLAIR_seginput.nii.gz
T1_IM=${IM_DIR}/HR_T1_seginput.nii.gz

# Remove small values so segmentation doesn't fail
epsilon=0.001
mrcalc ${IM_DIR}/rigreg_sT1_to_Flair_fast_transformed.nii.gz 0 ${epsilon} -replace ${T1_IM} -force -quiet
mrcalc $(ls ${IM_DIR}/*Flair_fast*_preproc.nii* | head -n 1 ) 0 ${epsilon} -replace ${FLAIR_IM} -force -quiet

# SAMSEG
OUT_DIR=${IM_DIR}/samseg
mkdir -p ${OUT_DIR}
run_samseg --input ${T1_IM} ${FLAIR_IM} --output ${OUT_DIR} \
--pallidum-separate --lesion --lesion-mask-pattern 0 1  \
--random-seed 22  --threads 8
rm ${OUT_DIR}/mode*_bias_*.mgz ${OUT_DIR}/template_coregistered.mgz
# Check unknowns within brain
BR_MASK=$(${IM_DIR}/*Flair_fast*_brainmask.nii* | head -n 1)
mrcalc ${OUT_DIR}/seg.mgz 0 -eq ${BR_MASK} -mult ${OUT_DIR}/unknownswithinbrain.nii.gz -datatype bit -force -quiet
mrstats ${OUT_DIR}/unknownswithinbrain.nii.gz -mask ${OUT_DIR}/unknownswithinbrain.nii.gz -quiet -output count > ${OUT_DIR}/count_unknownswithinbrain.txt

# LST
# threshold for segmentation
thLST=0.1
OUT_DIR=${IM_DIR}/LST
mkdir -p ${OUT_DIR}
INNII=${OUT_DIR}/input.nii
mrconvert ${FLAIR_IM} ${INNII} -quiet -force
PLES=${OUT_DIR}/ples_lpa_minput.nii
matlab -nodisplay -r "addpath('$SPM_DIR'); addpath('$SEGF_DIR'); cd '$OUT_DIR'; lst_lpa('$INNII', 0); lst_lpa_voi('$PLES', '$thLST'); exit"
mrconvert ${PLES} ${OUT_DIR}/ples_lpa.nii.gz 
mv ${OUT_DIR}/LST_tlv_${thLST}_*.csv ${OUT_DIR}/LST_lpa_${thLST}.csv
rm ${OUT_DIR}/input.nii ${OUT_DIR}/minput.nii ${OUT_DIR}/LST_lpa_minput.mat ${PLES}

# Remove adjusted images
rm ${IM_DIR}/HR_FLAIR_seginput.nii.gz ${IM_DIR}/HR_T1_seginput.nii.gz
