#!/usr/bin/zsh

# Example Process session with LR FLAIR (A, B, D) Using PRETTIER
# Diana Giraldo, Nov 2022
# Last update: Feb 2025

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

# PRETTIER directory for MRI super-resolution
PRETTIER_DIR=/home/vlab/PRETTIER

###################################################
# Image directory containing raw T1 and FLAIR images
CASE=0030403
DATE=20120507
IM_DIR=/home/vlab/MS_proj/MS_MRI/sub-${CASE}/ses-${DATE}/anat

###################################################################################
# Pre-process T1 and FLAIR images in IM_DIR
for RAW_IM in $(ls ${IM_DIR}/*(sT1|[Ff][Ll][Aa][Ii][Rr])*.nii*);
do
    zsh ${SCR_DIR}/preprocess.sh ${RAW_IM} ${IM_DIR} ${ANTS_DIR}
done

###################################################################################
# Super-Resolution: Use available LR FLAIR to obtain a HR image
slcth=2

# Copy LR FLAIR (and BRAIN masks) to subfolders
FL_DIR=${IM_DIR}/LR_FLAIR_preproc
BM_DIR=${IM_DIR}/LR_FLAIR_masks
mkdir -p ${FL_DIR} ${BM_DIR}
for IM in $(ls ${IM_DIR}/*[Ff][Ll][Aa][Ii][Rr]*preproc.nii.gz); 
do
    SLC=$( mrinfo ${IM} -spacing -config RealignTransform 0 | cut -d" " -f3 )
    if [[ ${SLC} > ${slcth} ]]; then
        cp ${IM} ${FL_DIR}/.
        cp ${IM%_preproc.nii.gz}_brainmask.nii.gz ${BM_DIR}/.
    fi
done

# Histogram matching
HM_DIR=${IM_DIR}/LR_FLAIR_hmatch
zsh ${SCR_DIR}/histmatch_folder.sh ${FL_DIR} ${BM_DIR} ${HM_DIR}
rm -r ${FL_DIR}

# Apply PRETTIER to each LR and obtain a brain mask for each one
SR_DIR=${IM_DIR}/LR_FLAIR_prettier
SRM_DIR=${IM_DIR}/LR_FLAIR_prettier_masks
mkdir -p ${SR_DIR} ${SRM_DIR}
for IM in $(ls ${HM_DIR}/*.nii.gz);
do
    IM_BN=$(basename ${IM} | sed 's/_preproc.nii.gz//')

    if [[ ! -f ${SR_DIR}/${IM_BN}_prettierEDSR.nii.gz ]]; 
    then
        ${PRETTIER_DIR}/prettier_mri.py --input ${IM} --model-name EDSR --output ${SR_DIR}/${IM_BN}_prettierEDSR.nii.gz --batch-size 12 --gpu-id 3 --quiet
        hd-bet -i ${SR_DIR}/${IM_BN}_prettierEDSR.nii.gz -o ${SRM_DIR}/${IM_BN}_bet.nii.gz -device cpu -mode fast -tta 0 > /dev/null
        rm ${SRM_DIR}/${IM_BN}_bet.nii.gz
        mv ${SRM_DIR}/${IM_BN}_bet_mask.nii.gz ${SRM_DIR}/${IM_BN}_brainmask.nii.gz
    else
        echo "${SR_DIR}/${IM_BN}_prettierEDSR.nii.gz already exists"
    fi
done

# Create a 1mm grid 
HR_grid=${IM_DIR}/HRgrid_1mm.nii.gz
zsh ${SCR_DIR}/get_HRgrid.sh ${BM_DIR} 1 ${HR_grid}

# Align and combine PRETTIER outputs
export PYTHONPATH=/home/vlab/mrtrix3/lib:$PYTHONPATH
${SCR_DIR}/align_combine.py $(ls ${SR_DIR}/*.nii.gz) ${HR_grid} ${IM_DIR}/HR_FLAIR_combined_prettierEDSR.nii.gz -masks $(ls ${SRM_DIR}/*.nii.gz) -iter 2 -force

###################################################################################
# Segmentation

# Prepare for segmentation
HRFLAIR=${IM_DIR}/HR_FLAIR_combined_prettierEDSR.nii.gz

# Create a brain mask for the HR image
HRFLMASK=${IM_DIR}/HR_FLAIR_prettierEDSR_brainmask.nii.gz 
if [[ ! -f ${HRFLMASK} ]]; then
    hd-bet -i ${HRFLAIR} -o ${IM_DIR}/HR_FLAIR_bet.nii.gz -device cpu -mode fast -tta 0 > /dev/null
    rm ${IM_DIR}/HR_FLAIR_bet.nii.gz
    mv ${IM_DIR}/HR_FLAIR_bet_mask.nii.gz ${HRFLMASK}
fi

# Remove small values so segmentation doesn't fail
epsilon=0.01    

FLAIR_IM=${IM_DIR}/HR_FLAIR_seginput.nii.gz
mrcalc ${HRFLAIR} 0 ${epsilon} -replace ${FLAIR_IM} -force -quiet

# LST
# threshold for segmentation
thLST=0.1
OUT_DIR=${IM_DIR}/prettier_LST
mkdir -p ${OUT_DIR}
INNII=${OUT_DIR}/input.nii
mrconvert ${FLAIR_IM} ${INNII} -quiet -force
PLES=${OUT_DIR}/ples_lpa_minput.nii
matlab -nodisplay -r "addpath('$SPM_DIR'); addpath('$SEGF_DIR'); cd '$OUT_DIR'; lst_lpa('$INNII', 0); lst_lpa_voi('$PLES', '$thLST'); exit"
mrconvert ${PLES} ${OUT_DIR}/ples_lpa.nii.gz -force -quiet
mv ${OUT_DIR}/LST_tlv_${thLST}_*.csv ${OUT_DIR}/LST_lpa_${thLST}.csv
rm ${OUT_DIR}/input.nii ${OUT_DIR}/minput.nii ${OUT_DIR}/LST_lpa_minput.mat ${PLES}

# SAMSEG
OUT_DIR=${IM_DIR}/prettier_samseg
mkdir -p ${OUT_DIR}
run_samseg --input ${FLAIR_IM} --output ${OUT_DIR} \
--pallidum-separate --lesion --lesion-mask-pattern 1 \
--random-seed 22 --threads 8
rm ${OUT_DIR}/mode*_bias_*.mgz ${OUT_DIR}/template_coregistered.mgz
# Check unknowns within brain
mrcalc ${OUT_DIR}/seg.mgz 0 -eq ${HRFLMASK} -mult ${OUT_DIR}/unknownswithinbrain.nii.gz -datatype bit -force -quiet
mrstats ${OUT_DIR}/unknownswithinbrain.nii.gz -mask ${OUT_DIR}/unknownswithinbrain.nii.gz -quiet -output count > ${OUT_DIR}/count_unknownswithinbrain.txt