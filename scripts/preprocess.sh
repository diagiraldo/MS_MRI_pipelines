#!/usr/bin/zsh

# Preprocess nifti images: denoise, brain extraction, N4
# Results are saver in folder of processed data
# Diana Giraldo, Nov 2022
# It requires: ANTs, HD-BET

############################################
# Inputs: 
# Raw Nifti image
RAW_IM=${1}
# Folder for processed MRI
PRO_DIR=${2}
# ANTs directory (bin)
ANTS_DIR=${3}
############################################

export ANTSPATH=${ANTS_DIR}

# Image basename
IM_BN=$(basename ${RAW_IM} | sed 's/.nii.gz//')

# .json file with DICOM info
JSFILE=$( echo ${RAW_IM} | sed 's/.nii.gz//').json

# Output directory
TMP_PTH=$(dirname ${RAW_IM} | rev | cut -d"/" --fields=1,2,3 | rev )
OUT_DIR=${PRO_DIR}/${TMP_PTH}
mkdir -p ${OUT_DIR}

# Denoise
DenoiseImage -d 3 -n Rician -i ${RAW_IM} -o ${OUT_DIR}/${IM_BN}_dn.nii.gz
# Calculate absolute value to remove negatives
ImageMath 3 ${OUT_DIR}/${IM_BN}_dnabs.nii.gz abs ${OUT_DIR}/${IM_BN}_dn.nii.gz
mv ${OUT_DIR}/${IM_BN}_dnabs.nii.gz ${OUT_DIR}/${IM_BN}_dn.nii.gz

# Brain Extraction
hd-bet -i ${OUT_DIR}/${IM_BN}_dn.nii.gz -o ${OUT_DIR}/${IM_BN}_bet.nii.gz -device cpu -mode fast -tta 0 > /dev/null
rm ${OUT_DIR}/${IM_BN}_bet.nii.gz
mv ${OUT_DIR}/${IM_BN}_bet_mask.nii.gz ${OUT_DIR}/${IM_BN}_brainmask.nii.gz

# Biasfield correction N4
N4BiasFieldCorrection -d 3 -i ${OUT_DIR}/${IM_BN}_dn.nii.gz -o ${OUT_DIR}/${IM_BN}_preproc.nii.gz -x ${OUT_DIR}/${IM_BN}_brainmask.nii.gz

# Remove denoised image
rm ${OUT_DIR}/${IM_BN}_dn.nii.gz

# Copy .json
cp ${JSFILE} ${OUT_DIR}/.

############################################
# Output:
# preprocessed image 
echo ${OUT_DIR}/${IM_BN}_preproc.nii.gz
############################################