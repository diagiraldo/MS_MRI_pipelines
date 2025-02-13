#!/usr/bin/zsh

# Intensity normalization across multiple images
# Diana Giraldo, July 2022
# It uses: 'mrhistmatch' from mrtrix3

############################################
# Inputs: 
# Folder with preprocessed images 
IN_DIR=${1}
# Folder with masks
INMSK_DIR=${2}
# Output folder to save intensity normalized images
OUT_DIR=${3}
############################################


# Make temporary directory
TMP_DIR=$(mktemp -d)

nLR=$( ls ${IN_DIR}/* | wc -l )
echo "-------------------------------------------------------------"
echo "Number of input images: ${nLR}"
echo "-------------------------------------------------------------"

IMG_DIR=${TMP_DIR}/img
MASK_DIR=${TMP_DIR}/mask
HMATCH_DIR=${TMP_DIR}/hmatch
mkdir -p ${IMG_DIR} ${MASK_DIR} ${HMATCH_DIR}

# Copy images to temporary directory
cp ${IN_DIR}/*.nii*  ${IMG_DIR}/.
cp ${INMSK_DIR}/*.nii* ${MASK_DIR}/.

# Select transversal/axial as reference 
nTRA=$(find ${IMG_DIR} -name "sub-*_ses-*_*TRA_preproc.nii*" -type f | wc -l)
if [ $nTRA -gt 0 ]; then
    IM_REF=$(ls ${IMG_DIR}/sub-*_ses-*_*TRA_preproc.nii* | head -n 1 )
else
    IM_REF=$(ls ${IMG_DIR}/sub-*_ses-*_*preproc.nii* | head -n 1 )
fi

FILENAME=${IM_REF##*/}
MASK_REF=${MASK_DIR}/${${FILENAME}%_preproc.nii.gz}_brainmask.nii.gz

# Move reference image to folder of matched images
mv ${IM_REF} ${HMATCH_DIR}/.
IM_REF=${HMATCH_DIR}/${FILENAME}

# Histogram matching
echo "Scaling histograms ..."
for IM in $( ls ${IMG_DIR}/*.nii* ); 
do
 FILENAME=${IM##*/}
 MASK=${MASK_DIR}/${${FILENAME}%_preproc.nii.gz}_brainmask.nii.gz
 mrhistmatch scale ${IM} ${IM_REF} ${HMATCH_DIR}/${FILENAME} -mask_input ${MASK} -mask_target ${MASK_REF} -force -quiet
done

# Move results
echo "-------------------------------------------------------------"
mkdir -p ${OUT_DIR}
mv ${HMATCH_DIR}/* ${OUT_DIR}/. -f
echo "Output in ${OUT_DIR} "
echo "-------------------------------------------------------------"
rm -r ${TMP_DIR}
