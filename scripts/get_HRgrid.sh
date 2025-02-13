#!/usr/bin/zsh

# Get HR grid in brain FOV
# Diana Giraldo, August 2022
# It uses commands from mrtrix3 

############################################
# Inputs: 
# Folder with LR brain masks
INMSK_DIR=${1}
# Isotropic voxel size of HR grid:
ISOVOX=${2}

# Outputs:
# HR grid file name
OUT_IM=${3}
############################################

# Make temporary directory
TMP_DIR=$(mktemp -d)

nLR=$( ls ${INMSK_DIR}/* | wc -l )
echo "-------------------------------------------------------------"
echo "Number of input masks: ${nLR}"

MASK_DIR=${TMP_DIR}/mask
mkdir -p ${MASK_DIR} 

# Copy images to temporary directory
cp ${INMSK_DIR}/*.nii* ${MASK_DIR}/.

# First reference grid/fov: Transversal image 
REF_MASK=$(ls ${MASK_DIR}/*TRA*.nii.gz | head -n 1 )
mv ${REF_MASK} ${TMP_DIR}/orig_LRmask_TRA.nii.gz
REF_MASK=${TMP_DIR}/orig_LRmask_TRA.nii.gz

# Regrid and pad (extend FOV just in case) reference mask
mrgrid ${REF_MASK} regrid -vox ${ISOVOX} -interp nearest - -quiet | mrgrid - pad ${TMP_DIR}/HRgrid1_TRA.nii -axis 2 20,20 -quiet

# Combine with other brain masks
for MASK in $(ls ${MASK_DIR}/*(SAG|COR)*.nii.gz)
do
    ORI=$(echo $(basename ${MASK}) | grep -o "SAG\|COR")
    mrgrid ${MASK} regrid -template ${TMP_DIR}/HRgrid1_TRA.nii -interp nearest ${TMP_DIR}/HRgrid1_${ORI}.nii -force -quiet
done

# Get brain FOV mask by dilating brain masks
mrcat ${TMP_DIR}/HRgrid1_TRA.nii ${TMP_DIR}/HRgrid1_*.nii - -quiet | mrmath - max - -axis 3 -quiet | maskfilter - dilate ${TMP_DIR}/HRgrid1.nii -npass 20 -force -quiet

# crop grid to FOV
mrgrid ${TMP_DIR}/HRgrid1.nii crop -mask ${TMP_DIR}/HRgrid1.nii ${TMP_DIR}/HRgrid_fov1.nii -force -quiet

# Pad HR fov grid such that all dimensions can be divided by a factor SF = slice spacing of reference
SF=$( mrinfo ${REF_MASK} -spacing -config RealignTransform 0 | cut -d" " -f3 )
SF=$( printf "%0.0f" $SF )
imdim=$(mrinfo ${TMP_DIR}/HRgrid_fov1.nii -size)
xod=$(echo $imdim | cut -d' ' -f1)
xnd=$(( ((${xod} + ${SF} - 1) / ${SF}) * ${SF} ))
xp1=$(( (${xnd} - ${xod})/2 ))
xp2=$(( (${xnd} - ${xod}) - ${xp1} ))
yod=$(echo $imdim | cut -d' ' -f2)
ynd=$(( ((${yod} + ${SF} - 1) / ${SF}) * ${SF} ))
yp1=$(( (${ynd} - ${yod})/2 ))
yp2=$(( (${ynd} - ${yod}) - ${yp1} ))
zod=$(echo $imdim | cut -d' ' -f3)
znd=$(( ((${zod} + ${SF} - 1) / ${SF}) * ${SF} ))
zp1=$(( (${znd} - ${zod})/2 ))
zp2=$(( (${znd} - ${zod}) - ${zp1} ))
mrgrid ${TMP_DIR}/HRgrid_fov1.nii pad ${TMP_DIR}/HRgrid_fov2.nii.gz -axis 0 ${xp1},${xp2} -axis 1 ${yp1},${yp2} -axis 2 ${zp1},${zp2} -force -quiet

# Copy HR grid
mv ${TMP_DIR}/HRgrid_fov2.nii.gz ${OUT_IM}
echo "Output in ${OUT_IM} "
echo "-------------------------------------------------------------"

rm -r ${TMP_DIR}
