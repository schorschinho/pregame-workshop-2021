#!/bin/sh
# Do the initial conversion of the SPAR and SDAT files using spec2nii
spec2nii philips -t DIM_USER_0 -f raw_synthetic -o nifti SyntheticData\^X\^\^\^_DWS_b9_nsa32_raw_act.SDAT SyntheticData\^X\^\^\^_DWS_b9_nsa32_raw_act.SPAR
spec2nii philips -t DIM_USER_0 -f raw_original -o nifti OriginalData\^X\^\^\^_DWS_b9_nsa32_raw_act.SDAT OriginalData\^X\^\^\^_DWS_b9_nsa32_raw_act.SPAR

# Run the dedicated python script to reshape and insert bvecs
# Requires FSL-MRS installation
python insert_bvecs.py

# Display results
mrs_tools info nifti/synthetic.nii.gz