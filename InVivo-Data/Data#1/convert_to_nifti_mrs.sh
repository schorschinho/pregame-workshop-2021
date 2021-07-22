#!/bin/sh
# Do the initial conversion of the SPAR and SDAT files using spec2nii
spec2nii philips -t DIM_USER_0 -f raw_metab -o nifti InVivoData1\^X\^\^\^_DWS_invivo_3b_+b0_Metab_raw_act_tot.SDAT InVivoData1\^X\^\^\^_DWS_invivo_3b_+b0_Metab_raw_act_tot.SPAR
spec2nii philips -t DIM_USER_0 -f raw_water -o nifti InVivoData1\^X\^\^\^_DWS_invivo_3b_+b0_Water_raw_act_tot.SDAT InVivoData1\^X\^\^\^_DWS_invivo_3b_+b0_Water_raw_act_tot.SPAR

# Run the dedicated python script to reshape and insert bvecs
# Requires FSL-MRS installation
python insert_bvecs.py

# Clean up
rm nifti/raw_metab.nii.gz
rm nifti/raw_water.nii.gz

# Display results
mrs_tools info nifti/metab.nii.gz
