from fsl_mrs.utils import mrs_io
from fsl_mrs.utils import nifti_mrs_tools as nmrs_tools

# Load raw data
raw_original = mrs_io.read_FID('nifti/raw_original.nii.gz')
raw_synthetic = mrs_io.read_FID('nifti/raw_synthetic.nii.gz')

# Reshape files (24 averages per conditions)
original = nmrs_tools.reshape(raw_original, (32, -1), d5='DIM_DYN', d6='DIM_USER_0')
synthetic = nmrs_tools.reshape(raw_synthetic, (32, -1), d5='DIM_DYN', d6='DIM_USER_0')

# Touch up headers with description
original.set_dim_info('DIM_USER_0', 'Diffusion weighting')
synthetic.set_dim_info('DIM_USER_0', 'Diffusion weighting')

# Add diffusion gradient information
dynamic_headers = {'Bval':[0, 1, 3, 6, 10, 20, 30, 40, 50]}
original.set_dim_header('DIM_USER_0', dynamic_headers)
synthetic.set_dim_header('DIM_USER_0', dynamic_headers)

original.save('nifti/original.nii.gz')
synthetic.save('nifti/synthetic.nii.gz')