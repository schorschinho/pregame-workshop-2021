from fsl_mrs.utils import mrs_io
from fsl_mrs.utils import nifti_mrs_tools as nmrs_tools

# Load raw data
raw_metab = mrs_io.read_FID('nifti/raw_metab.nii.gz')
raw_water = mrs_io.read_FID('nifti/raw_water.nii.gz')

# Reshape metab file (24 averages per conditions) and water (4 averages per condition)
metab = nmrs_tools.reshape(raw_metab, (24, 3, 3), d5='DIM_DYN', d6='DIM_USER_0', d7='DIM_USER_1')
water = nmrs_tools.reshape(raw_water, (4, 3, 3), d5='DIM_DYN', d6='DIM_USER_0', d7='DIM_USER_1')

# Touch up headers with description
metab.set_dim_info('DIM_USER_0', 'B values s/mm2')
metab.set_dim_info('DIM_USER_1', 'Diffusion gradient directions')

water.set_dim_info('DIM_USER_0', 'B values s/mm2')
water.set_dim_info('DIM_USER_1', 'Diffusion gradient directions')


# Add diffusion gradient information
dynamic_headers_0 = {'Bval':[614, 2113, 4524]}
dynamic_headers_1 = {'direction':[[1,1,-0.5], [1,-0.5,1], [-0.5,1,1]]}

metab.set_dim_header('DIM_USER_0', dynamic_headers_0)
water.set_dim_header('DIM_USER_0', dynamic_headers_0)
metab.set_dim_header('DIM_USER_1', dynamic_headers_1)
water.set_dim_header('DIM_USER_1', dynamic_headers_1)

metab.save('nifti/metab.nii.gz')
water.save('nifti/water.nii.gz')