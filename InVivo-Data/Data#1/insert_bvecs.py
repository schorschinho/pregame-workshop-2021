from fsl_mrs.utils import mrs_io
from fsl_mrs.utils import nifti_mrs_tools as nmrs_tools

# Load raw data
raw_metab = mrs_io.read_FID('nifti/raw_metab.nii.gz')
raw_water = mrs_io.read_FID('nifti/raw_water.nii.gz')

# Reshape metab file (24 averages per conditions) and water (4 averages per condition)
metab = nmrs_tools.reshape(raw_metab, (24, -1), d5='DIM_DYN', d6='DIM_USER_0')
water = nmrs_tools.reshape(raw_water, (4, -1), d5='DIM_DYN', d6='DIM_USER_0')

# Touch up headers with description
metab.set_dim_info('DIM_USER_0', 'Variable b values (4) and directions (3)')
water.set_dim_info('DIM_USER_0', 'Variable b values (4) and directions (3)')

# Add diffusion gradient information
dynamic_headers = {'Bval':[0,907,907,907,2155,2155,2155,3956,3956,3956],
                   'direction':[[0,0,0],
                                [1,1,-0.5],
                                [1,-0.5,1],
                                [-0.5,1,1],
                                [1,1,-0.5],
                                [1,-0.5,1],
                                [-0.5,1,1],
                                [1,1,-0.5],
                                [1,-0.5,1],
                                [-0.5,1,1]]}
metab.set_dim_header('DIM_USER_0', dynamic_headers)
water.set_dim_header('DIM_USER_0', dynamic_headers)

metab.save('nifti/metab.nii.gz')
water.save('nifti/water.nii.gz')