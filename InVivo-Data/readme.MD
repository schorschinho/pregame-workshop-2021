# In Vivo Datasets

## Dataset 1
Acquired in the human brain (parietal WM, VOI 30x15x20) with a DW-STEAM sequence (TE=45ms, TM=60ms) at 7 tesla (Philips scanner)

The DW-STEAM sequence is set up as follows:
- The diffusion gradients had a duration of 15 ms
- The duration between the start of the first and second diffusion gradients (diffusion time) was 79 ms.
- See [provided schematic](./Data%231/DW-STEAM.PNG) for sequence details.

The metabolite dataset was acquired with:
- 4 diffusion gradient amplitudes (G=0, 15, 25 and 35mT) and along 3 orthogonal directions ([1,1,-0.5],[-0.5,1,1],[1,-0.5,1])
- a total of 24 averages per condition (for each b-value and direction)
- outer-volume suppression
- VAPOR scheme (with squeezed delays) for water suppression
- PPU for cardiac triggering
- other acquisition parameters:
    - npoints = 1024
    - spectral width=3000 Hz

The water dataset was acquired with:
- exactly same parameters but without water suppresion and with a lower number of averages 
- a total of 4 averages per condition

### What we provide:
- A file with data for metabolites after coil-combination, e.g. a total of 24x(1+3x3)=240 spectra, each with 1024 points
- A file with data for water after coil-combination, e.g. a total of 4x(1+3x3)=40 spectra, each with 1024 points
- A file with all b-values (which account for the cross terms with imaging gradients)
- A file with sequence chronograms (for each gradient amplitude and direction)
- A basis-set for analysis of the spectra with LCModel or other preferred software
- A control file for analysis with LCModel (this gives help to get started, but feel free to adapt to your preferred parameters)
- Data converted to the proposed NIfTI-MRS standard. Data is stored as a 1x1x1x1024x24x10 NIfTI file. The b-values and vecotrs are encoded within a single diemension due to the unequal number of directions.
- A 3dT1 image (after fsl_deface). NB: because of the VOI location, we used a dielectric pad.


## Dataset 2
Acquired in the human brain (posterior cingulate cortex, VOI 20x20x20) with a DW-sLASER sequence (TE=125ms) at 7 tesla (Philips scanner).

The DW-sLASER sequence is set as follows:
- the diffusion gradients had a duration of 8 ms,
- the duration between the start of the first and second diffusion blocks (diffusion time) was 60 ms.
- See [provided schematic](./Data%232/DW-sLASER.PNG) for sequence details.

The metabolite dataset was acquired with:
- 3 diffusion gradient amplitudes (G=10, 18 and 26mT) and along 3 orthogonal directions ([1,1,-0.5],[-0.5,1,1],[1,-0.5,1])
- a total of 24 averages per condition (for each b-value and direction)
- VAPOR scheme (with squeezed delays) for water suppression
- PPU for cardiac triggering
- other acquisition parameters:
    - npoints = 1024
    - spectral width=3000 Hz

The water dataset was acquired with:
- exactly same parameters but without water suppresion and with a lower number of averages
- a total of 4 averages per condition

### We provide:
- A file with data for metabolites after coil-combination, e.g. a total of 24x3x3=216 spectra, each with 1024 points
- A file with data for water after coil-combination, e.g. a total of 4x3x3=36 spectra, each with 1024 points
- A file with all b-values (which account for the cross terms with imaging gradients)
- A file with sequence chronograms (for each gradient amplitude and direction)
- A basis-set for analysis of the spectra with LCModel or other preferred software
- A control file for analysis with LCModel (this gives help to get started, but feel free to adapt to your preferred parameters)
- Data converted to the proposed NIfTI-MRS standard. Data is stored as a 1x1x1x1024x24x3x3 NIfTI file. The final three diemensions encode the averages, b-values and gradient directions respectively.
- A 3dT1 image (after fsl_deface).
