# pregame-workshop-2021

This repository contains the material for the pregame of the “Best practices and tools for Diffusion MR Spectroscopy” workshop, taking place at the Lorentz center in Leiden (NL) in September 2021. 

The game consists of processing and fitting a simulated dataset (Synthetic-Data), allowing to compare the results to a ground truth, and two real datasets acquired on a human subject at 7T (InVivo-Data): one with a DW-STEAM, one with a DW-sLASER.

Datasets contain several b-values, and for each b-value, the individual repetitions. Real datasets contain 3 directions per b-value.
We provide the corresponding LCModel basis sets, but feel free to use your simulation tools. An example of a LCModel control file will be provided too. 
A fitted experimental macromolecule baseline is included in the LCModel basis set.  
The parameters of acquisition, the b-values and the files organisation are provided in the readme files of each folder. 


Deliverables:

 
For the Synthetic-Data:
* a .mat or .csv file  of the residuals that you obtain when computing the difference of your processed spectra at each b-value versus the ground truth at each b-value (before quantification). Please specify if you applied anything else than phase & frequency correction (or if you did not apply it). 
* a visualisation of the residuals after spectral quantification (i.e. data-fit) for each b-value. Please specify if you simulated your own basis set or if you used the one provided.
* a .mat or .csv file containing the signal decay for the provided b-values for metabolites that are (in your opinion) trustworthy. 
* if you are using an existing pipeline, and that you make substantial changes to optimise it, please send the initial and the last results, and write a short note with the changes you made that brought an improvement. 
any observation you’d like to make.
 
For the InVivo-Data, each b-value contains 3 directions, you can average them or treat them independently, just let us know what you decided. Please send us:
* a .mat or .csv of your processed spectra. Please specify if you applied anything else than phase & frequency correction (or if you did not apply it). 
* a visualisation of the residuals after spectral quantification (i.e. data-fit). Please specify if you simulated your own basis set or if you used the one provided.
* a .mat or .csv file containing the signal decay for the provided b-values for metabolites that are (in your opinion) trustworthy.
* if you are using an existing pipeline, and that you make substantial changes to optimise it, please send the initial and the last results, and write a short note with the changes you made that brought an improvement. 
* any observation you’d like to make.
