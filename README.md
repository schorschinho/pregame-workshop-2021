# pregame-workshop-2021

This repository contains the material for the pregame of the “Best practices and tools for Diffusion MR Spectroscopy” workshop, taking place at the Lorentz center in Leiden (NL) in September 2021. 

The game consists of processing and fitting a simulated dataset (Synthetic-Data), allowing to compare the results to a ground truth, and two real datasets acquired on a human subject at 7T (InVivo-Data): one with a DW-STEAM, one with a DW-sLASER.

Datasets contain several b-values, and for each b-value, the individual repetitions. Real datasets contain 3 directions per b-value.
We provide the corresponding LCModel basis sets, but feel free to use your simulation tools. An example of a LCModel control file will be provided too. 
A fitted experimental macromolecule baseline is included in the LCModel basis set.  
The parameters of acquisition, the b-values and the files organisation are provided in the readme files of each folder. 
