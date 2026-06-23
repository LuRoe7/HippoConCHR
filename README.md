# Intra-hippocampal connectivity in clinical high-risk for psychosis (CHR-P)
This repository contains code associated with a project investigating intra-hippocampal connectivity in individuals at clinical high risk for psychosis (CHR-P). A preprint describing this work is available at https://doi.org/10.64898/2026.03.10.26348090.  
The project uses demographic, clinical, and functional neuroimaging data from the North American Prodromal Longitudinal Study 3 (NAPLS-3). Data access is governed by the NIMH Data Archive and requires approval through an application led by a principal investigator affiliated with an approved institution. Due to these data sharing restrictions, the underlying data cannot be made publicly available. This repository therefore contains analysis code only.

## Software requirements
- R version 4.5.1
- Mplus version 8.6  
R package requirements are specified within individual scripts. Installation instructions for R and Mplus are available on their respective official websites

## Pipeline overview
Functional MRI data were preprocessed using fMRIPrep. Subsequent smoothing, denoising, and functional connectivity computations were performed in Python using Nilearn. The statistical analysis pipeline integrates R and Mplus-based modelling approaches applied to derived connectivity and clinical variables.

## Analysis workflow
The code is organised as a sequential pipeline of scripts: Scripts are numbered 01_* to 19_*. Scripts must be executed in order. Input and output paths are hard-coded within scripts. The full pipeline runs in approximately 10 minutes.

## Data and key variables
The analysis uses demographic, clinical, and neuroimaging data. Key variables include:
- Intra-hippocampal connectivity
- Negative symptoms
- Attenuated positive symptoms
- Depressive symptoms
- Psychosocial functioning
- Symbol coding performance
- Verbal memory performance

## Outputs
The scripts generate all figures and tables reporting the main results presented in the preprint.

## Citation
Preprint: https://doi.org/10.64898/2026.03.10.26348090

## Contact
Dr Lukas Roell  
Walter Benjamin Research Fellow (German Research Foundation)  
Department of Psychiatry, The University of Melbourne  
E-mail: lukas.roell@unimelb.edu.au  
