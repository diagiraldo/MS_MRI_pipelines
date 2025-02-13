# MRI processing pipelines for T2-W FLAIR of people with MS  

In this repository we share the scripts used to process T2-W FLAIR MRIs in the work "Leveraging Hand-Crafted Radiomics on Multicenter FLAIR MRI for Predicting Disability Progression in People with Multiple Sclerosis".

[Read the pre-print](https://www.medrxiv.org/content/10.1101/2025.01.23.25320971v1) 

## Pipeline for Low-Resolution (LR) MRI

Example script [`example_LR_ABD.sh`](example_LR_ABD.sh)

<img src="figures/pipeline_LR.png?raw=True" width="800px" style="margin:0px 0px"/>

## Pipeline for High-Resolution (HR) MRI

Example script [`example_HR_CsT1.sh`](example_HR_CsT1.sh)

<img src="figures/pipeline_HR.png?raw=True" width="800px" style="margin:0px 0px"/>

## Requirements

Pre-processing script depends on:
- [ANTs](https://github.com/ANTsX/ANTs)
- [HD-BET](https://github.com/MIC-DKFZ/HD-BET)

Segmentation is done with:
- [LST](https://www.applied-statistics.de/lst.html) for [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- [SAMSEG](https://surfer.nmr.mgh.harvard.edu/fswiki/Samseg) in [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)

Super-resolution (SR) and combination of SR outputs depend on:
- [PRETTIER-MRI](https://github.com/diagiraldo/PRETTIER)
- [MRtrix3](https://www.mrtrix.org/)

## Funding

This project received funding from the Flemish Government under the [“Flanders AI Research Program"](https://www.flandersairesearch.be/en).

## Contact

Diana L. Giraldo Franco [@diagiraldo](https://github.com/diagiraldo)

