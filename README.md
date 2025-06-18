# Sviatoslav Kharuk, Master's Thesis 2024
This repository contains the code used for my master thesis "Investigating batch effect and cell fate dynamics in the Drosophila melanogaster intestine using single-cell RNA velocity".

[Annotation of the master's thesis](https://drive.google.com/file/d/16_ereruw8RZpe3Ub4xuc2sNt1NO-4Dkj/view?usp=sharing)

Preprint: Proliferation and differentiation of intestinal stem cells depends on the zinc finger transcription factor BCL11/Chronophage - https://doi.org/10.1101/2024.09.08.611891
 

## Abstract
The study aims to investigate the impact of data integration approaches and cell dynamics in the Drosophila intestine using single-cell RNA sequencing and RNA velocity. The study compared two types of sequencing for single-cell RNA velocity, in particular, sequencing of the 5′- and 3′- ends of RNA. The quality of the RNA velocity data was assessed and how it affects the ability of the model to perform the analysis correctly. We also presented an approach to RNA velocity analysis and showed the impact of data integration on the output of the RNA velocity model. In addition, it was shown that the use of certain cell types for the model affects the trajectory inference, and it was shown how knocking out the Notch gene in the fruit fly gut affects cell differentiation.

## Structure
The repository is structured as follows:
* main_analyses.ipynb - contains the main analyses of the study
* functions.py - contains the functions used in the main analyses
* envs - conda environments with dependencies (from cluster and local machine)


Data is available upon reasonable request.