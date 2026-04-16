# UTPC-paper

This repository contains files to reproduce the figures and analyses described within the associated paper: Jean-François Arnoldi, Andrew L Jackson, Ignacio Peralta-Maraver & Nicholas Payne 2025 [A universal thermal performance curve arises in biology and ecology](https://doi.org/10.1073/pnas.2513099122). 

The version of the repository at the time of publication is available at [https://github.com/AndrewLJackson/UTPC-paper/tree/PNAS-2513099122]. Subsequent additions to the repository reflect on-going work and additional notes to aid understanding.

The main file is `Figures-UTPC-paper.Rmd` which reproduces all the figures in the paper. This file calls additional helper functions contained in `Analyses_thermal_functions.R`. To do so it loads three raw data files and one data file of parameter estimates from a cited paper Rezende & Bozinovic 2019: these files are contained in the folder `data/`

* Biochemistry_raw_data.xlsx
* Figness_raw_data.xlsx
* Running_raw_data.xlsx
* All_parameters_data.csv

The file `import-fit-merge.Rmd` loads these datasets, merges them together and adds additional model parameter estimates as fitted by our UTPC approach. The output from this process is a single tibble (data.frame) object saved in `output/combined_dat_params.rda` which is subsequently loaded by `Figures-UTPC-paper.Rmd`. 

A more detailed mathematical explanation of how Arroyo et al's model can be brought into line with our UTPC model by introducing a high temperature cutoff is provided in the file 'modifying-arroyo-et-al-model-with-cutoff.pdf'.
