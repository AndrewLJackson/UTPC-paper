# UTPC-paper

This repository contains files to reproduce the figures and analyses described within the associated paper _A universal thermal performance curve arises in biology and ecology_ by Jean-Fran&ccedil;ois Arnoldi, Andrew L Jackson, Ignacio Peralta-Maraver & Nicholas Payne.

The main file is `Figures-UTPC-paper.Rmd` which reproduces all the figures in the paper. This file calls additional helper functions contained in `Analyses_thermal_functions.R`. To do so it loads three raw data files and one data file of parameter estimates from a cited paper Rezende & Bozinovic 2019: these files are contained in the folder `data/`

* Biochemistry_raw_data.xlsx
* Figness_raw_data.xlsx
* Running_raw_data.xlsx
* All_parameters_data.csv

The file `import-fit-merge.Rmd` loads these datasets, merges them together and adds additional model parameter estimates as fitted by our UTPC approach. The output from this process is a single tibble (data.frame) object saved in `output/combined_dat_params.rda` which is subsequently loaded by `Figures-UTPC-paper.Rmd`. 