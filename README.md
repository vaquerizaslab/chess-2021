# Analysis code for Ing-Simmons, Machnik and Vaquerizas - response to Lee et al.

This repository contains code for analysing Hi-C data from Diaz et al. 2018 (`workflow/Snakefile`) and data from Rao et al. 2014 and
Wutz et al. 2017. obtained from the 4D Nucleome website (`workflow/4dn.smk`). These are set up as separate Snakemake workflows for ease of data processing.
The file `scripts/figures.Rmd` contains R code to generate most figures. Hi-C visualisations were generated using FAN-C and the python scripts in `scripts/`.

The .html output file with figures can be viewed [here](https://htmlpreview.github.io/?https://github.com/vaquerizaslab/chess-2021/blob/main/scripts/figures.html).
