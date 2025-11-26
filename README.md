# Trait Variation and Structural Bias in Conversion Factor Models of Soil Microbial Necromass Carbon

_Aubrey K. Fine, Fernanda Santos, Larry M. York_  
Environmental Sciences Division & Biosciences Division, Oak Ridge National Laboratory

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![renv](https://img.shields.io/badge/reproducibility-renv-blue)](https://rstudio.github.io/renv/)

## 🔬 Overview

This repository contains code, data, and figures for the manuscript:

> **Trait Variation and Structural Bias in Conversion-Factor Models of Soil Microbial Necromass Carbon**  
> Submitted to *Soil Biology & Biochemistry* (202X)

We evaluate the assumptions and uncertainty surrounding fixed conversion factors (CFs) used to estimate microbial necromass carbon (necC) from amino sugar biomarkers, muramic acid (MurA) and glucosamine (GlcN). The analysis combines:

- taxonomic variation in MurA and GlcN concentrations,  
- CF-based necromass models versus a direct carbon–mass approach, and  
- global sensitivity analysis (Sobol) to identify key drivers of uncertainty in the necromass fraction of SOC (fₙₑcC).

---

## 📁 Repository Structure

- `analysis/`: R and Python scripts and figure generation.
- `data/`: raw, processes, and derived datasets (+ metadata).
- `manuscript/`: manuscript text, journal figures, submission files.
- `supplement/`: supplementary figures and tables, plus a methods README mapping scripts -> output.
- `env/`: renv lock + conda environment. 
- `github/`: CI workflows.
---

## 📦 Software Requirements

### R

Analyses are written for **R ≥ 4.3**. Core packages include (but are not limited to):

- `tidyverse`
- `mgcv`
- `ggridges`
- `gratia`
- `patchwork`
- `forcats`
- `here`
- `lhs`, `boot`, `rstatix` (for uncertainty and statistics)

To install them:
```r
source("install_packages.R")
```

We recommend using `renv` for reproducible environments:
```r
install.packages("renv")  # if needed
renv::restore()

---

### Python

Some functionality (e.g., sensitivity analysis scaffolding) assumes a Python environment described in: 
- env/environment.yml

Create it with: 
```
conda env create -f env/environment.yml
conda activate necromass-cf



## 🧪 Reproducibility

To reproduce the key findings and figures:

1. Clone this repo and set your working directory
    ```r
    setwd("path/to/necromass-cf-uncertainty")
    ```
2. Run R setup
This loads packages, sets theme defaults, and checks environment:

Rscript analysis/R/01_setup.R

3. Generate figures (main paper)
Each script is self-contained and writes:
-preview panels → analysis/figures/

-publication-ready panels → manuscript/figures/

# Figure 2
Rscript analysis/R/06_Figure2_TraitDistributions.R

# Figure 3
Rscript analysis/R/07_Figure3_SobolIndices.R

# Figure 4
Rscript analysis/R/07_Figure4_Traits_vs_CF_v2.R

# Figure 5
Rscript analysis/R/08_Figure5_Diagnostics.R

# Figure 6
Rscript analysis/R/09_Figure6_ModelComparisons_v2.R

# Figure 7
Rscript analysis/R/10_Figure7_Feasibility.R

---

## 📊 Figures and Tables

- **Main figures:** 1-7
- **Supplementary figures:** S1-S11
- **Tables:** 1-2 and S1-S9

---

## 📜 License

This work is released under the MIT License.

---

## 📫 Contact

For questions, contact [Dr. Aubrey Fine](mailto:fineaubrey@gmail.com)
---

## 🔗 Citation

Please cite the following if you use this repository:

> Fine, A.K., Santos, F., York, L. (202X). *Trait Variation and Structural Bias in Conversion Factor Models of Soil Microbial Necromass Carbon*. Submitted to *Soil Biology and Biochemistry*.  
> DOI: https://doi.org/10.5281/zenodo.16905679



