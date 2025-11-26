Main Figures
Figure 2 — Variation in microbial amino sugar content

Script: analysis/R/06_Figure2_TraitDistributions.R
Inputs:

data/bactMurA.csv

data/fungGlcN.csv
Outputs:

analysis/figures/Fig2_panel.png

manuscript/figures/Fig2.pdf

Figure 3 — Global sensitivity of CF components (Sobol indices)

Script: analysis/R/07_Figure3_SobolIndices.R
Inputs:

processed Sobol index outputs (generated externally or within script)
Outputs:

analysis/figures/Fig3_panel.png

manuscript/figures/Fig3.pdf

Figure 4 — Trait–function relationships for CFB and CFF

Script: analysis/R/07_Figure4_Traits_vs_CF_v2.R
Inputs:

trait distributions and CF calculations (loaded or computed in-script)
Outputs:

analysis/figures/Fig4_panel.png

manuscript/figures/Fig4.pdf

Figure 5 — Diagnostics of the CF pipeline

Script: analysis/R/08_Figure5_Diagnostics.R
Inputs:

modeled necromass terms and diagnostic ratios
Outputs:

analysis/figures/Fig5_panel.png

manuscript/figures/Fig5.pdf

Figure 6 — Model comparisons (CF vs C-mass)

Script: analysis/R/09_Figure6_ModelComparisons_v2.R
Inputs:

model outputs for CFstd, CFHu2024, C_massMurA+GlcN, C_massGlcN
Outputs:

analysis/figures/Fig6_panel.png

manuscript/figures/Fig6.pdf

Figure 7 — Structural feasibility and SOC bounds

Script: analysis/R/10_Figure7_Feasibility.R
Inputs:

joint MS–GS–SOC distributions and CFF-based feasibility surfaces
Outputs:

analysis/figures/Fig7_panel.png

manuscript/figures/Fig7.pdf

Notes for Reviewers

All figure scripts write only their final panels; intermediate objects (e.g., LHS samples, Sobol tables) are created within scripts or loaded from data/ as needed.

Figure numbering in the repository corresponds directly to manuscript figure numbering.

Scripts do not require running upstream “01–05” pipelines; each is self-contained at the figure level.

R environment: see analysis/R/01_setup.R. Python smoke tests are handled in CI but are not required to regenerate figures.