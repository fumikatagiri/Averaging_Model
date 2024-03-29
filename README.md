# Averaging_Model

This repository is for supplemental datasets in "An averaging model for analysis and interpretation of high-order genetic interactions".

Supplemental Dataset 1. A folder containing four bacterial count data files (tab-delimited text) for “AvrRpt2_ETI”, “AvrRpm1_ETI”, “flg22_PTI”, and “elf18_PTI”. Each has columns of genotype, treatment, replicate, flat, pot, and colony. The colony column has log10-transformed colony counts (colony forming unit/cm2). Although the summarized data were originally reported in Tsuda et al. 2009: https://doi.org/10.1371/journal.pgen.1000772, these raw data were not published.

Supplemental Dataset 2. A .RData file (R workspace file) containing a list object, “ave.model.mats”, which contains the matrices for the averaging model for 2- to 7-gene systems (equivalents of matrix in Fig. TS1.1 in Text S1 for different number-gene systems).

Supplemental Dataset 3. An R script file (.r file), which was used to generate Fig. 5 from the data in Supplemental Dataset 1 (the data should be in the subfolder "Supplemental.Dataset1" as in this folder). In addition, an empty "outputs" subfolder is required for the outputs of this script (the "outputs" subfolder in this folder is already populated with the actual outputs of the script). The outputs of the script includes "reanalyzed.w.average.model....pdf", which is Fig. 5. The script includes algorithms for selecting significant genes for the analysis and the averaging model. In the script, the object of a 3-gene or 4-gene matrix included in Supplemental Dataset 2 is called “rec.mx”, which is generated by a function, “make.rec.mx”.


[updates made on Nov 12, 2023 for a revised version of the manuscript]
R scripts, "Fig1.ave.model.r", "Fig2.ave.model.r", and "Fig3.ave.model.r", and the "outputs.fig3" subfolder were added. 
- "Fig1.ave.model.r" generates the coefficient values in the tables in Fig. 1.
- "Fig2.ave.model.r" generates Fig. 2 in "outputs" subfolder.
- "Fig3.ave.model.r" uses the “AvrRpt2_ETI” data set and generates Fig. 3 in "outputs" subfolder and multiple intermediate and diagnostic plots in "outputs.fig3" subfolder. "outputs.fig3" subfolder is required to run this script.
