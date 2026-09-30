# P-TIP Stratified Care

Analysis code for the manuscript:

**Stratified Care for Child and Adolescent Conduct Problems: A Data-Driven Decision Rule for Initial Treatment Intensity**

## Overview

This repository contains the R code used for the analyses reported in the manuscript.

The analyses include:

- development and evaluation of prognostic models for lower-intensity treatment (LIT) and higher-intensity treatment (HIT)
- comparison of candidate prediction algorithms
- calculation of the personalized advantage index (PAI)
- derivation and evaluation of the treatment-intensity decision rule
- sensitivity analyses
- generation of tables and figures

## Data availability

The data analyzed in this study originate from multiple clinical trials and a routine-care dataset. Access to the original datasets is governed by the relevant study protocols and data-sharing agreements and can be requested as outlined in the corresponding publications.

The analysis code is provided to support transparency and reproducibility of the reported analyses.

## Software

Analyses were conducted in R.

Required R packages are specified in the individual analysis scripts.

## Repository structure

The repository contains the R analysis scripts used to reproduce the analyses reported in the manuscript.

The scripts are numbered in the order in which they are intended to be run:

1. data splitting and imputation
2. linear regression
3. lasso regression
4. decision tree
5. random forest
6. gradient boosting model
7. model-performance table
8. personalized predictions
9. ROC analysis
10. output tables
11. figures
12. treatment-duration sensitivity analysis
13. concordance and sensitivity analyses

The participant-level datasets are not included in this repository. The scripts therefore assume access to the required input data and use relative paths within the project structure.

## Citation

Please cite the associated manuscript when using this code.

## Contact

Felix Oswald  
University of Freiburg  
Email: felix.oswald@psychologie.uni-freiburg.de

Christopher Hautmann  
University Hospital Cologne  
Email: christopher.hautmann@uk-koeln.de