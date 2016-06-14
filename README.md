# assignPOP
An R package for genetic or genomic population assignment under supervised machine learning framework

## Description
This package is designed to perform population assignment tests using resampling (Monte-Carlo and K-fold) cross-validation in which samples are divided into training and test sets. A predicting model is built based on the training set and used to predict test individuals' source population. In addition to Monte-Carlo or K-fold resampling for training individuals, genetic loci can be subsampled (either random or based on locus Fst) for training loci. As such, one can evaluate effects of sample size on assignment results. Moreover, the package includes functions that allow user to (1) integrate genetic and non-genetic data for assignment, (2) evaluate which loci are more informative, and (3) generate publication-ready plots. 

## Install assignPOP to R/Rstudio from github
* step 1. Install package devtools by entering `install.packages("devtools")`
* step 2. Call the library, `library(devtools)`
* step 3. Call the function, `install_github("assignPOP", "alexkychen")` 