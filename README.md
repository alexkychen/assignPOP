[![Travis-CI Build Status](https://travis-ci.org/alexkychen/assignPOP.svg?branch=master)](https://travis-ci.org/alexkychen/assignPOP)
[![CRAN status](http://www.r-pkg.org/badges/version/assignPOP)](https://cran.r-project.org/package=assignPOP)
[![GitHub release](https://img.shields.io/github/release/alexkychen/assignPOP.svg)](https://github.com/alexkychen/assignPOP/releases)
[![license](https://img.shields.io/github/license/alexkychen/assignPOP.svg)](https://github.com/alexkychen/assignPOP/blob/master/LICENSE.md)

# assignPOP <img src="https://www.r-project.org/logo/Rlogo.svg" width="48">

Population Assignment using Genetic, Non-Genetic or Integrated Data in a Machine-learning Framework

## Description
This R package helps perform population assignment and infer population structure using a machine-learning framework. It employs supervised machine-learning methods to evaluate the discriminatory power of your data collected from source populations, and is able to analyze **large genetic, non-genetic, or integrated (genetic plus non-genetic) data** sets. This framework is designed for solving the upward bias issue discussed in previous studies. Main features are listed as follows.

- Use principle component analysis (PCA) for dimensionality reduction (or data transformation)
- Use Monte-Carlo cross-validation to estimate mean and variance of assignment accuracy
- Use *K*-fold cross-validation to estimate membership probability
- Allow to resample various sizes of training datasets (proportions or fixed numbers of individuals and proportions of loci)
- Allow to choose from various proportions of training loci either randomly or based on locus *Fst* values
- Provide several machine-learning classification algorithms, including LDA, SVM, naive Bayes, decision tree, and random forest, to build tunable predictive models.
- Output results in publication-quality plots that can be modified using ggplot2 functions

## Install assignPOP
You can install the released version from CRAN or the up-to-date version from this Github respository.

- To install from CRAN
  * Simply enter `install.packages("assignPOP")` in your R console

- To install from Github
  * step 1. Install devtools package by entering `install.packages("devtools")`
  * step 2. Import the library, `library(devtools)`
  * step 3. Then enter `install_github("alexkychen/assignPOP")` 

Note: When you install the package from Github, you may need to install additional packages before the assignPOP can be successfully installed. Follow the hints that R provided and then re-run `install_github("alexkychen/assignPOP")`.

## Package tutorial
Please visit our tutorial website for more infomration
* [http://alexkychen.github.io/assignPOP/](http://alexkychen.github.io/assignPOP/)

## What's new
Changes in ver. 1.1.5  (2018.3.23)
- Slightly modify levels method in accuracy.kfold
- Update assign.MC & assign.kfold to detect pop size and train.inds/k.fold setting
- Update accuracy.MC & assign.matrix to handle test individuals not from every pop
- Slightly modify levels method in accuracy.kfold
- fix bugs in accuracy.plot for K-fold results
- fix membership.plot title positioning and set text size to default

<details>
<summary>History</summary>

Changes in ver. 1.1.4  (2018.3.8)
- Fix missing assign.matrix function

Changes in ver. 1.1.3  (2017.6.15)
- Add unit tests (using package testthat)

Changes in ver. 1.1.2  (2017.5.13)
- Change function name read.genpop to read.Genepop; Add function read.Structure.
- Update read.genpop function, now can read haploid data
</details>

## Cite this package
Chen K-Y, Marschall EA, Sovic MG, Fries AC, Gibbs HL, Ludsin SA. assignPOP: An R package for population assignment using genetic, non-genetic, or integrated data in a machine-learning framework. *Methods in Ecology and Evolution*. 2018;9:439–446. https://doi.org/10.1111/2041-210X.12897

## Previous version
Previous packages can be found and downloaded at [archive branch](https://github.com/alexkychen/assignPOP/tree/archive)
