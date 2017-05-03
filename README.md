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
- Output results in publication-quality plots that can be edited using ggplot2 functions

## Install assignPOP
You can install the package from CRAN or this Github respository. The lastest version should be on the Github, whereas the CRAN hosts the stable version. See the version number on the badges. 

- To install from CRAN [![CRAN status](http://www.r-pkg.org/badges/version/assignPOP)](https://cran.r-project.org/package=assignPOP)
  * Simply enter `install.packages("assignPOP")` in your R console

- To install from Github [![GitHub release](https://img.shields.io/github/release/alexkychen/assignPOP.svg)](https://github.com/alexkychen/assignPOP/releases)
  * step 1. Install devtools package by entering `install.packages("devtools")`
  * step 2. Import the library, `library(devtools)`
  * step 3. Then enter `install_github("alexkychen/assignPOP")` 

Note: When you install the package from Github, you may need to install additional packages before the assignPOP can be successfully installed. Follow the hints that R provided and then re-run `install_github("alexkychen/assignPOP")`.

## Package tutorial
Please visit our tutorial website for more infomration
* [http://alexkychen.github.io/assignPOP/](http://alexkychen.github.io/assignPOP/)

## What's new

- 2017.5.2 Update read.genpop function, now can read haploid data

## Citation
Chen, K-Y., Marschall, E.A. Sovic M.G., Fries, A.C., Gibbs, H.L., Ludsin, S.A. (2017) assignPOP: Population Assignment using Genetic, Non-Genetic or Integrated Data in a Machine-learning Framework. R package version 1.1.1. https://CRAN.R-project.org/package=assignPOP 

## Previous version
Previous packages can be found and downloaded at [archive branch](https://github.com/alexkychen/assignPOP/tree/archive)
