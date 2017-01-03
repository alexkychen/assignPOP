[![Travis-CI Build Status](https://travis-ci.org/alexkychen/assignPOP.svg?branch=master)](https://travis-ci.org/alexkychen/assignPOP)
[![CRAN status](http://www.r-pkg.org/badges/version/assignPOP)](https://cran.r-project.org/package=assignPOP)
[![GitHub release](https://img.shields.io/github/release/alexkychen/assignPOP.svg)](https://github.com/alexkychen/assignPOP/releases)

# assignPOP <img src="https://www.r-project.org/logo/Rlogo.svg" width="48">

Population Assignment using Genomic, Non-Genetic or Integrated Data in a Machine Learning Framework

## Description
This R package helps perform population assignment using a machine learning framework. It employs supervised machine learning methods to evaluate the discriminatory power of your known data set, and is capable of analyzing **large genetic, non-genetic, or integrated (genetic plus non-genetic) data** sets. This framework is also designed for solving the upwardly biased issue that was discussed in previous studies. Other features are listed below.

- Use principle component analysis (PCA) for dimensionality reduction (or data transformation)
- Use Monte-Carlo cross-validation to evaluate the variation of assignment accuracy
- Use *K*-fold cross-validation to estimate membership probability
- Allow to resample training individuals with various proportions or numbers
- Allow to resample training loci with various proportions either randomly or based on locus *Fst* value within training data
- Provide several machine learning classifiers, including LDA, SVM, naive Bayes, decision tree, and random forest to build tunable predictive models.
- Output results in publication-quality plots while being editable using ggplot2 library

## Install assignPOP
You can install the package from CRAN or Github respository here. The lastest version should be on the Github, whereas the CRAN hosts the stable version. Check the version on the badges. 

- To install from CRAN [![CRAN status](http://www.r-pkg.org/badges/version/assignPOP)](https://cran.r-project.org/package=assignPOP)
  * Simply enter `install.packages("assignPOP")` in your R console

- To install from Github [![GitHub release](https://img.shields.io/github/release/alexkychen/assignPOP.svg)](https://github.com/alexkychen/assignPOP/releases)
  * step 1. Install devtools package by entering `install.packages("devtools")`
  * step 2. Import the library, `library(devtools)`
  * step 3. Then enter `install_github("alexkychen/assignPOP")` 

Note: When you install the package from Github, you may need to install additional packages before the assignPOP can be successfully installed. Follow the hints that R provided and then re-run `install_github("alexkychen/assignPOP")`.

## Package tutorial
Please visit the following site for more infomration
* [http://alexkychen.github.io/assignPOP/](http://alexkychen.github.io/assignPOP/)

## Citation
Chen, K-Y., Marschall, E.A. Sovic M.G., Fries, A.C., Gibbs, H.L., Ludsin, S.A. (2017) assignPOP: Population Assignment using Genomic, Non-Genetic or Integrated Data in a Machine Learning Framework. R package version 1.1.2. https://CRAN.R-project.org/package=assignPOP 

## Previous version
Previous packages can be found and downloaded at [archive branch](https://github.com/alexkychen/assignPOP/tree/archive)
