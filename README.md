[![Travis-CI Build Status](https://travis-ci.org/alexkychen/assignPOP.svg?branch=master)](https://travis-ci.org/alexkychen/assignPOP)
[![CRAN status](http://www.r-pkg.org/badges/version/assignPOP)](https://cran.r-project.org/package=assignPOP)
[![GitHub release](https://img.shields.io/github/release/alexkychen/assignPOP.svg)](https://github.com/alexkychen/assignPOP/releases)

# assignPOP <img src="https://www.r-project.org/logo/Rlogo.svg" width="48">

Population Assignment using Genomic, Non-Genetic or Integrated Data in a Machine-Learning Framework

## Description
This R package helps perform population assignment using a machine learning framework. It employs supervised machine learning methods to evaluate the discriminatory power of your known data set, and is capable of analyzing **large genetic, non-genetic, or integrated (genetic plus non-genetic) data** sets. This framework is also designed for solving the upwardly biased issue that was discussed in previous studies. Other features are listed below.

- Use principle component analysis (PCA) for dimensionality reduction (or data transformation)
- Use Monte-Carlo cross-validation to evaluate the variation of assignment accuracy
- Use *K*-fold cross-validation to estimate membership probability
- Allow to resample training individuals with various proportions or numbers
- Allow to resample training loci with various proportions either randomly or based on locus *Fst* value
- Provide several machine learning classifiers, including LDA, SVM, naive Bayes, decision tree, and random forest to build tunable predictive models.
- Output results in publication-quality plots while being editable using ggplot2 library

## Install assignPOP
- Install from CRAN [![CRAN status](http://www.r-pkg.org/badges/version/assignPOP)](https://cran.r-project.org/package=assignPOP)
  * Simply enter `install.packages("assignPOP")` in your R console

- Install from Github [![GitHub release](https://img.shields.io/github/release/alexkychen/assignPOP.svg)](https://github.com/alexkychen/assignPOP/releases)
  * step 1. Install devtools package by entering `install.packages("devtools")`
  * step 2. Import the library, `library(devtools)`
  * step 3. Then enter `install_github("alexkychen/assignPOP")` 

## Package tutorial
Please visit the following site for more infomration
* [http://alexkychen.github.io/assignPOP/](http://alexkychen.github.io/assignPOP/)