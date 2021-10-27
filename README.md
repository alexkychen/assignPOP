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
Changes in ver. 1.2.4 (2021.10.27)
- Update membership.plot - add argument 'plot.k' and 'plot.loci' to skip related question prompt.

<details>
<summary>History</summary>

Changes in ver. 1.2.3 (2021.8.17)
- Update assign.X - (1)Add argument 'common' to specify whether stopping the analysis when inconsistent features between data sets were found. (2)Add argument 'skipQ' to skip data type checking on non-genetic data. (3)Modify argument 'mplot' to handle membership probability plot output.

Changes in ver. 1.2.2 (2020.11.6)
- Update read.Genepop and read.Structure - locus has only one allele across samples will be kept. Use reduce.allele to remove single-allele or low variance loci.
- In ver. 1.2.1, errors might be generated when running assign.MC (and other assignment test functions) due to existence of single-allele loci. (fixed in ver. 1.2.2)

Changes in ver. 1.2.1 (2020.8.24)
- Update read.Genepop to increase file reading speed (~40 times faster)
- Update read.Structure to increase file reading speed (~90 times faster)
- read.Structure now also can handle triploid and tetraploid organisms (see arg. ploidy)
- fix bug in allele.reduce to handle small p threshold across all loci

Changes in ver. 1.2.0 (2020.7.24)
- Add codes to check model name in assign.MC, assign.kfold, assign.X
- Add text to SVM description
- Fix cbind/stringsAsFactors issues in several places for R 4.0
- Able to inject arugments used in models (e.g., gamma in SVM) 

Changes in ver. 1.1.9 (2020.3.16)
- Fix input non-genetic data (x1) error in assign.X

Changes in ver. 1.1.8  (2020.2.28)
- update following functions to work with R 4.0.0
- accuracy.MC, accuracy.kfold, assign.matrix, compile.data, membership.plot
- add stringsAsFactor=T to read.table and read.csv
- temporarily turn off testthat due to its current failure to pass test in Debian system

Changes in ver. 1.1.7  (2019.8.26)
- add broken-stick method for principal component selection in assign.MC, assign.kfold, and assign.X functions
- update accuracy.MC, accuracy.kfold, assign.matrix to handle missing levels of predicted population in test results
- update assign. and accuracy. functions to handle numeric population names

Changes in ver. 1.1.6  (2019.6.8)
- fix multiprocess issue in assign.kfold function

Changes in ver. 1.1.5  (2018.3.23)
- Update assign.MC & assign.kfold to detect pop size and train.inds/k.fold setting
- Update accuracy.MC & assign.matrix to handle test individuals not from every pop
- Slightly modify levels method in accuracy.kfold
- fix bugs in accuracy.plot for K-fold results
- fix membership.plot title positioning and set text size to default

Changes in ver. 1.1.4  (2018.3.8)
- Fix missing assign.matrix function

Changes in ver. 1.1.3  (2017.6.15)
- Add unit tests (using package testthat)

Changes in ver. 1.1.2  (2017.5.13)
- Change function name read.genpop to read.Genepop; Add function read.Structure.
- Update read.genpop function, now can read haploid data
</details>

## Cite this package
Chen, K. Y., Marschall, E. A., Sovic, M. G., Fries, A. C., Gibbs, H. L., & Ludsin, S. A. (2018). assign POP: An R package for population assignment using genetic, non-genetic, or integrated data in a machine-learning framework. *Methods in Ecology and Evolution*. 9(2)439-446. https://doi.org/10.1111/2041-210X.12897

[Papers citing our package](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=14878258167162189944&as_sdt=5)

## Previous version
Previous packages can be found and downloaded at the [releases page](https://github.com/alexkychen/assignPOP/releases)

## Version compatibility (2020.7.24)
assignPOP version 1.1.9 and earlier are not fully compatible with newly released R 4.0.0. 
If you're using R 4.0.0 (or newer), please update your assignPOP to 1.2.0.  

