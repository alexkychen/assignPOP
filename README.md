# assignPOP <img src="https://www.r-project.org/logo/Rlogo.svg" width="48">
Population Assignment using Genomic, Non-Genetic or Integrated Data in a Machine-Learning Framework

### Version 1.1.1

## Description
This R package helps perform population assignment using a machine learning framework. It employs supervised machine learning methods to evaluate the discriminatory power of your known data set, and is capable of analyzing **large genetic, non-genetic, or integrated (genetic plus non-genetic) data** sets. This framework is also designed for solving the upwardly biased issue that was discussed in previous studies. Other features are listed below.

- Use principle component analysis (PCA) for dimensionality reduction (or data transformation)
- Use Monte-Carlo cross-validation to evaluate the variation of assignment accuracy
- Use *K*-fold cross-validation to estimate membership probability
- Allow to resample training individuals with various proportions or numbers
- Allow to resample training loci with various proportions either randomly or based on locus *F~ST~* value
- Provide several machine learning classifiers, including LDA, SVM, naive Bayes, decision tree, and random forest to build tunable predictive models.
- Output results in publication-quality plots while being editable using ggplot2 library

## Install assignPOP to R/Rstudio from github
In your R/Rstudio console,
* step 1. Install devtools package by entering `install.packages("devtools")`
* step 2. Import the library, `library(devtools)`
* step 3. Call the function, `install_github("alexkychen/assignPOP")` 

## Package tutorial
Please visit the following page for more infomration
* [http://rpubs.com/alexkychen/assignPOP](http://rpubs.com/alexkychen/assignPOP)
