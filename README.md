# assignPOP
An R package for population assignment using genomic or integrated data in a machine-learning framework

## Description
This package is designed to perform population assignment tests using large genomic data or integrated (genetic plus non-genetic) data. It employs resampling (Monte-Carlo and K-fold) cross-validation approaches coupled with multiple machine-learning classification algorithms to evaluate perfomance of predictive models of population assignment. It allows users to subset varying sample sizes of individuals and loci, either randomly or based on locus Fst, for cross-validation. It uses PCA for data dimension reduction and outputs results in publication-quality graphs.

## Install assignPOP to R/Rstudio from github
In your R/Rstudio console,
* step 1. Install devtools package by entering `install.packages("devtools")`
* step 2. Import the library, `library(devtools)`
* step 3. Call the function, `install_github("assignPOP", "alexkychen")` 

## Package web page
Please visit our website for more infomration
* [http://alexkychen.github.io/assignPOP](http://alexkychen.github.io/assignPOP)

## Note
This package was tested under R-3.2.4