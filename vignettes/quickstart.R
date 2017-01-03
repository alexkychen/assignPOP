## ---- eval=FALSE---------------------------------------------------------
#  install.packages("assignPOP")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("alexkychen/assignPOP")

## ---- warning=FALSE------------------------------------------------------
library(assignPOP)

## ---- eval=FALSE---------------------------------------------------------
#  varName1 <- read.genpop("YourGenepopFileName.txt", pop.names=c("pop_A","pop_B","pop_C"))
#  

## ---- eval=FALSE---------------------------------------------------------
#  varName2 <- compile.data(varName1, "YourNonGeneticFile.csv")

## ---- eval=FALSE---------------------------------------------------------
#  assign.MC(varName1, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1), loci.sample="fst", iterations=30, model="lda", dir="ResultFolder/" )

## ---- eval=FALSE---------------------------------------------------------
#  assign.kfold(varName2, k.fold=c(3, 4, 5), train.loci=c(0.1, 0.25, 0.5, 1), loci.sample="random", model="svm", dir="ResultFolder2/" )

## ---- eval=FALSE---------------------------------------------------------
#  #Import and concatenate known data
#  genetic_known <- read.genpop("YourGenepop_Kn.txt", pop.names=c("A","B","C"))
#  integrate_known <- compile.data(genetic_known, "NonGeneticData_Kn.csv")
#  
#  #Import and concatenate unknown data
#  genetic_unknown <- read.genpop("YourGenepop_Un.txt")
#  integrate_uknown <- compile.data(genetic_unknown, "NonGeneticData_Un.csv")
#  

## ---- eval=FALSE---------------------------------------------------------
#  assign.X(x1=integrate_known, x2=integrate_uknown, model="naiveBayes", dir="ResultFolder3/")
#  

## ---- eval=FALSE---------------------------------------------------------
#  #For Monte-Carlo results
#  accuRes_MC <- accuracy.MC(dir="ResultFolder/")
#  
#  #For K-fold results
#  accuRes_KF <- accuracy.kfold(dir="ResultFolder2/")

## ---- echo=FALSE---------------------------------------------------------
accuRes_MC <- read.table("../inst/extdata/Rate.txt", header = T)

## ----  fig.align='center', fig.width=7, fig.asp=0.6----------------------
#Import ggplot2 library to fine tune the plot
library(ggplot2)

#Make an assignment accuracy box plot and set Y-axis range between 0 and 1
accuracy.plot(accuRes_MC, pop=c("all", "pop_A", "pop_B", "pop_C")) + ylim(0,1)

## ---- eval=FALSE---------------------------------------------------------
#  membership.plot(dir="ResultFolder2/", style = 1)
#  

## ---- echo=FALSE, out.width="700px"--------------------------------------
knitr::include_graphics("Membership-genetic-morph-style1.png")

## ---- eval=FALSE---------------------------------------------------------
#  # Check out top 20 high Fst loci
#  check.loci(dir="ResultFolder/", top.loci = 20)
#  

