context("Test assign.kfold")

genin <- read.Genepop("testData/GenepopEx1.txt")

test_that("Perform K-fold cross-validation, genetic data",{
  assign.kfold(genin, k.fold=3, train.loci=1, loci.sample="random", dir="ResKF/", multiprocess=F, skipQ=T)
  expect_true(file.exists("ResKF/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResKF/AnalysisInfo.txt")), 19)
})
unlink("ResKF", re = T)

comin <- compile.data(genin, "testData/varDummy1.csv", skipQ = T)

test_that("Perform K-fold cross-validation, integrated data",{
  assign.kfold(comin, k.fold=3, train.loci=1, loci.sample="random", dir="ResKF/", multiprocess=F, skipQ=T)
  expect_true(file.exists("ResKF/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResKF/AnalysisInfo.txt")), 19)
})
unlink("ResKF", re = T)

nongen <- read.csv("testData/varDummy1.csv", header=T)
pop_label <- c(rep("A",8), rep("B",10), rep("C",6))
nongen <- cbind(nongen, pop_label)

test_that("Perform K-fold cross-validation, non-genetic data",{
  assign.kfold(nongen, k.fold=3, dir="ResKF/", multiprocess=F, skipQ=T)
  expect_true(file.exists("ResKF/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResKF/AnalysisInfo.txt")), 18)
})
unlink("ResKF", re = T)
