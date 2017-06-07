context("Test assign.MC")

genin <- read.Genepop("testData/GenepopEx1.txt")

test_that("Perfrom Monte-Calro cross-validation - loci.sample=random",{
  assign.MC(genin, dir="ResMC/", train.inds=0.5, train.loci=1, loci.sample="random", iterations=3, multiprocess=F)
  expect_true(file.exists("ResMC/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResMC/AnalysisInfo.txt")), 20)
})
unlink("ResMC", re = T)

test_that("Perfrom Monte-Calro cross-validation - loci.sample=Fst",{
  assign.MC(genin, dir="ResMC/", train.inds=0.5, train.loci=1, loci.sample="Fst", iterations=3, multiprocess=F)
  expect_true(file.exists("ResMC/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResMC/AnalysisInfo.txt")), 20)
})
unlink("ResMC", re = T)

test_that("Perfrom Monte-Calro cross-validation - model=lda",{
  assign.MC(genin, dir="ResMC/", model="lda", train.inds=3, train.loci=c(0.5,1), loci.sample="random", iterations=3, multiprocess=F)
  expect_true(file.exists("ResMC/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResMC/AnalysisInfo.txt")), 20)
})
unlink("ResMC", re = T)

test_that("Perfrom Monte-Calro cross-validation - model=naiveBayes",{
  assign.MC(genin, dir="ResMC/", model="naiveBayes", train.inds=c(0.5,0.6), train.loci=c(0.5,1), loci.sample="random", iterations=3, multiprocess=F)
  expect_true(file.exists("ResMC/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResMC/AnalysisInfo.txt")), 20)
})
unlink("ResMC", re = T)

test_that("Perfrom Monte-Calro cross-validation - model=tree",{
  assign.MC(genin, dir="ResMC/", model="tree", train.inds=c(3,4), train.loci=c(0.5,1), loci.sample="random", iterations=3, multiprocess=F)
  expect_true(file.exists("ResMC/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResMC/AnalysisInfo.txt")), 20)
})
unlink("ResMC", re = T)

test_that("Perfrom Monte-Calro cross-validation - model=randomForest",{
  assign.MC(genin, dir="ResMC/", model="randomForest", train.inds=0.5, train.loci=1, loci.sample="random", iterations=3, multiprocess=F)
  expect_true(file.exists("ResMC/AnalysisInfo.txt"))
  expect_equal(length(readLines("ResMC/AnalysisInfo.txt")), 20)
})
unlink("ResMC", re = T)
