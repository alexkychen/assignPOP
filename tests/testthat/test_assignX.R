context("Test assign.X")

genin <- read.Genepop("testData/GenepopEx1.txt")
genin_unknown <- read.Genepop("testData/GenepopUnk.txt")

test_that("Perform one-time assignment test on unknown individuals",{
  assign.X(x1=genin, x2=genin_unknown, dir="AssignRes/", mplot=T)
  expect_true(file.exists("AssignRes/AnalysisInfo.txt"))
  expect_true(file.exists("AssignRes/AssignmentResult.txt"))
})
unlink("AssignRes", re = T)