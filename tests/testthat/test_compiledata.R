context("Test compile.data")

test_that("Compile genetic and non-genetic data",{
  expect_output(str(compile.data(read.Genepop("testData/GenepopEx1.txt"), "testData/varDummy1.csv", skipQ = T)),"List of 5")
  expect_output(str(compile.data(read.Structure("testData/StructureEx1.txt"), "testData/varDummy2.csv", skipQ = T)),"List of 5")
})