context("Test reduce.allele")

test_that("reduce.allele() removes loci with minor alleles",{
  #data <- read.Genepop("GenepopEx1.txt")
  expect_output(str(reduce.allele(read.Genepop("testData/GenepopEx1.txt"))),"List of 3")
  #data <- read.Structure("StructureEx1.txt")
  expect_output(str(reduce.allele(read.Structure("testData/StructureEx1.txt"), p = 0.8)),"List of 3")
})