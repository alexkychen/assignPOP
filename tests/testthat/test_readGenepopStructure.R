context("Test read.Genepop and read.Structure")

test_that("read.Genepop imports Genepop file as a list", {
  expect_output(str(read.Genepop("GenepopEx1.txt")), "List of 3")
  expect_output(str(read.Genepop("GenepopEx1.txt", pop.names = c("A","B","C"))), "List of 3")
  expect_output(str(read.Genepop("GenepopEx1.txt", pop.names = c("A","B","C"), haploid = T)), "List of 3")
  expect_output(str(read.Genepop("GenepopEx2.txt")), "List of 3")
})

test_that("read.Structure imports Structure file as a list",{
  expect_output(str(read.Structure("StructureEx1.txt")),"List of 3")
  expect_output(str(read.Structure("StructureEx1.txt", haploid = T)),"List of 3")
  expect_output(str(read.Structure("StructureEx2.txt")),"List of 3")
  expect_output(str(read.Structure("StructureEx2.txt", haploid = T)),"List of 3")
})