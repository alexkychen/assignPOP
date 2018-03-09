context("Test assign.matrix")

test_that("run assign.matrix",{
  expect_output(str(assign.matrix(dir = "ResMCtest/")),"xtabs")
})