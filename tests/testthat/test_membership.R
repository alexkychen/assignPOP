context("Test membership.prob.plot")

test_that("Plot membership probability",{
  plot <- membership.plot(dir="ResKFtest/", style=1, non.genetic=T)
  expect_s3_class(plot, "ggplot")
})