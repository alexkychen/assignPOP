context("Test accuracy.xxx")

test_that("Calculate assignment accuracy for Monte-Carlo results",{
  AccuMC <- accuracy.MC(dir="ResMCtest/")
  expect_output(str(AccuMC),"data.frame")
  expect_true(file.exists("ResMCtest/Rate_of_3_tests_3_pops.txt"))
  plot <- accuracy.plot(AccuMC)
  expect_type(plot, "list")
})

unlink("ResMCtest/Rate_of_3_tests_3_pops.txt")

test_that("Calculate assignment accuracy for K-fold results",{
  AccuKF <- accuracy.kfold(dir="ResKFtest/")
  expect_output(str(AccuKF),"data.frame")
  expect_true(file.exists("ResKFtest/Rate_of_3_tests_3_pops.txt"))
  plot <- accuracy.plot(AccuKF)
  expect_type(plot, "list")
})

unlink("ResKFtest/Rate_of_3_tests_3_pops.txt")

