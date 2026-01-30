testDataBreg <- readRDS("dataTest-boundedReg.rds")

tol <- 1e-3

test_that("Bounded regression with lambda2 = 0 - test on the documentation example", {
  
  res <- bounded.reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 0
    )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_0$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 5 - test on the documentation example", {
  
  res <- bounded.reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 5
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_5$coefficients, tolerance = tol)
  
})

test_that("Bounded regression with lambda2 = 10 + S  - test on the documentation example", {
  
  res <- bounded.reg(
    testDataBreg$x,
    testDataBreg$y,
    lambda2 = 10, 
    struct = testDataBreg$S
  )
  expect_equal(res$coefficients, testDataBreg$breg_lambda2_10_S$coefficients, tolerance = tol)
  
})
