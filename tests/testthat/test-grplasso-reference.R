context("Consistency of the Group-Lasso solution path")

testData <- readRDS("dataTest-GroupLasso.rds")

tol <- 1e-3

get_grplasso <- function(x, y, group, lambda, intercept, normalize) {
  
  out_quadr   <- quadrupen::group.lasso(x, y, group, type = "l2", lambda1 = lambda, 
                                        intercept = intercept, normalize = normalize, lambda2 = 0,
                                        control = list(maxfeat = ncol(x)))
  coef_quadr  <-  as.matrix(out_quadr$coefficients)
  group_quadr <- rowsum(coef_quadr^2, group)
  inter_quadr <- out_quadr$intercept
  
  res <- list(coef = coef_quadr, group = group_quadr, intercept = inter_quadr, lambda = lambda)
  res
}

test_that("Group-Lasso with lambda2 = 0, intercept and normalization, FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_inter_norm$lambda, TRUE, TRUE)
  
  expect_equal(quad, testData$grplasso_inter_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, intercept and no normalization - FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_inter_nonorm$lambda, TRUE, FALSE)
  
  expect_equal(quad, testData$grplasso_inter_nonorm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, no intercept and normalization - FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_norm$lambda, FALSE, TRUE)
  
  expect_equal(quad, testData$grplasso_nointer_norm, check.attributes = FALSE, tolerance = tol)
}
)

test_that("Group-Lasso with lambda2 = 0, no intercept and no normalization - FISTA - test on the documentation example", {
  quad <- get_grplasso(testData$x, testData$y, testData$group, testData$grplasso_nointer_nonorm$lambda, FALSE, FALSE)
  
  expect_equal(quad, testData$grplasso_nointer_nonorm, check.attributes = FALSE, tolerance = tol)
}
)
