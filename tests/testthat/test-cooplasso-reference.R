context("Consistency of the Cooperative-Lasso solution path (package 'gglasso')")

library(Matrix)
library(scoop)
# library(adelie)

get_cooplasso_scoop <- function(x,y,group,intercept,normalize=TRUE) {
  
  out_scoop   <- scoop::coop.lasso(x, y, group, intercept = intercept, normalize = normalize, lambda.min = 1e-2, n.lambda = 50)
  if (intercept) {
    coef_scoop  <- t(as.matrix(out_scoop@coefficients[, -1]))
    inter_scoop <- out_scoop@coefficients[, 1]
  } else {
    coef_scoop  <- t(as.matrix(out_scoop@coefficients))
    inter_scoop <- rep(0,50)
  }
  group_scoop <- rowsum(coef_scoop^2, group)
  
  out_quadr   <- quadrupen::group.lasso(x, y, group, type = "coop", lambda1 = out_scoop@lambda, intercept = intercept, normalize = normalize, lambda2 = 0)
  coef_quadr  <-  as.matrix(out_quadr$coefficients)
  group_quadr <- rowsum(coef_quadr^2, group)
  inter_quadr <- out_quadr$intercept
  
  return(
    list(quadr = list(coef = coef_quadr, group = group_quadr, intercept = inter_quadr), 
         scoop = list(coef = coef_scoop, group = group_scoop, intercept = inter_scoop))
  )
}

test_that("Coop-Lasso is correct w.r.t a reference solution", {
  
  bardet <- readRDS("bardet.rds")
  group <- rep(1:20, each = 5)
  
  tol <- 1e-2
  
  x <- as.matrix(bardet$x)
  y <- bardet$y
  
  ## Run the tests...
  with.intercept <- get_cooplasso_scoop(x,y,group,intercept=TRUE)
  expect_equal(with.intercept$quadr,
               with.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  
  without.intercept <- get_cooplasso_scoop(x,y,group,intercept=FALSE)
  expect_equal(without.intercept$quadr,
              without.intercept$scoop, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get_cooplasso_scoop(x,y,group,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quadr,
              with.intercept$scoop, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get_cooplasso_scoop(x,y,group,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quadr,
              without.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  
})
