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
  
  bardet <- readRDS("tests/testthat/bardet.rds")
  group <- rep(1:20, each = 5)
  
  tol <- 1e-2
  
  x <- as.matrix(bardet$x)
  y <- bardet$y
  
  ## Run the tests...
  with.intercept <- get_cooplasso_scoop(x,y,group,intercept=TRUE)
  expect_equal(with.intercept$quadr,
               with.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  
  # without.intercept <- get_cooplasso_scoop(x,y,group,intercept=FALSE)
  # expect_equal(without.intercept$quadr,
  #             without.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  # 
  # with.intercept <- get_cooplasso_scoop(x,y,group,intercept=TRUE,normalize=FALSE)
  # expect_equal(with.intercept$quadr,
  #             with.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  # 
  # without.intercept <- get_cooplasso_scoop(x,y,group,intercept=FALSE,normalize=FALSE)
  # expect_equal(without.intercept$quadr,
  #             without.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  
  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
  set.seed(seed)
  
  beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
  group <- rep(1:5, c(25,10,25,10,25)) 
  n <- 100
  p <- length(beta)
  
  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1
  
  x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)
  
  ## Run the tests...
  with.intercept <- get_cooplasso_scoop(x,y,group,intercept=TRUE)
  expect_equal(with.intercept$quadr,
               with.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  
  # without.intercept <- get_cooplasso_scoop(x,y,group,intercept=FALSE)
  # expect_equal(without.intercept$quadr,
  #             without.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  # 
  # with.intercept <- get_cooplasso_scoop(x,y,group,intercept=TRUE,normalize=FALSE)
  # expect_equal(with.intercept$quadr,
  #             with.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  # 
  # without.intercept <- get_cooplasso_scoop(x,y,group,intercept=FALSE,normalize=FALSE)
  # expect_equal(without.intercept$quadr,
  #             without.intercept$scoop, check.attributes = FALSE, tolerance = tol)
  
})
