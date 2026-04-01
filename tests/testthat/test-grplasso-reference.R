context("Consistency of the Group-Lasso solution path")

library(Matrix)
library(scoop)
# library(adelie)

get_grplasso_scoop <- function(x,y,group,intercept,normalize=TRUE) {
  
  out_scoop   <- scoop::group.lasso(x, y, group, intercept = intercept, normalize = normalize, lambda.min = 1e-2, n.lambda = 50)
  if (intercept) {
    coef_scoop  <- t(as.matrix(out_scoop@coefficients[, -1]))
    inter_scoop <- out_scoop@coefficients[, 1]
  } else {
    coef_scoop  <- t(as.matrix(out_scoop@coefficients))
    inter_scoop <- rep(0,50)
  }
  group_scoop <- rowsum(coef_scoop^2, group)
  
  out_quadr   <- quadrupen::group.lasso(x, y, group, lambda1 = out_scoop@lambda, intercept = intercept, normalize = normalize, lambda2 = 0)
  coef_quadr  <-  as.matrix(out_quadr$coefficients)
  group_quadr <- rowsum(coef_quadr^2, group)
  inter_quadr <- out_quadr$intercept
  
  return(
    list(quadr = list(coef = coef_quadr, group = group_quadr, intercept = inter_quadr), 
         scoop = list(coef = coef_scoop, group = group_scoop, intercept = inter_scoop))
  )
}

# get_grplasso <- function(x,y,group,intercept,normalize=TRUE) {
# 
#   group_adelie <- c(1, 1 + cumsum(tail(tabulate(group), n = -1)))
#   
#   out_adelie <- adelie::grpnet(x, adelie::glm.gaussian(y), groups = group_adelie, alpha = 1, intercept = intercept, standardize = normalize, lmda_path_size = 50)
#   lambda <- coefficients(out_adelie)$lambda
#   coef_adelie <- t(as.matrix(predict(out_adelie, x, lambda, "coefficients")$betas))
#   group_adelie <- rowsum(coef_adelie^2, group)
#   inter_adelie <- as.numeric(coefficients(out_adelie)$intercepts)
#   if (!intercept) inter_adelie <- rep(0, length(inter_adelie))
# 
#   out_quadrupen <- group.lasso(x, y, group, lambda1 = lambda * nrow(x), intercept = intercept, normalize = normalize, lambda2 = 0)
#   coef_quad <-  as.matrix(out_quadrupen$coefficients)
#   group_quad <- rowsum(coef_quad^2, group)
#   inter_quad <- out_quadrupen$intercept
# 
#   return(
#     list(quadru = list(coef = coef_quad, group = group_quad, intercept = inter_quad), 
#          adelie = list(coef = coef_adelie, group = group_adelie, intercept = inter_adelie))
#      )
# }

test_that("Group-Lasso is correct w.r.t a reference solution", {

  bardet <- readRDS("bardet.rds")
  group <- rep(1:20, each = 5)
  
  tol <- 1e-2

  x <- as.matrix(bardet$x)
  y <- bardet$y

  ## Run the tests...
  with.intercept <- get_grplasso_scoop(x,y,group,intercept=TRUE)
  expect_equal(with.intercept$quadr,
              with.intercept$scoop, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get_grplasso_scoop(x,y,group,intercept=FALSE)
  expect_equal(without.intercept$quadr,
               without.intercept$scoop, check.attributes = FALSE, tolerance = tol)
   
  with.intercept <- get_grplasso_scoop(x,y,group,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quadr,
              with.intercept$scoop, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get_grplasso_scoop(x,y,group,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quadr,
              without.intercept$scoop, check.attributes = FALSE, tolerance = tol)

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
  with.intercept <- get_grplasso_scoop(x,y,group,intercept=TRUE)
  expect_equal(with.intercept$quadr,
              with.intercept$scoop, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get_grplasso_scoop(x,y,group,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quadr,
              with.intercept$scoop, check.attributes = FALSE, tolerance = tol)

})
