context("Consistency of the Group-Lasso solution path (package 'gglasso')")

require(gglasso)
library(Matrix)
library(scoop)
library(adelie)

data(bardet)

# x <- scale(as.matrix(bardet$x), FALSE, TRUE)
x <- scale(as.matrix(bardet$x))
y <- bardet$y

group <- rep(1:20,each=5)
group_adelie <- seq(1,100,by=5)

wk <- rep(1,20)

out_adelie <- adelie::grpnet(x, adelie::glm.gaussian(y), groups = group_adelie, alpha = 1, standardize = FALSE, penalty = wk, lmda_path_size = 50)
coef_adelie <- t(as.matrix(coefficients(out_adelie)$betas))
coef_adelie_grp <- rowsum(coef_adelie^2, group)
lambda <- coefficients(out_adelie)$lambda

out_gglasso  <- gglasso::gglasso(x, y, group, lambda = lambda, pf = wk)
coef_gglasso <- out_gglasso$beta
coef_gglasso_grp <- rowsum(coef_gglasso^2, group)

out_quadrupen <- quadrupen::group.lasso(x,y,group,lambda1 = lambda * nrow(x), normalize = FALSE, lambda2=0)
coef_quad <-  as.matrix(out_quadrupen$coefficients)
coef_quad_grp <- rowsum(coef_quad^2, group)

out_scoop <- scoop::group.lasso(x,y,group,lambda = lambda * nrow(x), wk = wk, normalize = FALSE, optim.method = "fista")
coef_scoop <- t(out_scoop@coefficients[,-1])
coef_scoop_grp <- rowsum(coef_scoop^2, group)

matplot(lambda, t(coef_adelie_grp), log="x", type = "l")
matplot(lambda, t(coef_gglasso_grp), log="x", type = "l")
matplot(lambda, t(coef_scoop_grp), log="x", type = "l")
matplot(lambda, t(coef_quad_grp), log="x", type = "l")


plot(lambda, colSums(coef_adelie_grp - coef_gglasso_grp), log="x")
plot(lambda, colSums(coef_adelie_grp - coef_quad_grp), log="x")
plot(lambda, colSums(coef_adelie_grp - coef_scoop_grp), log="x")
plot(lambda, colSums(coef_gglasso_grp - coef_quad_grp), log="x")
plot(lambda, colSums(coef_gglasso_grp - coef_scoop_grp), log="x")
plot(lambda, colSums(coef_quad_grp - coef_scoop_grp), log="x")

plot(lambda, colSums(coef_adelie - coef_gglasso), log="x")
plot(lambda, colSums(coef_adelie - coef_quad), log="x")
plot(lambda, colSums(coef_adelie - coef_scoop), log="x")
plot(lambda, colSums(coef_gglasso - coef_quad), log="x")
plot(lambda, colSums(coef_gglasso - coef_scoop), log="x")
plot(lambda, colSums(coef_quad - coef_scoop), log="x")

image(Matrix(coef_adelie - coef_gglasso))
image(Matrix(coef_adelie - coef_quad))
image(Matrix(coef_adelie - coef_scoop))
image(Matrix(coef_scoop  - coef_quad))

expect_equal(coef_quad,
             coef_adelie, check.attributes = FALSE, tolerance = 1e-2)

sum(out_gglasso$beta - as.matrix(out_quadrupen$coefficients)^2)

get_grplasso <- function(x,y,group,intercept,normalize=TRUE) {

  out_gglasso <- gglasso(x, y, group, loss="ls", pf = rep(1,length(unique(group))), intercept=intercept, nlambda = 50, lambda.factor = 1e-2)
  lambda <- out_gglasso$lambda
    
  out_quadrupen <- group.lasso(x,y,group,lambda1 = lambda * nrow(x), intercept=intercept,normalize=FALSE, lambda2=0, control = list(threshold=1e-6))

  return(list(quad = quad, enet = t(enet)))
}

test_that("Elastic-net is correct w.r.t a reference solution", {

  tol <- 1e-2
  
  data(bardet)
  group <- rep(1:20,each=5)
  
  x <- scale(as.matrix(bardet$x), FALSE, TRUE)
  y <- bardet$y

  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
  set.seed(seed)

  beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
  n <- 100
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

})


test_that("Elastic-net is correct w.r.t a reference solution - FISTA", {

  tol <- 1e-1
  
  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)
  
  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE, method="fista")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)
  
  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE, method="fista")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## Run the tests...
  without.intercept <- get.enet(x,y,intercept=FALSE, method="fista")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)
  
  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE, method="fista")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)
  

})
