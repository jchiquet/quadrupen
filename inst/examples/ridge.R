library(quadrupen)

beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)
p <- ncol(x)


## neighborhood prior
C <- bandSparse(p,k=0:1,diagonals=list(rep(1,p),rep(-1,p-1)))
L1 <-  t(C) %*% C

## clique prior
Ioo <- matrix(1,25,25)
Iww <- matrix(1,10,10)
A <- bdiag(Ioo,Iww,Ioo,Iww,Ioo)
diag(A) <- 0
L2 <- -A
diag(L2) <- colSums(A) + 1e-2

ridge.classical <- ridge(x,y,intercept=F)
plot(ridge.classical)
ridge.neighborhood <- ridge(x,y, struct=L1, lambda_max=1000)
plot(ridge.neighborhood)
ridge.block <- ridge(x,y, struct=L2, lambda_max=1000)
plot(ridge.block)


beta_ridge <- as.numeric(ridge(x,y,lambda=1)$coefficients)
beta_ridge_neighbor <- as.numeric(ridge(x,y,lambda=10,struct=L1)$coefficients)
beta_ridge_block <- as.numeric(ridge(x,y,lambda=10,struct=L2)$coefficients)
beta_enet_block  <- as.numeric(elastic_net(x,y,lambda2=5,lambda1=20,struct=L2)$coefficients)

par(mfrow=c(2,2))
plot(beta_ridge, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_ridge_neighbor, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_ridge_block, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
plot(beta_enet_block, type="b", ylim=c(-1,1))
points(beta, pch="+",col="red")
