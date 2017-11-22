## Mutivariate Gaussian data
mu <- 3
sigma <- 10
n <- 100
beta <- rep(c(0,-2,2),c(80,10,10))
cor <- 0.8
eps <- 0.25 # correlation between relevant and irrelevant variable
S11 <- toeplitz(cor^(0:(80-1))) ## Toeplitz correlation for irrelevant variable
S22 <- matrix(cor,10,10) ## bloc correlation bewteen relevant variables
diag(S22) <- 1
Sigma <- bdiag(S11,S22,S22) + eps

## Build a vector of label for true nonzeros
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- c("relevant")
labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))

## uncomment if you want insights on the correlation structure
## image(Sigma)
x <- as.matrix(matrix(rnorm(100*n),n,100) %*% chol(Sigma))
y <- mu + x %*% beta + rnorm(n, 0, sigma)

## Test simple and double cross-validation
cv.simple <- crossval(x, y, "lasso")
plot(cv.simple)
fit <- lasso(x,y)

## plot the solution path
plot(fit)
## a quick summary of the fit
print(fit)
## AIC and BIC
criteria(fit)

## Call to stability selection function, 200 subsampling
stab <- stability(x,y, "lasso", subsamples=200, min.ratio=1e-2)
## plot the stability path
plot(stab, labels=labels, nvar=20)
## a quick summary of the fit
print(stab)
