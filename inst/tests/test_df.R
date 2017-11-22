rm(list=ls())
library(quadrupen)
library(reshape2)
source("crit-class.R")

## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
n <- 80
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

vec.lambda2 <- 10^seq(.5,-3,len=50)
crit <- sapply(vec.lambda2, function(lambda2) {
  pen.criteria(elastic.net(x,y,lambda2=lambda2,nlambda1=50),pen=log(nrow(x)),sigma=10)@criterion
})
colnames(crit) <- vec.lambda2
crit <- melt(crit, varnames=c("lambda1", "lambda2"))

d <- ggplot(data=crit, aes(x=lambda1, y=lambda2, z=value))
d <- d + geom_tile(aes(fill=value)) + stat_contour(size=0.2, binwidth=diff(range(crit$value))/10)
d <- d + scale_x_continuous(trans=log10_trans()) + xlab(expression(log[10](lambda[1])))
d <- d + scale_y_continuous(trans=log10_trans()) + ylab(expression(log[10](lambda[2])))
d <- d + annotation_logticks()
print(d)

vec.lambda2 <- 10^seq(.5,-3,len=50)
df <- sapply(vec.lambda2, function(lambda2) {
  elastic.net(x,y,lambda2=lambda2,nlambda1=50)@df
})
colnames(df) <- vec.lambda2
rownames(df) <- elastic.net(x,y,nlambda1=50)@lambda1
df <- melt(df, varnames=c("lambda1", "lambda2"))

d <- ggplot(data=df, aes(x=lambda1, y=lambda2, z=value))
d <- d + geom_tile(aes(fill=value)) + stat_contour(size=0.2, binwidth=diff(range(df$value))/10)
d <- d + scale_x_continuous(trans=log10_trans()) + xlab(expression(log[10](lambda[1])))
d <- d + scale_y_continuous(trans=log10_trans()) + ylab(expression(log[10](lambda[2])))
d <- d + annotation_logticks()
print(d)

## randomized BIC
## BIC.samples <- replicate(100, {
##   train=sample(1:n,floor(n/2))
##   out <- elastic.net(x[train, ],y[train],lambda2=lambda2,lambda1=ref@lambda1)
##   return(pen.criteria(out, crit="BIC"))
## })

## library(plyr)
## data <- melt(BIC.samples)
## colnames(data) <- c("lambda1", "replicate", "BIC")
## BIC.ave <- ddply(data, .(lambda1), function(x) data.frame(mean=mean(x$BIC), sd=sd(x$BIC)))
## BIC.ave[, 1] <- as.numeric(as.character(BIC.ave[, 1]))
## data[, 1] <- as.numeric(as.character(data[, 1]))
## plot.BIC <- ggplot(data=BIC.ave,aes(x=lambda1,y=mean)) +
##     geom_line(aes(x=lambda1, y=BIC, group=replicate, alpha=0.5, colour="red"), data=data) +
##        geom_smooth(aes(ymin=mean-sd, ymax=mean+sd), data=BIC.ave, alpha=0.2, stat="identity", width=3) +
##            labs(title="", y="", x="") + theme(legend.direction = "horizontal") + scale_x_log10() +
##               theme(legend.position="none")

## print(plot.BIC)
