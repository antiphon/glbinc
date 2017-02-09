# Dev quadratic approximation v2 blocks, like in H&Y 2009, offset now in
library(devtools)
load_all(".")


set.seed(2)
source("tests/0-make-test-set-2.R")

verb <- 1
lmin <- 0.001
eps <- 1e-3

X <- Xo
y <- yo
o <- o


if(!exists("f0"))
  t0 <- system.time(  f0 <- glbin_bcgd(X, y, offset=o, eps=eps, nlambda=100, index = index, verb=verb, lambda.min = lmin, SSR = FALSE)  )

t1 <- system.time(f1 <- glbin_lcd(X, y, offset=o, index=index, eps=1e-2,
                                  lambda.min=lmin, nlambda=100, verb=verb, std=T,
                                  alpha = 1,
                                  stability_eps = 1e-5)  )

if(!exists("f2")){
  t2 <- system.time(f2 <- grpreg::grpreg(X[,-1], group = index[-1], y, offset=o, family="binomial", lambda = f1$lambda))
  f2v <- f1
  f2v$beta <- f2$beta
  f2v$lambda <- f2$lambda
}
#
k <- which.min(f0$aic)/2
k1 <- which.min(f1$aic)/2
print(c(aic0=k, aic1=k1))
print(rbind(bcgd=f0$beta[, k], iwls=f1$beta[,k], grpreg=f2v$beta[,k], true=beta))

#
par(mfrow=c(3,2))
plot(beta, f0$beta[,k], main="BCGD", asp=1); abline(0,1); plot(f0, log="x")
plot(beta, f1$beta[,k], main="IWLS", asp=1); abline(0,1); plot(f1, log="x")
plot(beta, f2$beta[,k], main="grpreg", asp=1); abline(0,1); plot(f2v, log="x")
#
print(rbind(c=t0,R=t1, grpreg=t2))
