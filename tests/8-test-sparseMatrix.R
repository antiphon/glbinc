# test the c-side lcd sparse version
library(devtools)
load_all(".")
set.seed(2)
source("tests/0-make-test-set-1-sparse.R")

offset <- o
X <- Xo
y <- yo

# drop intercept
X <- X[,-1]
index <- index[-1]

# as sparse
library(Matrix)
Xs <- Matrix(X)
#
# offset <- o*0
# X<-X
# y<-y

# Fits

eps <- 0.001

if(!exists("f2")) {
  t2 <- system.time( f2 <- glbin_lcd_c(X, y, index) )
}

# my c
t1 <- system.time( f1 <- glbin_lcd_c(X, y,
                                     offset = offset, index = index, verb=1,
                                     eps=eps, nlambda = 100, AIC_stop=5) )

# Checks

f2v <- f1
f2v$beta <- f2$beta
f2v$lambda <- f2$lambda

k1 <- 20#which.min(f1$aic)
k2 <- k1#which.min(AIC(f2))
k3 <- k2#which.min(f3$aic)
m <- min(10, ncol(X))

par(mfrow=c(3,3))
plot(f1, main="c"); plot(beta, f1$beta[,k1]); abline(0,1); plot(f1$logLik)
plot(f2v, main="grpreg"); plot(beta, f2$beta[,k2]); abline(0,1); plot(-f2$loss)
plot(f3, main="R"); plot(beta, f3$beta[,k3]); abline(0,1); plot(f3$logLik)

e <- rbind(grpreg=f2$beta[1:m,k2],
           myc=f1$beta[1:m,k1],
           myR=f3$beta[1:m,k3],
           true=beta[1:m])
print(e)

cat("\n**TIME\n")
print(rbind(grep=t2,
            myc=t1,
            myR=t3))
