# test the c-side lcd sparse
library(devtools)
library(Matrix)
load_all(".")

set.seed(2)
source("tests/0-make-test-set-2-sparse.R")
offset <- o
X <- Xo
y <- yo
# drop intercept
X <- X[,-1]
index <- index[-1]
# Fits

eps <- 0.001
verb <- 1
nlambda <- 100
lmin <- 0.001
#
t1 <- system.time( f1 <- glbin_lcd_c(as.matrix(X), y,
                                     offset = offset, index = index, verb=verb,
                                     eps=eps, nlambda = nlambda, lambda.min = lmin, AIC_stop=3) )

t2 <- system.time( f2 <- glbin_lcd_c_sparse(Matrix(X, sparse=T), y,
                                            offset = offset, index = index, verb=verb,
                                            eps=eps, nlambda = nlambda, lambda.min = lmin, AIC_stop=3) )

k1 <- 4
m <- 10
e <- rbind(dense=f1$beta[1:m,k1],
           sparse=f2$beta[1:m,k1],
           true=beta[1:m])
print(e)

# check penalties same
#print(rbind(lam_d=head(f1$lambda), lam_s=head(f2$lambda)))

cat("\n**TIME\n")
print(rbind(dense=t1, sparse=t2))
