# test the c-side glm for computing the lambda max
library(Matrix)
library(devtools)
load_all(".")
set.seed(10)
source("tests/0-make-test-set-1.R")

offset <- o
X <- Xo
y <- yo

# drop intercept
X <- X[,-1]
# Fits

eps <- 0.00001

center <- apply(X,2, mean)
scale <- apply(X,2, sd)

if(!exists("f2")) {
f0 <- glm(y~1+X, offset=offset, family="binomial")
}


f1 <- glm_binom_std_c(X, y, offset, center, scale, TRUE, 100, eps=eps, 1000, w_limit = 0.1)
f2 <- glm_binom_std_sparse_c(Matrix(X, sparse=TRUE), y, offset, center, scale, TRUE, 100, eps=eps, 1000)
# backscale
beta_res <- c(f1$beta0, f1$beta)
beta_res[-1] <- beta_res[-1]/scale
beta_res[1] <- beta_res[1] - sum(center * beta_res[-1])

beta_res2 <- c(f2$beta0, f2$beta)
beta_res2[-1] <- beta_res2[-1]/scale
beta_res2[1] <- beta_res2[1] - sum(center * beta_res[-1])

# Checks

z <- rbind(glm=f0$coefficients, c=beta_res, cs=beta_res2, true=beta)
print(z)

