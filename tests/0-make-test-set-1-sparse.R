# Make test set 1 with sparse X

P <- 5
N <- 200
beta <- rnorm(P+1, 1, 2) # inc. intercept
# set first 30% coefficient to zero
ind0 <- 1+1:(.3*P)
beta[ind0] <- 0
beta0 <- beta
# covariates With intercept
X <- matrix(rnorm(N*P), ncol=P)
# sparsify X
si <- sample(length(X), length(X)/5)
X[si] <- 0
# add intersept
X <- cbind(1,X)

prob <- 1/( 1+ exp(-X%*%beta))
# with offset
o <- rgamma(N,1,1)
Xo <- X
Xo[,-1] <- X[,-1]-mean(o)/9
# sparse!
XX <- (Xo[,-1])
XX[si] <- 0
Xo[,-1] <- XX
probo <- 1/(1 + exp(-Xo%*%beta - o))

y <- rbinom(N, 1, prob)
table(y)
# with offset
yo <-rbinom(N, 1, probo)
print(table(yo))

# grouping index vector
n0 <- length(ind0)
index <- c(0, 0, rep(1, n0-1), 1+1:(P-n0))
index <- c(0, 2:ncol(X)-1)
# add bit more non-penalised
index <- index - 3
index[1:4] <- 0
