# Make test set 1

P <- 9
N <- 1000
beta <- rnorm(P+1, 1, 2) # inc. intercept
# set first 30% to zero
ind0 <- 1+1:(.3*P)
beta[ind0] <- 0
beta0 <- beta
# With intercept
X <- cbind(1, matrix(rnorm(N*P), ncol=P))

prob <- 1/( 1+ exp(-X%*%beta))
# with offset
o <- rgamma(N,5,1)
Xo <- X
Xo[,-1] <- X[,-1]-mean(o)/9
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
