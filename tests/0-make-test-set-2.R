# Test set 2: big one with offset

if(!exists("Xbig")){
  P <- 1000
  N <- 9000
  beta <- rnorm(P+1, 1, 2) # inc. intercept
  # set first P*0.5 to zero
  ind0 <- 1+1:(P*.5)
  beta[ind0] <- 0
  # With intercept
  X <- cbind(1, matrix(rnorm(N*P), ncol=P))
  prob <- 1/( 1+ exp(-X%*%beta))
  y <- rbinom(N, 1, prob)
  # with offset
  o <- rgamma(N, 1,2)
  Xo <- X
  Xo[,-1] <- Xo[,-1]-mean(o)/9
  probo <- 1/( 1+ exp(-Xo%*%beta - o))
  yo <- rbinom(N, 1, probo)
  Xbig <- 1
  # grouping index vector
  n0 <- length(ind0)
  index <- c(0, rep(1, n0), 1+1:(P-n0))

}
