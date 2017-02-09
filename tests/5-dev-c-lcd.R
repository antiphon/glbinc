# dev the c-side lcd
library(devtools)
load_all(".")
source("tests/0-make-test-set-1.R")

offset <- rep(0,N)

# drop intercept
X <- X[,-1]
index <- index[-1]
# check for global intercept

lambda <- makeLambda(X=X,y=y, offset=offset, index=index, nlambda=10)

un_i <- which(index==0)
pen_i <- which(index!=0)
# penalised groups
gr_i <- table(index)
# sqrt df


# NumericMatrix X, NumericVector y, NumericVector offset
# IntegerVector G0,
# IntegerVector G1,
# NumericVector group_weight,
# NumericVector lambda,
# bool intercept_in_X,
# double alpha,
# double eps,
# int dfmax,
# int maxiter
K <- table(index)
G0 <- as.integer( if(min(index)==0) sum(index==0) else 0 )
G1 <- as.integer( if(min(index)==0) cumsum(K) else c(0, cumsum(K)) )
group_weight <- sqrt( table(index[index!=0]) )
#
out <- glbin_lcd_cpp(X, y, offset,
                     G0, # unpenalised 0-G0
                     G1, # penalised G0-G1
                     group_weight,
                     lambda,
                     TRUE,
                     1,
                     1e-3,
                     ncol(X)+1,
                     100
                     )
