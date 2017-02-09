# debuggin
library(devtools)
load_all(".")
Q <- readRDS("~/Dropbox/work/joint-analysis-of-rainforest-interactions/examples-v2/experiment-2/exp2-50/data/test-Q.rds")
ok <- Q$subset
X <- Q$X[ok,]
y <- Q$y[ok]
o <- Q$offset[ok]
index <- Q$pen_index
# if we add the intercept, drop the first to reduce singularities
X <- X[,-1]
index <- index[-1]


# go
if(1){ # problem
fit <- glbin_lcd_c_sparse(X=X, y=y, index, offset=o,
                          lambda.min = 0.001, nlambda = 100,
                          eps = 1e-3,
                          dfmax = ncol(X)+1, verb = 1, AIC_stop = 3,
                          add_intercept = TRUE)

}

# problem was in in lambdamax
