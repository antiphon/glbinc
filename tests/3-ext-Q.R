# Dev the c-version of bcgd SSR
library(devtools)
load_all(".")


#Q <- readRDS("~/Dropbox/work/joint-analysis-of-rainforest-interactions/examples-v2/fit-bci/bci-v1-Q.rds")
Q <- readRDS("~/Dropbox/work/joint-analysis-of-rainforest-interactions/examples-v2/example-fit-1/testQ.rds")
# X <- Q$X
# Xz <- as.matrix(X)
# yo <- Q$y
# o <- Q$offset
# index <- Q$pen_index

ok <- Q$subset
X <- as.matrix( Q$X[ok,] )
y <- Q$y[ok]
o <- Q$offset[ok]
index <- Q$pen_index

# drop the intercept
X <- X[,-1]
index <- index[-1]

verb <- 1
lmin <- 0.001
eps <- 1e-3
#t0 <- system.time(  f0 <- glbin_bcgd_2(X, yo, offset = o, eps=eps, nlambda=50, index = index, verb=verb, lambda.min = lmin, SSR=TRUE)  )
Xs <- scale(X)
lambda_vec <- makeLambda(nlambda = 100, lambda.min = 0.01,
                         X=Xs, y=y, offset =o, index = index)

fit <- glbin_lcd_c(X=X, y=y, index, offset=o,
                   #lambda.min = lambda.min, nlambda = nlambda,
                   AIC_stop = 5,
                   verb = 1,
                   add_intercept = T)

#
lmax1 <- max_penalty_lcd(as.matrix(X), y, index, o)
lmax <- max_penalty(X, y, index, o, verb = 10, eps = 1e-3)
