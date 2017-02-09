# Tets the lambda.max function

library(devtools)
library(Matrix)
load_all(".")

set.seed(4)
source("tests/0-make-test-set-1-sparse.R")
off <- o
X <- Xo
y <- yo

# Q <- readRDS("~/Dropbox/work/joint-analysis-of-rainforest-interactions/examples-v2/fit-bci/bci-v1-Q.rds")
# X <- Q$X
# off <- Q$offset
# y <- Q$y
# index <- Q$pen_index

Xs <- scale(as.matrix(X))
Xs[,1] <- 1
# drop intercept
X <- X[,-1]
index <- index[-1]

center <- colMeans(X)
scale <- apply(X, 2, sd)


t1 <- system.time(lmax <- max_penalty_lcd(Xs, y, index, off, TRUE))
t2 <- system.time(lold <- max_penalty_lcd_old(X, y, index, off, TRUE))
t3 <- system.time(lnew <- max_penalty(Xs[,-1], y, index, off, add_intercept = FALSE, verb=0))
t4 <- system.time(lnews <- max_penalty(Matrix(X, sparse=TRUE), y, center = center, scale = scale, index, off, verb=0))
t4 <- system.time(lnewalt <- max_penalty(X, y, index, off, verb=0, use_glm = TRUE))

print(c(glm=lmax, glmold=lold$lambda_max, dense=lnew, sparse=lnews, glm2 = lnewalt))
#print(c(lmax, lnew, lnews))
#print(c(dense=lnew, sparse=lnews))
print(rbind(t1,t2,t3,t4))
#print(rbind(t1,t3,t4))

#lv <- makeLambda(100, 0.001, X=X[,-1], y=y, offset=off, index=index[-1], alpha = 1)

