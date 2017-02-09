# Dev quadratic approximation v2 blocks, like in H&Y 2009
library(devtools)
load_all(".")


set.seed(20)
source("tests/0-make-test-set-1.R")

verb <- 1
lmin <- 0.001
eps <- 1e-3

t0 <- system.time(  f0 <- glbin_lcd(X, y, eps=eps, nlambda=100, index = index, verb=verb, lambda.min = lmin)  )
