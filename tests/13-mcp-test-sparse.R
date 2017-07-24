# MCP sparse
#devtools::load_all(".")
#library(devtools)
#load_all(".")
set.seed(1)
source("tests/0-make-test-set-1.R")

library(grpreg)
pen <- c("grLasso", "grMCP", "grSCAD"); names(pen)<-pen
f <- lapply(pen, function(v) grpreg(as.matrix(X)[,-1], y, group = index[-1], family = "binomial", penalty = v))
g <- lapply(pen, function(v) glbin_lcd_c_sparse(Matrix(X, sparse=TRUE), y, index = index, add_intercept = F, verb = 1, penalty=v, lambda = f[[v]]$lambda))

par(mfrow=c(2,3))
lapply(pen, function(n)plot(f[[n]], main=paste("grpreg", n)))
lapply(pen, function(n)plot(g[[n]], main=paste("glbin", n)))
