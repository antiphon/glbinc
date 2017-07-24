# MCP std version
devtools::load_all(".")
#library(devtools)
#load_all(".")
set.seed(1)
source("tests/0-make-test-set-1.R")

library(grpreg)
pen <- c("grLasso", "grMCP", "grSCAD"); names(pen)<-pen


alpha <- 0.99

if(!exists("tf"))tf <- system.time( f <- lapply(pen, function(v) grpreg(as.matrix(X)[,-1], y, group = index[-1], family = "binomial", penalty = v, alpha = alpha)) )

tg <- system.time( g <- lapply(pen, function(v) glbin_lcd_c_std(X[,-1], y, index = index[-1], add_intercept = T, verb = 0, penalty=v, lambda = f[[v]]$lambda, alpha = alpha)) )

tgs <- system.time( gs <- lapply(pen, function(v) glbin_lcd_c_sparse(Matrix(X[,-1], sparse=TRUE), y, index = index[-1], add_intercept = T, verb = 0, penalty=v, lambda = f[[v]]$lambda, alpha = alpha)) )

par(mfrow=c(3,3))
lapply(pen, function(n)plot(f[[n]], main=paste("grpreg", n)))
lapply(pen, function(n)plot(g[[n]], main=paste("glbin std", n)))
lapply(pen, function(n)plot(gs[[n]], main=paste("glbin sparse", n)))


print(rbind(grpreg=tf, tg_std=tg, tg_sparse=tgs))
