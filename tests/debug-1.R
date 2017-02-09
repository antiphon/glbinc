# debug something encountered with multipaper experiment 1, beta0 = NA

library(devtools)
load_all(".")

a <- readRDS("../../../joint-analysis-of-rainforest-interactions/examples-v4-right-Q/experiment-1/perkele.rds")
X <- a$X
y <- a$y
o <- a$o
index <- a$index
lambda <- a$lvec

f <- glbin_lcd_c(X, y, o, index, lambda = lambda)

head(f$beta)
