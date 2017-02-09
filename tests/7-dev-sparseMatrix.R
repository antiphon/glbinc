# Dev the c-version of bcgd SSR
library(devtools)
load_all(".")

library(Matrix)
#Q <- readRDS("~/Dropbox/work/joint-analysis-of-rainforest-interactions/examples-v2/fit-bci/bci-v1-Q.rds")
#x <- Q$X
set.seed(1)
i <- sample(1:10, 20, T)
j <- sample(1:20, 20,T)
v <- round(runif(20),2)
x <- sparseMatrix(i = i, j=j, x = v, dims = c(10,20))

o <- sparse_c(x)

