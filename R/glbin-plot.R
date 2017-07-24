#' Plot method
#'
#' @export

plot.glbin <- function(x, color=TRUE, only_pen = FALSE, ...) {
  l <- x$lambda
  nl <- length(l)
  beta <- x$beta
  if(only_pen) beta <- beta[x$index != 0, ]
  g <- ncol(beta)
  ry <- range(beta, na.rm=TRUE)
  plot(NA, xlim=rev(range(l)), ylim=ry, xlab="penalty", ylab="beta", ...)
  for(i in 1:nrow(beta)){
    v <- rep(NA, nl)
    v[1:g] <- beta[i,]
    lines(l, v, type="l", col=if(color) i else 1)
  }
}
