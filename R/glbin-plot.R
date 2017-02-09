#' Plot method
#'
#' @export

plot.glbin <- function(x, color=TRUE, ...) {
  l <- x$lambda
  nl <- length(l)
  g <- ncol(x$beta)
  ry <- range(x$beta, na.rm=TRUE)
  plot(NA, xlim=rev(range(l)), ylim=ry, xlab="penalty", ylab="beta", ...)
  for(i in 1:nrow(x$beta)){
    v <- rep(NA, nl)
    v[1:g] <- x$beta[i,]
    lines(l, v, type="l", col=if(color) i else 1)
  }
}
