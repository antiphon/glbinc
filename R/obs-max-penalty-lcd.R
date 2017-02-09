#' Compute the maximum penalty
#'
#' Copied from grpreg
#'
#' @param X design matrix without intercept
#' @param y response 0/1
#' @param index grouping
#' @param offset offset
#' @param alpha lasso vs ridge weight
#' @param add_intercept TRUE
#' @details
#' See glbin for details of the parameters
#' @export

max_penalty_lcd <- function(X, y, index, offset, add_intercept =TRUE, alpha=1){
  warning("obsolete, use max_penalty")
  p <- ncol(X)
  un_i <- which(index==0)
  pen_i <- which(index!=0)
  #
  if(any(index==0)){
    if(add_intercept)
      fit <- glm(y ~ X[,un_i], family=binomial, offset=offset)
    else fit <- glm(y ~ -1 + X[,un_i], family=binomial, offset=offset)
  }
  else{ # no other choice
    if(!add_intercept) stop("model has no terms unpenalised terms or intercept.")
    fit <- glm(y ~ 1, family=binomial, offset=offset)
  }
  w <- fit$weights
  if (max(w) < 1e-4) stop("Unpenalized portion of model is already saturated; exiting...")
  r <- residuals(fit, "working") * w
  n <- length(y)
  # penalised groups
  gr_i <- split((1:p)[pen_i], index[pen_i])
  # sqrt df
  sdfg <- sqrt( sapply(gr_i, length) )
  # groups
  ngroups <- length(gr_i)
  # the gradients per group
  gest <- sapply(1:ngroups, function(g){
    ig <- gr_i[[g]]
    z <- crossprod(X[,ig],r)
    sqrt(sum(z^2))/sdfg[g]/n
  })
  max(gest)/alpha
}

