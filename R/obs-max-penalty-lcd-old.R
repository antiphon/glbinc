#' Compute the maximal penalty for lcd
#'
#'
#' @param X design matrix
#' @param y response 0/1
#' @param index grouping
#' @param offset offset
#'
#' @details
#' See glbin for details of the parameters
#' @export

max_penalty_lcd_old <- function(X, y, index, offset, alpha=1, eps = 1e-4){
  p <- ncol(X)
  un_i <- which(index==0)
  # assume first is intercept
  ff <- glm(y ~ -1 + X[,un_i], family=binomial, offset=offset)
  bo <- ff$coefficients
  bo[is.na(bo)] <- 0
  # penalised groups
  gr_i <- split((1:p)[-un_i],index[-un_i])
  # group counts
  sdfg <- sqrt( sapply(gr_i, length) )
  # groups
  ngroups <- length(gr_i)

  # IRLS fit with only the unpenalised:
  loop <- TRUE
  beta <- c(bo)
  tX <- t(X[,un_i])
  while(loop){
    bold <- beta
    eta <- c(crossprod(tX, beta)) + offset
    mu <- 1.0 / (1.0 + exp(-eta))
    # stability
    #big <- mu < 1e-5 | mu > (1-1e-5)
    w <- 0.25#mu * (1-mu)
    #
    v <- y - mu
    for(i in un_i) {
      xi <- X[,i]
      sx <- sum(w * xi * xi)
      beta[i] <- sum(xi * v)/sx + beta[i]
    }
    loop <- ( eps < max(abs(bold-beta)) )
  }
  # ok:
  # residuals
  #  tX <- t(X[,un_i])
  eta <- c(crossprod(tX, beta)) + offset
  mu <- 1/(1+exp(-eta))
  # stability
  w <- 0.25
  #
  v <- y - mu
  r <- v/w
  #
  gest <- rep(0, ngroups)
  n <- length(y)
  for(j in 1:ngroups){
    gi <- gr_i[[j]]
    inner <- crossprod((X[,gi]) , w * r)
    gest[j] <- sqrt(sum(inner^2))/(n * sdfg[j])
  }

  lambda_max <- max(gest)
  list(lambda_max=lambda_max, beta_un_i = bo)
}
