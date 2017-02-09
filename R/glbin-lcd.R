#' Block coordinate gradient descent for logistic group lasso
#'
#' Like in Breheny and Huang 2009
#'
#' @param X covariate/design matrix. intercept is taken to be the first column!
#' @param y response 0/1 vector
#' @param index grouping. set 0 or NA for unpenalised terms. if not given, assumes X[,1] is intercept
#' @param offset offset terms
#' @param eps convergence criterion
#' @param lambda vector of lambdas
#' @param lambda.min fraction of max lambda to go down to
#' @param nlambda number of penalties
#' @param dfmax max df, stop if reached
#' @param verb verbosity
#' @export

glbin_lcd <- function(X,
                      y,
                      offset,
                      index,
                      eps = 1e-4,
                      lambda,
                      lambda.min = 0.01,
                      nlambda = 50,
                      dfmax = ncol(X)+1,
                      verb=1,
                      std=TRUE,
                      alpha = 1,
                      maxit = 1000,
                      stability_eps = 1e-5
) {
  p <- ncol(X)
  n <- length(y)

  # the grouping
  if(missing(index)) {
    index <- c(0, 1:(p-1)) # assuming first column is intercept
  }
  else{
    index[is.na(index)] <- 0
  }

  # two categories: penalised and unpenalised.
  un_i <- which(index==0)
  pen_i <- which(index!=0)
  # check for global intercept
  int_i <- all(X[,1]==1)
  # drop it from
  if(int_i) un_i <- setdiff(un_i,1)
  # penalised groups
  gr_i <- split((1:p)[pen_i], index[pen_i])
  # sqrt df
  sdfg <- sqrt( sapply(gr_i, length) )
  # groups
  ngroups <- length(gr_i)

  # standardize X, not unpenalised
  if(std){
    std_i <- 1 #un_i #un_std_i
    A <- X[,-std_i]
    A <- scale(A)
    X[,-std_i] <- A
    center <- attr(A, "scaled:center")
    scale <- attr(A, "scaled:scale")
    # if some are constant scale
    scale[is.na(scale) | scale == 0] <- 1
  }

  # offset
  if(missing(offset)){
    offset <- rep(0, n)
  }
  #
  # the penalty vector lambda
  if(missing(lambda)){
    # The maximal lambda.
    # check the intercept issue
    X0 <- if(int_i) X[,-1] else X
    index0 <- if(int_i) index[-1] else index
    lambda_vec <- makeLambda(nlambda = nlambda, lambda.min = lambda.min,
                         X=X0, y=y, offset =offset, index = index0)
  }
  else{
    lambda_vec <- sort(lambda, decreasing = TRUE)
    nlambda <- length(lambda_vec)
  }
  # soft threshold
  S <- function(a,b){
    if(a > b) a-b
    else if(a < -b) a+b
    else 0
  }
  #
  cat2 <- if(verb) cat else function(...) NULL

  # Use a quadratic approximation and IRLS

  # initials etc.
  beta <- rep(0, p)
  ybar <- mean(y)
  beta[1] <- log(ybar/(1-ybar))# - mean(offset)
  lik <- df <- NULL
  beta_res <- NULL
  beta_star <- beta
  it1 <- 0
  tX <- t(X)
  # lets go:
  for(lambda in lambda_vec){
    # IRLS
    loop_Q <- TRUE
    it2 <- 0
    odiff <- Inf
    while(loop_Q){
      df1 <- 1*int_i + length(un_i)
      beta1 <- beta
      # quadr. approx.
      #browser()
      eta <- c(crossprod(tX , beta)) + offset
      mu <- 1/(1 + exp(-eta))
      w <- .25 #mu * (1-mu)

      # working residual
      r <- (y - mu)/w
      # update intercept
      if(int_i) {
        shift <- sum(r)/n
        beta[1] <- beta[1] + shift
        r <- r - shift
        df1 <- df1 + 1
      }
      # update unpenalised
      shifts <- 0
      for(i in un_i){
        xi <- X[,i]
        shift <- sum(xi * w * r)/sum(xi * xi * w) # unscaled here
        #shifts <- shifts + shift * xi
        beta[i] <- shift + beta[i]
        beta_star[i] <- beta[i]
        r <- r - shift * xi
      }

      # the parameter updates in groups
      # for(j in 1:ngroups){
      #   ig <- gr_i[[j]]
      #   # check threshold
      #   #inner <- crossprod(X[,ig], v + crossprod(t(w*X[,ig]), beta[ig]))
      #   z <- crossprod(X[,ig], r)/n + beta[ig]
      #   z_norm <- sqrt( sum( z * z ) )
      #   lg <- lambda * sdfg[j]
      #   len <- S(z_norm * 0.25, lg) / 0.25
      #   # if(len == 0 & beta[ig[1]]==0){
      #   #   #beta[ig] <- 0
      #   # }
      #   # else{
      #   if(len !=0 | beta[ig[1]] != 0){
      #     beta[ig] <- len * z / z_norm
      #     #bnorm <- sqrt( sum(beta[ig] * beta[ig]) )
      #     # individual shrink
      #     #lt <- lg / (bnorm + 5e-3 + (1-alpha) * lg)
      #     #
      #     #for(k in ig){
      #       #for(k in 1:length(ig)){
      #       #xk <- X[,k]
      #       #sw <- sum(w * xk * xk)/n
      #       #newb <- sum(xk * v)/n + sw * beta[k]
      #       #beta[ig[k]] <- len * z[k] / z_norm
      #       #beta[k]<- newb / (sw + lt) # shrink
      #       #beta_star[ig[k]] <- 1 #newb / sw  # not shrink
      #     #}
      #   }
      #
      # }
      for(j in 1:ngroups){
        ig <- gr_i[[j]]
        z <- crossprod(X[,ig], r)/n + beta[ig]
        z_norm <- sqrt( sum( z * z ) )
        lg <- lambda * sdfg[j]
        len <- S(z_norm * 0.25, lg) / 0.25
        if(len !=0 | beta[ig[1]] != 0) {
          bn <- len * z / z_norm
          shift <- crossprod(rbind(tX[ig,]), beta[ig] - bn)
          beta[ig] <- bn
          r <- r + shift
        }
        if(len > 0) df1 <- df1 + length(ig) * len / z_norm
      }
#      eta <- c(crossprod(tX , beta)) + offset
#      mu <- 1/(1 + exp(-eta))
#      w <- 0.25
      # working residual
#      r <- (y - mu)/w
      # for(j in 1:ngroups){
      #   ig <- gr_i[[j]]
      #   z <- crossprod(X[,ig], r)/n + beta[ig]
      #   z_norm <- sqrt( sum( z * z ) )
      #   lg <- lambda * sdfg[j]
      #   len <- S(z_norm * 0.25, lg) / 0.25
      #   if(len !=0 | beta[ig[1]] != 0)
      #     beta[ig] <- len * z / z_norm
      # }

      diff <- max(abs(beta1-beta))
      #diff <- max(abs(beta1-beta)[-un_i])
      loop_Q <- diff > eps
      it2 <- it2 + 1
      if(it2 > maxit & loop_Q){
        cat2("[convergence problem", diff,"]\n")
        loop_Q <- FALSE
      }
      odiff <- diff
    }
    #
    df <- c(df, df1)
    # log-likelihood
    eta <- crossprod(tX, beta) + offset
    lik1 <- sum(y * eta - log(1+exp(eta)))
    lik <- c(lik, lik1)
    #
    #browser()
    beta_res <- cbind(beta_res, beta)
    it1 <- it1+1
    cat2("\r** outer loop [", it1, "/", nlambda,":", lambda,"] df", df1)
    if(df1 >= dfmax) {cat2(" -> dfmax reached\n"); break}
  }
  cat2("\n")

  # unstandardize
  if(std){
    beta_res[-std_i,] <- beta_res[-std_i,]/scale
    shift <- c( crossprod(center, beta_res[-std_i,, drop=FALSE]) )
    beta_res[std_i,] <- beta_res[std_i,] - shift
  }
  #
  rownames(beta_res) <- colnames(X)

  #
  z <- list(beta=beta_res, lambda=lambda_vec[1:it1], df = df, logLik = lik, aic=-2*lik + 2 * df)
  class(z) <- c("glbin", is(z))
  z
}



