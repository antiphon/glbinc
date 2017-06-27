#' Block coordinate gradient descent for logistic group lasso (c-version)
#'
#' Like in Breheny and Huang 2009
#'
#' @param X covariate/design matrix. Global intercept should not be here.
#' @param y response 0/1 vector
#' @param index grouping. set 0 or NA for unpenalised terms. if not given, assumes X[,1] is intercept
#' @param offset offset terms
#' @param eps convergence criterion
#' @param lambda vector of lambdas
#' @param lambda.min fraction of max lambda to go down to
#' @param nlambda number of penalties
#' @param lambda.log create log-equidistant lambda vec. default: TRUE
#' @param dfmax max df, stop if reached
#' @param add_intercept should the global intercept be added. default: TRUE
#' @param AIC_stop default 0. After aic has increased this many steps, halt. 0: go to the end of lambda. If used, should be more than 1.
#' @param verb verbosity
#' @param gc_force FALSE. force garbage collection before going c-side? Might free memory for large X.
#'
#' @details
#' Like in Breheny and Huang 2009 except this version includes offset terms.
#'
#' @useDynLib glbinc
#' @export

glbin_lcd_c <- function(X,
                        y,
                        offset,
                        index,
                        eps = 1e-4,
                        lambda,
                        lambda.min = 0.001,
                        nlambda = 100,
                        lambda.log = TRUE,
                        dfmax = ncol(X)+1,
                        verb=1,
                        add_intercept = TRUE,
                        std = TRUE,
                        alpha = 1,
                        maxit = 1000,
                        AIC_stop = 0,
                        gc_force = FALSE
) {

  # check for global intercept
  int_i <- all(X[,1]==1)
  # drop it from X if exists
  if(int_i) {
    X <- X[,-1]
    index <- index[-1]
  }
  n <- length(y)
  p <- ncol(X)
  # the grouping
  if(missing(index)) {
    index <- 1:p
  }
  else{
    index[is.na(index)] <- 0 # for compatibility with grplasso
  }
  # two categories: penalised and unpenalised.
  un_i <- which(index==0)
  pen_i <- which(index!=0)
  # penalised groups
  gr_i <- table(index)
  # sqrt df
  group_weights <- sqrt( gr_i )
  # standardize X. this makes X dense
  if(std){
    X <- scale(X)
    center <- attr(X, "scaled:center")
    scale <- attr(X, "scaled:scale")
    # if some are constant scale
    scale[is.na(scale) | scale == 0] <- 1
  }
  if(gc_force) gc()
  # offset
  if(missing(offset)){
    offset <- rep(0, n)
  }
  #
  # the penalty vector lambda
  if(missing(lambda)){
    lambda_vec <- makeLambda(nlambda = nlambda, lambda.min = lambda.min, log.scale = lambda.log,
                         X=X, y=y, offset =offset, index = index)
  }
  else{
    lambda_vec <- sort(lambda, decreasing = TRUE)
    nlambda <- length(lambda_vec)
  }
  # lets go
  K <- table(index)
  G0 <- as.integer( if(min(index)==0) sum(index==0) else 0 )
  G1 <- as.integer( if(min(index)==0) cumsum(K) else c(0, cumsum(K)) )
  if(G0) lambda_vec[1] <- lambda_vec[1] + 1e-5

  group_weight <- sqrt( table(index[index!=0]) )

  if(gc_force) gc()
  #

  out <- glbin_lcd_cpp(X, y, offset,
                       G0, # unpenalised 0-G0
                       G1, # penalised G0-G1
                       group_weight,
                       lambda_vec,
                       add_intercept,
                       alpha,
                       verb,
                       eps,
                       dfmax,
                       maxit,
                       AIC_stop
  )

  beta_res <- out$beta
  rownames(beta_res) <- colnames(X)
  df <- out$df
  #
  lik <- out$lik
  aic <- out$AIC
  # did we get max df
  its <- sum(lik!=0)
  ok <- 1:its
  if(its < nlambda){
    # set uncomputed things to NA
    beta_res[,-ok] <- NA
    lik[-ok] <- NA
    df[-ok] <- NA
    aic[-ok] <- NA
  }

  # unstandardize
  std_i <- 1:p
  if(add_intercept){
    beta_res <- rbind("(Intercept)"=out$beta0, beta_res)
    std_i <- 1+std_i
  }

  if(std){
    beta_res[std_i, ok] <- beta_res[std_i, ok, drop=FALSE]/scale
    shift <- c( crossprod(center, beta_res[std_i, ok, drop=FALSE]) )
    beta_res[-std_i, ok] <- beta_res[-std_i, ok] - shift
  }
  #
  #
  index_out <- index
  if(add_intercept) index_out <- c(0, index_out)

  z <- list(beta=beta_res, lambda=lambda_vec, df = df, logLik = lik, aic=aic, index = index_out)
  class(z) <- c("glbin", is(z))
  z
}



#' Block coordinate gradient descent for logistic group lasso (c-version) separate std
#'
#' Like in Breheny and Huang 2009
#'
#' @param X covariate/design matrix. Global intercept should not be here.
#' @param y response 0/1 vector
#' @param index grouping. set 0 or NA for unpenalised terms. if not given, assumes X[,1] is intercept
#' @param offset offset terms
#' @param eps convergence criterion
#' @param lambda vector of lambdas
#' @param lambda.min fraction of max lambda to go down to
#' @param nlambda number of penalties
#' @param dfmax max df, stop if reached
#' @param add_intercept should the global intercept be added. default: TRUE
#' @param AIC_stop default 0. After aic has increased this many steps, halt. 0: go to the end of lambda. If used, should be more than 1.
#' @param verb verbosity
#' @param gc_force FALSE. force garbage collection before going c-side? Might free memory for large X.
#'
#' @useDynLib glbinc
#' @export

glbin_lcd_c_std <- function(X,
                        y,
                        offset,
                        index,
                        eps = 1e-4,
                        lambda,
                        lambda.min = 0.001,
                        nlambda = 100,
                        dfmax = ncol(X)+1,
                        verb=1,
                        add_intercept = TRUE,
                        std = TRUE,
                        alpha = 1,
                        maxit = 1000,
                        AIC_stop = 0,
                        gc_force = FALSE
) {

  # check for global intercept
  int_i <- all(X[,1]==1)
  # drop it from X if exists
  if(int_i) {
    X <- X[,-1]
    index <- index[-1]
  }
  n <- length(y)
  p <- ncol(X)
  # the grouping
  if(missing(index)) {
    index <- 1:p
  }
  else{
    index[is.na(index)] <- 0 # for compatibility
  }
  # two categories: penalised and unpenalised.
  un_i <- which(index==0)
  pen_i <- which(index!=0)
  # penalised groups
  gr_i <- table(index)
  # sqrt df
  group_weights <- sqrt( gr_i )
  # standardize X
  if(std){
    std_i <- 1:p
    center <- colMeans(X[,std])
    scale <- apply(X[,std], 2, sd)
    # if some are constant scale
    scale[is.na(scale) | scale == 0] <- 1
  }
  if(gc_force) gc()
  # offset
  if(missing(offset)){
    offset <- rep(0, n)
  }
  #
  # the penalty vector lambda
  if(missing(lambda)){
    lambda_vec <- makeLambda(nlambda = nlambda, lambda.min = lambda.min,
                             X=X, y=y, offset =offset, index = index)
  }
  else{
    lambda_vec <- sort(lambda, decreasing = TRUE)
    nlambda <- length(lambda_vec)
  }
  # lets go
  K <- table(index)
  G0 <- as.integer( if(min(index)==0) sum(index==0) else 0 )
  G1 <- as.integer( if(min(index)==0) cumsum(K) else c(0, cumsum(K)) )
  if(G0) lambda_vec[1] <- lambda_vec[1] + 1e-5

  group_weight <- sqrt( table(index[index!=0]) )

  if(gc_force) gc()
  #

  out <- glbin_lcd_cpp(X, y, offset,
                       G0, # unpenalised 0-G0
                       G1, # penalised G0-G1
                       group_weight,
                       lambda_vec,
                       add_intercept,
                       alpha,
                       verb,
                       eps,
                       dfmax,
                       maxit,
                       AIC_stop
  )

  beta_res <- out$beta
  rownames(beta_res) <- colnames(X)
  df <- out$df
  #
  lik <- out$lik
  aic <- out$AIC
  # did we get max df
  its <- sum(lik!=0)
  ok <- 1:its
  if(its < nlambda){
    # set uncomputed things to NA
    beta_res[,-ok] <- NA
    lik[-ok] <- NA
    df[-ok] <- NA
    aic[-ok] <- NA
  }


  # unstandardize
  if(add_intercept){
    beta_res <- rbind("(Intercept)"=out$beta0, beta_res)
    std_i <- 1+std_i
  }

  if(std){
    beta_res[std_i,ok] <- beta_res[std_i,ok]/scale
    shift <- c( crossprod(center, beta_res[std_i,ok, drop=FALSE]) )
    beta_res[-std_i,ok] <- beta_res[-std_i,ok] - shift
  }
  #
  #
  index_out <- index
  if(add_intercept) index_out <- c(0, index_out)
  z <- list(beta=beta_res, lambda=lambda_vec, df = df, logLik = lik, aic=aic, index = index_out)
  class(z) <- c("glbin", is(z))
  z
}



