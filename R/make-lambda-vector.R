#' Create the lambda vector
#'
#' @param nlambda default 100
#' @param lambda.min default 0.001
#' @param ... passed on to max_penalty to compute the max lambda if lambda.max is missing
#' @param lambda.max no default
#' @param logscale use log-scale equidistant? Default: true
#' @export

makeLambda <- function(nlambda = 100, lambda.min = 0.001, ..., lambda.max, log.scale=TRUE){
  if(missing(lambda.max)) lambda.max <- max_penalty(...)

  if(log.scale){
    f1 <- log
    f2 <- exp
  }
  else{
    f1 <- f2 <- identity
  }
  if(lambda.min == 0) {
    lambda <- c(f2( seq( f1(lambda.max), f1(0.001 * lambda.max), length=nlambda-1) ), 0)
  }
  else{
    lambda <- f2( seq( f1(lambda.max), f1(lambda.min * lambda.max), length=nlambda) )
  }
  lambda
}
