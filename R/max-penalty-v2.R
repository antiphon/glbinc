#' Compute the maximum penalty, possibly sparse
#'
#'
#'
#' @param X design matrix without intercept
#' @param y response 0/1
#' @param index grouping
#' @param offset offset
#' @param center the column means for standardization
#' @param scale sd vector of columns for standardization.
#' @param alpha lasso vs ridge weight
#' @param add_intercept TRUE
#' @param use_glm use the glm-function?
#' @param ... passed on to c-methods
#' @details
#' See glbin for details of the parameters
#' @import Matrix
#' @export

max_penalty <- function(X, y, index, offset, add_intercept = TRUE,
                        center, scale, alpha = 1, eps=1e-3, maxiter=1000, verb=0, use_glm = FALSE, ...){
  p <- ncol(X)
  un_i <- which(index==0)
  pen_i <- which(index!=0)
  is_sparse <- is(X, "sparseMatrix")
  #
  if((missing(center)| missing(scale))){
    center <- colMeans(X)
    scale <- apply(X, 2, sd)
  }


  if(length(un_i)){
    Xun <- cbind( X[, un_i] )

    if(use_glm){
      Xun <- t((t(as.matrix(Xun))-center[un_i])/scale[un_i])
      if(add_intercept) Xun <- cbind(1, Xun)
      fit0 <- glm(y~ -1 + Xun, offset = offset, family = "binomial", control = glm.control(eps, maxiter, verb), ... )
      fit <- list(r=fit0$residuals, w=fit0$weights)
    }
    else{ # own
      is_sparse_un <- is(Xun, "sparseMatrix")
      callf <- if(is_sparse_un) glm_binom_std_sparse_c else glm_binom_std_c
      fit <- callf(Xun, y, offset,
                   add_intercept = add_intercept, eps = eps,
                   center = center[un_i], scale = scale[un_i],
                   maxiter = maxiter, verb = verb, ...
      )
    }
  }
  else{ # no covariates. no need to go sparse
    if(!add_intercept) stop("model has no terms unpenalised terms or intercept.")
    x0 <- matrix(0, nrow=1, ncol=0)
    fit <- glm_binom_std_c(x0, y, offset,
                           add_intercept = add_intercept, eps = eps,
                           center = 0, scale = 1,
                           maxiter = maxiter, verb = verb, ...
    )
  }
  w <- fit$w

  if (max(w) < 1e-4) stop("Unpenalized portion of model is already saturated; exiting...")
  r <- fit$r * w
  n <- length(y)
  # penalised groups
  gr_i <- split((1:p)[pen_i], index[pen_i])
  # sqrt df
  sdfg <- sqrt( sapply(gr_i, length) )
  # groups
  ngroups <- length(gr_i)
  # the gradients per group
  cp <- if(is_sparse) Matrix::crossprod else base::crossprod

  gest <- sapply(1:ngroups, function(g){
    ig <- gr_i[[g]]
    Xg <- t((t(X[,ig]) - center[ig])/scale[ig] )
    z <- cp(Xg, r)
    sqrt(sum(z^2))/sdfg[g]/n
  })

  max(gest)/alpha
}

