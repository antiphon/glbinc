#' Group Penalised Logistic Regression
#'
#'
#'
#' @details
#' Several Lasso packages in CRAN. Needed one with 1) Group penalties 2) Support for sparse Matrices 3) Offset terms!
#'
#' None existed that fullfilled all three, so wrote this.
#'
#' We use local coordinate descent (LCD), as described in the reference.
#'
#' Main function is \code{\link{glbin_lcd_c}}, sparse version \code{\link{glbin_lcd_c_sparse}}.
#'
#'
#' @references
#' Patrick Breheny and Jian Huang, Penalized methods for bi-level variable selection, (2009) Statistics and Its Interface
#'
#' @docType package
"_PACKAGE"
