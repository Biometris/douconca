# unweighted least-squares (OLS) functions 
#' @noRd
#' @keywords internal
unweighted_lm_pq_fit <- function(Y, 
                                 X) {
  # multivariate multiple regression of Y on X
  # using qr decomposition
  # value Y_fit
  Y_fit <- qr.fitted(qr(X), Y)
  return(Y_fit)
}

#' @noRd
#' @keywords internal
unweighted_lm_Orthnorm_fit <- function(Y, 
                                       X = numeric(0)) {
  # multivariate multiple regression of Y on orthonormal X
  # value Y_residual
  beta <- t(X) %*% Y
  Yfit <- X %*% beta
  return(Yfit)
}

#' @noRd
#' @keywords internal
SVD <- function(Y) {
  svdY <- svd(Y)
  Ustar <- svdY$u
  id <- which(svdY$d > 1.e-6)
  return(Ustar[, id, drop = FALSE])
}