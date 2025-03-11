#' @title Performs a weighted redundancy analysis
#'
#' @description
#' \code{wrda} is formula-based implementation of weighted redundancy analysis.
#' 
#' @inheritParams cca0
#'
#' @param response matrix or data frame of the abundance data (dimension 
#' \emph{n} x \emph{m}). Rownames of \code{response}, if any, are carried 
#' through.
#' @param weights row weights (a vector). If not specified unit weights are 
#' used.
#
#' @details
#' The algorithm is a modified version of published R-code for weighted 
#' redundancy analysis (ter Braak, 2022).
#'
#' Compared to  \code{\link[vegan]{rda}}, \code{wrda} does not have residual 
#' axes, \emph{i.e.} no SVD or PCA of the residuals is performed.
#'
#' @returns
#' All scores in the \code{wrda} object are in scaling \code{"sites"} (1): 
#' the scaling with \emph{Focus on Case distances}.
#'
#' @references
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-4.
#' \url{https://CRAN.R-project.org/package=vegan}.
#'
#' @seealso \code{\link{scores.wrda}}, \code{\link{anova.wrda}},
#' \code{\link{print.wrda}}
#' 
#' @example demo/dune_wrda.R
#' 
#' @export
wrda <- function(formula, 
                 response, 
                 data, 
                 weights = rep(1 / nrow(data), nrow(data)),
                 traceonly = FALSE,
                 cca_object = NULL,
                 object4QR = NULL) {
  call <- match.call()
  if (!is.null(cca_object)) {
    # check whether eY works with and without trait covariates
    Yw <- cca_object$Ybar
    transp <- all(dim(Yw) !=dim(response))
    if (nrow(Yw) == ncol(Yw)) {
      if (rownames(Yw)[1] == colnames(response)[1]) transp <- TRUE else
        if (rownames(Yw)[1] == rownames(response)[1]) transp <- FALSE else {
          warning("cca_object not usable.\n")
          cca_object <- NULL
        }
    }
  }
  if (is.null(cca_object)) {
    Yn <- as.matrix(response) 
    Wn <- weights / sum(weights)
    sWn <- sqrt(Wn)
    # transform to the unweighted case
    Yw <- as.matrix(response) * sWn
    # center
    Yw <- unweighted_lm_Orthnorm(Yw, matrix(sWn)) * sqrt(nrow(Yn) / (nrow(Yn) - 1))
    total_variance <- sum(Yw ^ 2) #total_inertia
  } else {  # use cca_object
    if (transp) {
      Yw <- t(Yw) 
      cca_object$weights <- rev(cca_object$weights)
    }
    K <- cca_object$weights[[1]]
    Wn <- cca_object$weights[[2]]
    sWn <-sqrt(Wn)
    sumY <- cca_object$sumY
    total_variance <- cca_object$tot.chi
  } # cca_object
  msqr <- msdvif(formula, data = data, weights = Wn, XZ = TRUE,
                 object4QR = object4QR)
  eY <- qr.resid(msqr$qrZ, Yw)
  Yfit_X <- qr.fitted(msqr$qrXZ, eY)
  ssY_gZ <- sum(eY ^ 2) 
  ssY_XgZ <- sum(Yfit_X ^ 2)
  if (traceonly) {
    return(c(Condition_inertia =total_variance-ssY_gZ, Fit_inertia = ssY_XgZ))
  }
  svd_Yfit_X <- SVDfull(Yfit_X)
  biplot <- NULL
  eig <- svd_Yfit_X$d ^ 2
  names(eig) <- paste0("wRDA", seq_along(eig))
  CCA <- with(svd_Yfit_X, list(eig = eig, poseig = eig, u = u, v = v, 
                               rank = rank, qrank = msqr$qrXZ$rank,
                               tot.chi = ssY_XgZ, QR = msqr$qrXZ, 
                               biplot = biplot, envcentre = NULL, 
                               centroids = NULL))
  rank_CA <- min(nrow(response)-1, ncol(response)) 
  eig_CA <- rep(NA, rank_CA)
  names(eig_CA) <- paste0("wPCA", seq_len(rank_CA))
  if (ncol(msqr$Zw) == 1) {
    pCCA <- NULL 
  } else {
    pCCA <- list(rank = min(ncol(Yw), msqr$qrZ$rank), 
                 tot.chi = total_variance - ssY_gZ,
                 QR = msqr$qrZ, envcentre = NULL)
  }
  if (length(svd_Yfit_X$d) == 1) {
    diagd <- matrix(svd_Yfit_X$d)
  } else {
    diagd <- diag(svd_Yfit_X$d)
  }
  # need to be orthogonalized w.r.t Z
  site_axes <- list(
    site_scores = list(
      site_scores_unconstrained = qr.resid(msqr$qrZ, eY %*% svd_Yfit_X$v) / sWn,
      lc_env_scores = (svd_Yfit_X$u %*% diagd) /sWn
    )
  )
  species_axes <- list(species_scores = 
                         list(species_scores_unconstrained = svd_Yfit_X$v))
  object <- list(call = call, method = "wrda", tot.chi = total_variance,
                 formula = formula, site_axes = site_axes, 
                 species_axes = species_axes, Nobs = nrow(Yw), eigenvalues = eig,
                 weights = list(columns = rep(1 / ncol(eY), ncol(eY)), rows = Wn),
                 data = data, Ybar = Yw, pCCA = pCCA, CCA = CCA, 
                 CA = list(tot.chi = ssY_gZ - ssY_XgZ, rank = rank_CA, eig = eig_CA),
                 inertia = "weighted variance"
  )
  class(object) <- "wrda"
  return(object)
}
