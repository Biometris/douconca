#' @title Performs (weighted) double constrained correspondence analysis (dc-CA)
#'
#' @description
#' \code{dc_CA} allows three versions of double constained correspondence analysis (dc-CA)
#' that differ in their preprocessing.
#' See the argument \code{preprocessY} below.
#' All versions use
#' the two-step algorithm of ter Braak et al. (2018). This algorithm
#' combines and extends community- (sample-) and species-level analyses.
#' The first step uses \code{\link[vegan]{cca}} (Oksanen et al. 2022)
#' to regress the transposed abundance data on to the traits
#' and (weighted) redundancy analysis
#' to regress the community-weighted means (CWMs) of the ortho-normalized traits,
#' obtained from the first step, on to the environmental predictors.
#' The sample total of the abundance data are used as weights.
#' The redundancy analysis is carried out using \code{vegan} \code{\link[vegan]{rda}}
#' if sites have equal weights (after division of the rows by their total) or,
#' in the general weighted case, using \code{\link{wrda}}.
#'
#' @param formulaEnv formula or one-sided formula for the rows (samples) with row predictors in \code{dataEnv}.
#' When two-sided, the left hand side of the formula is not used.
#' Specify row covariates (if any ) by adding \code{+ Condition(covariate-formula)}
#' to \code{formulaEnv} as in \code{\link[vegan]{rda}}. The \code{covariate-formula} should not contain a \code{~} (tilde).
#' Default: \code{~.}, i.e. all variables in \code{dataEnv} are predictor variables.
#' @param formulaTraits  formula or one-sided formula for the columns (species) with colum predictors in \code{dataTraits}.
#' When two-sided, the left hand side of the formula is not used. Specify column covariates (if any ) by  adding \code{+ Condition(covariate-formula)}
#' to \code{formulaTraits} as in \code{\link[vegan]{cca}}. The \code{covariate-formula} should not contain a \code{~} (tilde).
#' Default: \code{~.}, i.e. all variables in \code{dataTraits} are predictor traits.
#' @param response matrix, data frame of the abundance data (dimension \emph{n} x \emph{m}) or
#' object from \code{\link{fCWM_SNC}}.
#' Rownames of \code{response}, if any, are carried through. See Details for analyses starting from
#' community weighted means.
#' @param dataEnv matrix or data frame of the row predictors, with rows corresponding to those in \code{response}.
#' (dimension \emph{n} x \emph{p}).
#' @param dataTraits matrix or data frame of the column predictors,
#'  with rows corresponding to the columns in \code{response}.(dimension \emph{m} x \emph{q}).
#' @param divide.by.site.totals logical;
#'  default \code{TRUE} for closing the data by dividing the rows in the \code{response} by their total.
#' @param dc_CA_object  optional object from an earlier run of this function. Useful if the
#' same formula for the columns (\code{formulaTraits}), \code{dataTraits} and \code{response} are used
#' with a new formula for the rows. If set, the data of the previous run is used and the result of its first step
#' is taken for the new analysis.
#' The \code{data$Y} element in the object can be set to \code{NULL}, which is useful,
#' In the case of non-public abundance data and a wish for reproducibility of the analysis
#' use the function \code{fCWM_SNC} to create an object for use as \code{response} argument
#' in a new call to \code{dc_CA}. From this object, the analysis can be reproduced from CWMs
#' (and, for full analysis, SNCs) instead of from abundance data.
#'
#'
#' @param verbose logical for printing a simple summary (default: TRUE)
#
#' @details
#' Empty (all zero) rows and columns in \code{response} are removed from the \code{response} and the corresponding
#' rows from \code{dataEnv} and \code{dataTraits}. Subsequently, any columns with missing values
#' are removed from  \code{dataEnv} and \code{dataTraits}. It gives an error (object 'name_of_variable' not found),
#' if variables with missing entries are specified in \code{formulaEnv} and \code{formulaTraits}.
#'
#' The algorithm follows the two-step algorithm of ter Braak et al. (2018).
#' and consits of two steps. First, the transpose of the \code{response}
#' is regressed on to the traits (the column predictors) using \code{\link[vegan]{cca}}
#' with \code{formulaTraits}.
#' The column scores of this analysis (in scaling 1) are community weigthed means (CWM) of the
#' orthonormalized traits.
#' These are then regressed on the environmental (row) predictors using \code{\link{wrda}} with
#' with \code{formulaEnv}.
#'
#' A dc-CA can be carried out on, what statisticians call, the sufficient statistics of the methods
#' by another specification
#' of the \code{response}. In this case, \code{response} should be a list with as first element
#' community weighted means (CWMs with respect to the orthonormalized traits to result in a dc-CA)
#' and, optionally, further elements, for functions related to \code{dc_CA}.
#' If the CWMs are from standardized traits, the function does a CWM-RDA, as it is called in
#' Kleyer et al. 2012. Note that CWM-RDA is an ad-hoc method
#' which not have the mathematical rigour of dc-CA, \emph{e.g.} it does not optimize
#' the fourth-corner correlation between the trait and environment axes.
#' The function \code{\link{fCWM_SNC}} helps
#' to set the appropriate \code{response} in these non-standard applications of dc-CA.
#'
#' The statistics and scores in the example \code{dune_dcCA.r},
#'  have been checked against the results in Canoco 5.15 (ter Braak & Smilauer, 1918).
#'
#' In the current implementation, \code{formulaEnv} and \code{formulaTraits} should
#' contain variable names as is, \emph{i.e.} transformations of variables in the formulas gives
#' an error ('undefined columns selected') when the \code{\link{scores}} function is applied.
#'
#' @returns
#' A list of \code{class} \code{dccav}; that is a list with elements
#' \describe{
#' \item{CCAonTraits}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{cca}} analysis
#' of the transpose of the closed \code{response} using formula \code{formulaTraits}. }
#' \item{formalaTraits}{the argument \code{formulaTraits}. If the
#' formula was \code{~.}, it was changed to explicit trait names.}
#' \item{data}{a list of \code{Y} (response data after removing empty rows and columns and after closure)
#' and \code{dataEnv} and \code{dataTraits}.}
#' \item{weights}{a list of unit-sum weights of row and columns.
#' The names of the list are \code{c("row","columns")}, in that order.}
#' \item{Nobs}{number of sites (rows).}
#' \item{CWMs_orthonormal_traits}{Community weighted means w.r.t. orthonormalized traits.}
#' \item{RDAonEnv}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{rda}} analysis
#' of the column scores of the \code{cca}, which are the CWMs of orthonormalized traits,
#' using formula \code{formulaEnv}. }
#' \item{formalaEnv}{the argument \code{formulaEnv}. If the
#' formula was \code{~.}, it was changed to explicit environmental variable names. }
#' \item{eigenvalues}{the dc-CA eigenvalues (same as those of the \code{\link[vegan]{rda}} analysis)}
#' \item{inertia}{a one-column matrix with four inertias (weighted variances):
#' \itemize{
#' \item total: the total inertia.
#' \item conditionT: the inertia explained by the condition in \code{formulaTraits}
#' if present (neglecting row constraints).
#' \item traits_explain: the inertia explained by the traits (neglecting the row predictors and any
#' condition in \code{formulaTraits}).
#' This is the maximum that the row predictors could explain in dc-CA
#' (the sum of the following two items is thus less than this value).
#' \item conditionE: the trait-constrained inertia explained by the condition in \code{formulaEnv}.
#' \item constraintsTE: the trait-constrained inertia explained by the predictors (without the row covariates).
#' }
#'  }
#' }
#' If \code{verbose} is \code{TRUE} (or after \code{out <-print(out)} is invoked )
#' there are three more items (in this version).
#' \itemize{
#' \item \code{c_traits_normed}: mean, sd, VIF and (regression) coefficients of
#'  the traits that define the dc-CA trait axes (composite traits), and their optimistic t-ratio.
#' \item \code{c_env_normed}:  mean, sd, VIF and (regression) coefficients of the environmental variables that define the dc-CA axes
#'  in terms of the environmental variables (composite gradients), and their optimistic t-ratio.
#' \item \code{species_axes}: a list with four items
#'  \itemize{
#'  \item \code{species_scores}: a list with names \code{c("species_scores_unconstrained", "lc_traits_scores")} with the
#'  matrix with species niche centroids along the dc-CA axes (composite gradients) and
#'  the matrix with linear combinations of traits.
#'  \item \code{correlation}: a matrix with inter-set correlations of the traits with their SNCs.
#'  \item \code{b_se}: a matrix with (unstandardized) regression coefficients for traits and their optimistic standard errors.
#'  \item \code{R2_traits}: a vector with coefficient of determination (R2) of the SNCs on to the traits.
#'  The square-root thereof could be called the species-trait correlation in analogy with
#'  the species-environment correlation in CCA.
#'  }
#'  \item \code{sites_axes}: a list with four items
#'  \itemize{
#'  \item \code{site_scores}: a list with names \code{c("site_scores_unconstrained", "lc_env_scores")} with the
#'  matrix with community weighted means (CWMs) along the dc-CA axes (composite gradients) and
#'  the matrix with linear combinations of environmental variables.
#'  \item \code{correlation}: a matrix with inter-set correlations of the environmental variables with their CWMs.
#'  \item \code{b_se}: a matrix with (unstandardized) regression coefficients for environmental
#'  variables and their optimistic standard errors.
#'  \item \code{R2_env}: a vector with coefficient of determination (R2) of the CWMs on to the environmental variables.
#'  The square-root thereof has been called the species-environmental correlation in CCA.
#'  }
#'
#' }
#' All scores in the \code{dccav} object
#' are in scaling \code{"sites"} (1): the scaling with \emph{Focus on Case distances} .

#'
#' @references
#' Kleyer, M., Dray, S., Bello, F., Lepš, J., Pakeman, R.J., Strauss, B., Thuiller,
#' W. & Lavorel, S. (2012) Assessing species and community functional responses to
#' environmental gradients: which multivariate methods?
#' Journal of Vegetation Science, 23, 805-821.
#' http://dx.doi.org/10.1111/j.1654-1103.2012.01402.x
#'
#' ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' https://doi.org/10.1007/s10651-017-0395-x or
#' http://rdcu.be/ETPh
#'
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-4.
#' http://CRAN.R-project.org/package=vegan.
#'
#' @seealso \code{\link{plot_dcCA}}, \code{\link{scores.dcca}}, \code{\link{print.dcca}} and \code{\link{anova.dcca}}
#' @example demo/dune_dcCA.R
#' @export

dc_CA <- function(formulaEnv = NULL, formulaTraits = NULL,
        response =NULL, dataEnv=NULL, dataTraits= NULL, divide.by.site.totals = TRUE, dc_CA_object  = NULL, verbose = TRUE) {
  # response matrix or data frame, dataEnv and dataTraits data frames in which formualaE and formulaT are evaluated
  #dc_CA_object = result (value) of a previous run, can be used to save computing time for
  # runs that modify the formula for samples (step2: RDAonEnv) only
  # The step1 (CCAonTraits and the data and formulaTraits) are taken from dc_CA_object into the new result.
  # If set, formulaTraits, response, dataEnv, dataTraits are not used at all and have no efffect on the result
  call <- match.call()
  if (!is.null(response)) {
    if (is.list(response) && is.matrix(response[[1]])) {
        # response is a list of CWMs_orthonormal_traits and a weights list
        if (any(is.na(response[[1]])))stop("The CWMs should not have missing entries")
        # create a sufficient dc_CA_object object
       dc_CA_object <- response
    }
  }
  if (is.null(dc_CA_object)){
    #  check and amend: make sure there are no empty rows or columns -----------------------------------------------------------------------
    if (any(is.na(response)))stop("The response should not have missing entries")
    if (any(response <0)) stop("The response should not have negative values")
    if (is.null(dataTraits)) stop("dataTraits must be specified in dc_CA")
    if (!is.matrix(response)) response <- as.matrix(response)

    id0 <-1
    while(length(id0)){
      TotR <- rowSums(response)
      id0 <- which(TotR == 0)
      if (length(id0)){
        response <- response[-id0,]
        dataEnv  <- dataEnv[-id0,]
      }
      TotC <- colSums(response)
      id0 <- which(TotC == 0)
      if (length(id0)){
        response <- response[,-id0]
        dataTraits <- dataTraits[-id0, ]
      }
    }

    # delete columns with missing data
    id = rep(FALSE, ncol(dataEnv))
    for (ii  in seq_along(id)){
      id[ii] <- sum(is.na(dataEnv[,ii]))==0
      if (!id[[ii]]) warning(
        paste("variable", names(dataEnv)[ii], "has missing values and is deleted from the environmental data")
      )
    }
    dataEnv <- dataEnv[, id]

    id = rep(FALSE, ncol(dataTraits))
    for (ii  in seq_along(id)){
      id[ii] <- sum(is.na(dataTraits[,ii]))==0
      if (!id[[ii]]) warning(
        paste("variable", names(dataTraits)[ii], "has missing values and is deleted from trait data")
      )
    }
    dataTraits <- dataTraits[, id]

    dataEnv <- as.data.frame(lapply(dataEnv, function(x){if (is.character(x)) x<- as.factor(x) else x; return(x) } ))
    dataTraits <- as.data.frame(lapply(dataTraits, function(x){if (is.character(x)) x<- as.factor(x) else x; return(x) } ))

    rownames(dataEnv) <- rownames(response)
    rownames(dataTraits) <- colnames(response)


    # end of check ----------------------------------------------------------------------
    if (divide.by.site.totals){
      response <- sweep(response, 1, STATS = rowSums(response), FUN = '/')
    }
    TotR <- rowSums(response)
    TotC <- colSums(response)
    tY <- t(response)
    if (is.null(formulaTraits)) formulaTraits <- ~.
    if (is.null(formulaEnv)) formulaEnv <- ~.
    formulaTraits <- change_reponse(formulaTraits, "tY", dataTraits)
    environment(formulaTraits)<- environment()
    step1 <-vegan::cca(formulaTraits, data = dataTraits)
    data= list(Y = response, dataEnv = dataEnv, dataTraits = dataTraits)
    n <- nrow(data$Y)
    CWMs_orthonormal_traits <- vegan::scores(step1, display= "species",
                                             scaling = "species",
                                             choices = seq_len(Rank_mod(step1)) ) * sqrt((n-1)/n)
    if (rownames(CWMs_orthonormal_traits)[1]=="col1") rownames(CWMs_orthonormal_traits) <- paste("Sam", seq_len((nrow(dataEnv))),sep="")
    out1 <- list(CCAonTraits = step1,
                 formulaTraits= formulaTraits,
                 data = list(Y = response, dataEnv = dataEnv, dataTraits = dataTraits),
                 call = call,
                 weights = list(rows = TotR/sum(TotR), columns = TotC/sum(TotC)),
                 Nobs = n,
                 CWMs_orthonormal_traits = CWMs_orthonormal_traits
    )
  } else {
    #step1 <- dc_CA_object$CCAonTraits
    if ("dcca" %in% class(dc_CA_object))
      out1 <- dc_CA_object[c("CCAonTraits", "formulaTraits","data","call","weights","Nobs","CWMs_orthonormal_traits")]
    else if (is.list(response) && is.matrix(response[[1]])) {
      out1 <- dc_CA_object
      out1$Nobs <- nrow(out1$CWM)
      if(is.null(dataEnv)) dataEnv <- out1$data$dataEnv else dataEnv <- as.data.frame(lapply(dataEnv, function(x){if (is.character(x)) x<- as.factor(x) else x; return(x) } ))
      if (is.null(dataTraits)) {
        if (!is.null(out1$dataTraits)) dataTraits <- out1$dataTraits else if (!is.null(out1$data$dataTraits))
                               dataTraits <- out1$data$dataTraits else warning(" No trait data supplied to the dc_CA function.")
      } else { warning(paste(" With CWM as first element in response in dc_CA, the trait data",
         "used to obtain the CWMs are best supplied as response$dataTraits or response$data$dataTraits.",
         "Use the default dataTraits argument, which is NULL."))
          dataTraits <- as.data.frame(lapply(dataTraits, function(x){if (is.character(x)) x<- as.factor(x) else x; return(x) } ))
      }
      CWM2ortho <-  f2_orth(out1$CWM,out1$formulaTraits,dataTraits,out1$weights$columns)
      out1$CWMs_orthonormal_traits <- CWM2ortho$CWMs_orthonormal_traits * sqrt((out1$Nobs-1)/(out1$Nobs))
     # out1$traits_explain <- sum(out1$CWMs_orthonormal_traits^2*out1$weights$rows)*(out1$Nobs)/(out1$Nobs-1)
      if (is.null(formulaEnv)){
        if(!is.null(out1$formulaEnv)) formulaEnv <- out1$formulaEnv else formulaEnv<- ~.
      }
      #print(formulaEnv)
      out1$formulaEnv <- NULL
       if (!is.null(out1$SNC)&& !is.null(out1$weights$rows)){
        # print(out1$formulaEnv)
        #  print(str(out1$dataEn))
        SNC2ortho <-  f2_orth(out1$SNC,formulaEnv,dataEnv,out1$weights$rows)
        out1$SNCs_orthonormal_env <- SNC2ortho$CWMs_orthonormal_traits
      }
       out1$data <- list(dataEnv = dataEnv, dataTraits = dataTraits)

      out1$data <- list(CWM = out1$CWM, dataEnv = dataEnv, dataTraits = dataTraits)
      out1$CWM <- NULL
      out1$CWM2CWM_ortho <- CWM2ortho$CWM2CWM_ortho

    } else{
       stop(paste(" the class of dc_CA_object should be dcca, whereas it is now:",
                      class(dc_CA_object)[1],class(dc_CA_object)[2]))
    }
  }

  formulaEnv <- change_reponse(formulaEnv, "out1$CWMs_orthonormal_traits", data = out1$data$dataEnv)
  environment(formulaEnv)<- environment()
  if (diff(range(out1$weights$rows)) < 1.0e-4/out1$Nobs ){
    step2 <- vegan::rda(formulaEnv, data = dataEnv)
    eigenvalues <-  vegan::eigenvals(step2, model = "constrained")
  } else{
    step2 <- wrda(formulaEnv, response =out1$CWMs_orthonormal_traits, weights = out1$weights$rows ,data = dataEnv)
    eigenvalues <-  step2$CCA$eig
  }

  out <- c(out1, list(RDAonEnv = step2,
                      formulaEnv = formulaEnv,
                      eigenvalues =  eigenvalues
                      )
  )

  out$c_traits_normed0 <- try(f_canonical_coef_traits2(out))

  inertia <- try(f_inertia(out))
  if("try-error"%in% class(inertia)) {warning("could not obtain inertia's"); print(inertia)}

  out$inertia <- inertia
  if ("rda"%in% class(out$RDAonEnv))  class(out) <- c("dccav", "dcca" ,"list") else
                                      class(out) <- c("dcca", "list")

  if (verbose) {
    out<-print(out)
  }




  return(out)
}