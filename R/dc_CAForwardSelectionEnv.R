#' Forward selection of environmental variables using dc-CA  with a
#' fixed trait formula
#' 
#' @param mod initial dA-CAmodel with at least on environmental variables, 
#' and the a trait formula to be used throughout the selection description.
#' @param consider character vector of names in \code{mod$data$dataEnv}
#' to consider of addition; default \code{NULL} for all variables in 
#' \code{mod$data$dataEnv}.
#' @param permutations a list of control values for the permutations as 
#' returned by the function \code{\link[permute]{how}}, or the number of 
#' permutations required (default 999), or a permutation matrix where each 
#' row gives the permuted indices. 
#' If \code{test = 2}, a two-vector with the number for the
#' species-level test first and that for the sites-level second), or
#' a list of two permutation matrices (again, species first, sites second)
#' where each row gives the permuted indices.
#' @param n_axes number of dc-CA eigenvalues to select upon.
#' The sum of \code{n_axes} eigenvalues is taken as criterion.
#' Default \code{"full"} for selection without dimension reduction to 
#' \code{n_axes}.
#' If \code{n_axes =1}, selection is on the first eigenvalue of the dc-CA 
#' for selection of variables that form an optimal one-dimensional model.
#' @param test integer value for statistical testing; default: \code{2}. 
#' For \code{test = 0,1,2}, no test, site-based permutation test and site- 
#' & species-based permutation tests, respectively.
#' With \code{test = 1}, the test is likely liberal.
#' @param threshold_P  significance level, after adjustment for testing 
#' multiplicity, for addition of a variable to the model.
#' @param initial_env_model character specifying what should be inside 
#' Condition(). Default = \code{"1"}, \emph{i.e.} nothing.
#' @param permutations a list of control values for the permutations as 
#' returned by the function \code{\link[permute]{how}}, or the number of 
#' permutations required (default 999), or a permutation matrix where each 
#' row gives the permuted indices. 
#' If \code{test = 2}, a two-vector with the number for the
#' species-level test first and that for the sites-level second), or
#' a list of two permutation matrices (again, species first, sites second)
#' where each row gives the permuted indices.
#' @param PvalAdjustMethod method for correction for multiple testing 
#' in  \code{\link[stats]{p.adjust}}, default \code{"holm"}, which
#' is an improved version Bonferroni.
#' @param max_step maximal number of variables selected. 
#' @param verbose show progress, default: \code{TRUE}.
#' @details
#' 
#' The selection is on the basis of the additional fit (inertia) of a variable 
#' given the variables already in the model.
#' 
#' @returns list with three elements: \code{final...} with selected variables
#'  and \code{model_final} and \code{process} with account of the selection process 
#' @example demo/dune_forward_slct.R
#' 
#' @export
dc_CAForwardSelectionEnv <- function(mod,  
                                     consider = NULL, 
                                     permutations = 999,
                                     n_axes = "all",
                                     initial_env_model = "1",
                                     test = 2, 
                                     threshold_P = 0.10, 
                                     PvalAdjustMethod = "holm",
                                     max_step = 10, 
                                     verbose = TRUE) {
  # function for forward selection of environmental variables using dc-CA with a
  # fixed trait formula.
  divide <- diff(range(mod$weights$rows)) < 1.0e-4 / mod$Nobs
  veganperm <- divide && is.character(n_axes)
  formulaTraits <- mod$formulaTraits
  if (!length(consider)) consider <- names(mod$data$dataEnv)
  ## step 1
  considerk <- consider
  names(consider) <- consider
  considered <- initial_env_model
  fit_measures <- matrix(0, nrow = length(consider), ncol = 3)
  fit_measuresL <- list()
  attr(n_axes, which = "Names_initial_model") <- 
    colnames(model.matrix(make_formula(1, " ", initial_env_model),
                          data = mod$data$dataEnv))
  ### step 2 - end
  stepk <- 0
  while (fit_measures[1, 3] < threshold_P && 
         length(considerk) && stepk < max_step) {
    stepk <- stepk + 1
    fit_measures <- matrix(nrow = length(considerk), ncol = 2)
    colnames(fit_measures) <- c("variance", "pval1")
    rownames(fit_measures) <- considerk
    for (k in seq_along(considerk)) {
      out_FS <- try(dc_CA(formulaEnv = make_formula(n_axes, considered, 
                                                    considerk[k]),
                          dc_CA_object = mod,
                          verbose  = FALSE))
      crit <- if (inherits(out_FS, "try-error") || 
                  is.null(out_FS$eigenvalues)) 0 else 
                    selection_crit(out_FS,n_axes)
      if (test && crit > 0) {
        if (is.numeric(n_axes)){
          out_FS<- try(dc_CA(formulaEnv = make_formula("all", considered,  considerk[k]),
                             dc_CA_object = mod, verbose  = FALSE))
        }
        if (inherits(out_FS, "try-error") ||is.null(out_FS$eigenvalues)) pval <- NA else{
          if (test == 1) {
            an <- if (veganperm) anova(out_FS$RDAonEnv, permutations = permutations) else 
              anova_sites(out_FS, permutations = permutations, n_axes = n_axes)$table
            pval <- an$`Pr(>F)`[1]
          } else {
            an <- if (veganperm) anova(out_FS, permutations = permutations) else
              anova(out_FS, permutations = permutations, n_axes = n_axes)
            pval <- an$maxP$`Pr(>F)`[1]
          }
        }
      } else pval <- NA
      fit_measures[k,] <- c(variance=crit, pval1 = pval)
    }
    
    pvaladj <- if (test) stats::p.adjust(fit_measures[,"pval1"], 
                                         method = PvalAdjustMethod) else 0
    fit_measures <- cbind(fit_measures, pvaladj)
    id <- order(fit_measures[,"variance"], decreasing = TRUE)
    fit_measures <- fit_measures[id, , drop = FALSE]
    if (verbose) print(round(fit_measures, 5))
    fit_measuresL[[stepk]] <- fit_measures
    best_trait <- rownames(fit_measures)[1]
    considered <- c(considered, consider[best_trait])
    if (verbose) {
      print("considered") 
      print(considered)
    }
    considerk <-  rownames(fit_measures)[-1]
    fit_measuresL[[stepk]] <- fit_measures
  }
  names(fit_measuresL) <- considered[-1]
  # selected variables with last non-significant value
  tab <- t(sapply(fit_measuresL, function(x) {
    x[1, , drop = FALSE]
  }))
  colnames(tab) <- colnames(fit_measures)
  selected_vars <- if (nrow(tab) == 1) rownames(tab) else rownames(tab)[-nrow(tab)]
  formulaEnv2 <- make_formula(n_axes = "all", initial_env_model, 
                              paste(selected_vars, collapse = "+"))
  mod2 <- dc_CA(formulaEnv = formulaEnv2, 
                dc_CA_object = mod, verbose = FALSE)
  ret <- list(finalWithOneExtra = tab, model_final = mod2, 
              process = fit_measuresL)
  if (nrow(tab) == 1 && max_step > 1) {
    names(ret)[1] <- 
      if (length(consider) > 1) "final_Non_Significant" else "final" 
  }
  return(ret)
}

#' returns the selection_criterion.
#' 
#' @noRd
#' @keywords internal
selection_crit <- function(out_FS, 
                           n_axes = "full"){
  if (is.character(n_axes)) {
    crit <- out_FS$inertia["constraintsTE", "weighted variance"] 
  } else {
    eig <- if (inherits(out_FS$RDAonEnv, "rda"))
      vegan::eigenvals(out_FS$RDAonEnv, model = "constrained") else 
        sum(out_FS$eigenvalues)
    n_eig <- length(eig)
    crit <- sum(eig[seq_len(min(n_eig, n_axes))])
  }
  return(crit)
}

#' @param n_axes number of dc-CA eigenvalues. The sum of \code{n_axes} 
#' eigenvalues is taken  as criterion. Default \code{"full"} for selection 
#' without n_axes reduction.
#' If \code{n_axes =1}, selection is on the first eigenvalue of the dc-CA with 
#' all variables in the model.
#' 
#' @param considered variables currently in the model
#' @param considerkk variable_to_consider
#' 
#' @noRd
#' @keywords internal
make_formula <- function(n_axes, 
                         considered, 
                         considerkk) {
  # make formulas of all variables  if(is.numeric(n_axes )) else
  # a formula with Condition(variables currently in the model) + variable_to_consider.
  fconsidered <- paste(considered, collapse = "+")
  if (is.character(n_axes)) {
    formula_FS <-
      paste("~", considerkk, "+ Condition(", fconsidered , ")", collapse = " ")
  } else {
    formula_FS <- paste("~", fconsidered, " + ", considerkk , collapse = " ")
  }
  formula_FS <- as.formula(formula_FS)
  return(formula_FS)
}
