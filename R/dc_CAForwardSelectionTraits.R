#' Forward selection of traits using dc-CA  with a fixed environmental formula
#' 
#' @inheritParams dc_CAForwardSelectionEnv
#' 
#' @param mod initial dA-CAmodel with at least on one trait, and an
#' environmental formula to be used throughout the selection.
#' @param consider character vector of names in \code{mod$data$dataTraits} to
#' consider of addition. Default \code{NULL} for all variables in 
#' \code{mod$data$dataTraits}.
#' @param initial_env_model character specifying what should be inside 
#' Condition(). Default = "1" (nothing)
#' 
#' @details 
#' The selection is on the basis of the additional fit (inertia) of a variable 
#' given the variables already in the model. 
#' 
#' @returns list with three elements: \code{final...} with selected variables
#' and \code{model_final} and \code{process} with account of the selection 
#' process 
#' 
#' @example demo/dune_forward_slct.R
#' @export
dc_CAForwardSelectionTraits <- function(mod,  
                                        consider = NULL, 
                                        permutations = 999, 
                                        n_axes = "all",
                                        initial_trait_model = "1", 
                                        test = 2, 
                                        threshold_P = 0.10,
                                        PvalAdjustMethod = "holm",
                                        max_step = 10, 
                                        verbose = TRUE) {
  # function for forward selection of traits using dc-CA with a
  # fixed trait formula.
  divide <- diff(range(mod$weights$rows)) < 1.0e-4 / mod$Nobs
  veganperm <- divide & is.character(n_axes)
  formulaEnv <- mod$formulaEnv
  if (!length(consider)) consider <- names(mod$data$dataTraits)
  ### step 1
  considerk <- consider
  names(consider) <- consider
  considered <- initial_trait_model
  fit_measures <- matrix(0, nrow = length(consider), ncol = 3)
  fit_measuresL <- list()
  attr(n_axes, which = "Names_initial_model") <- 
    colnames(model.matrix(make_formula(1, " " ,initial_trait_model),
                          data = mod$data$dataTraits))
  ### step 2 - end
  stepk <- 0
  while (fit_measures[1, 3] < threshold_P && length(considerk) && 
         stepk < max_step) {
    stepk <- stepk + 1
    fit_measures <- matrix(nrow = length(considerk), ncol = 2)
    colnames(fit_measures) <- c("variance", "pval1")
    rownames(fit_measures) <- considerk
    for (k in seq_along(considerk)) {
      out_FS <- try(dc_CA(
        formulaEnv = mod$formulaEnv,
        formulaTraits = make_formula(n_axes, considered,  considerk[k]),
        response = mod$data$Y, 
        dataEnv = mod$data$dataEnv, 
        dataTraits = mod$data$dataTraits, divideBySiteTotals = divide,
        verbose  = FALSE))
      crit <- if (inherits(out_FS, "try-error") || is.null(out_FS$eigenvalues))
        0 else selection_crit(out_FS,n_axes)
      if (test && crit > 0) {
        if (is.numeric(n_axes)) {
          out_FS <- try(dc_CA(
            formulaEnv = mod$formulaEnv,
            formulaTraits = make_formula("all", considered,  considerk[k]),
            response = mod$data$Y, 
            dataEnv = mod$data$dataEnv, 
            dataTraits = mod$data$dataTraits, divideBySiteTotals = divide,
            verbose  = FALSE))
        }
        if (inherits(out_FS, "try-error") || is.null(out_FS$eigenvalues)) {
          pval <- NA
        } else {
          if (test == 1) {
            an <- if(veganperm) anova(out_FS$RDAonEnv, permutations = permutations) else 
              anova_species(out_FS, permutations = permutations, n_axes = n_axes)$table
            pval <- an$`Pr(>F)`[1]
          } else {
            an <- if(veganperm) anova(out_FS, permutations = permutations) else
              anova(out_FS, permutations = permutations, n_axes = n_axes)
            pval <- an$maxP$`Pr(>F)`[1]
          }
        }
      } else pval <- NA
      fit_measures[k, ] <- c(variance = crit, pval1 = pval)
    }
    pvaladj <- if (test) stats::p.adjust(fit_measures[,"pval1"], 
                                         method = PvalAdjustMethod) else 0
    fit_measures <- cbind(fit_measures, pvaladj)
    id <- order(fit_measures[, "variance"], decreasing = TRUE)
    fit_measures <- fit_measures[id, , drop = FALSE]
    if (verbose) print(round(fit_measures, 5))
    fit_measuresL[[stepk]] <- fit_measures
    best_trait <- rownames(fit_measures)[1]
    considered <- c(considered, consider[best_trait])
    if (verbose) {
      print("considered")
      print(considered)
    }
    considerk <- rownames(fit_measures)[-1]
    fit_measuresL[[stepk]] <- fit_measures
  }
  names(fit_measuresL) <- considered[-1]
  fit_measuresL
  # selected variables with last non-significant value
  tab <- t(sapply(fit_measuresL, function(x) {
    x[1, , drop = FALSE]
  }))
  colnames(tab) <- colnames(fit_measures)
  selected_vars <- if (nrow(tab) == 1) rownames(tab) else rownames(tab)[-nrow(tab)]
  formulaTrait2 <- make_formula(n_axes= "all", initial_trait_model, 
                                paste(selected_vars, collapse = "+"))
  mod2 <- dc_CA(formulaEnv = mod$formulaEnv,
                formulaTraits = formulaTrait2,
                response = mod$data$Y,
                dataEnv = mod$data$dataEnv,
                dataTraits = mod$data$dataTraits,
                divideBySiteTotals = divide,
                verbose = FALSE)
  ret <- list(finalWithOneExtra = tab, model_final = mod2, 
              process = fit_measuresL)
  if (nrow(tab) == 1 && max_step > 1) {
    names(ret)[1] <- if (length(consider) > 1) "final_Non_Significant" else "final"
  }
  return(ret)
}
