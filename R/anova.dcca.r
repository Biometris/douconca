#' @title Community- and Species-Level Permutation Test in Double Constrained Correspondence Analysis (dc-CA)
#'
#' @description
#' \code{anova.dcca} performs the community- and species-level permutation tests of dc-CA
#' and combines these with the 'max test', which takes the maximum of the \emph{P}-values.
#' The function arguments are similar to (but more restrictive than) those of \code{\link[vegan]{anova.cca}}.
#
#' @param  object  an object from \code{\link{dc_CA}}.
#' @param permutations a list of control values for the permutations for species and sites (species first, sites second,
#' for traits and environment) as
#'  returned by the function \code{\link[permute]{how}}, or
#'  the number of permutations required (default 999,
#'  or a two-vector with the number for the species-level test first and
#'  that for the sites-level second), or
#'  a list of two permutation matrices (again, species first, sites second)
#'  where each row gives the permuted indices.
#' @param  by character \code{"axis"} which performs a series of test, one for each axis,
#'  with the eigenvalue of the axis as  test statistic.
#'  Default: \code{NULL} which set the test statistic to the inertia named \code{constraintsTE}
#'  in the inertia element of \code{\link{dc_CA}}). The interpretation of this inerita is,
#'  at the species-level, the environmentally constrained inertia explained by the traits (without trait covariates) and,
#'  at the community-level, the trait-constrained inertia
#'  explained by the environmental predictors (without covariates).
#'  The default (\code{NULL}) is computationally quicker as it avoids computation of an svd of permuted data sets.
#' @details
#' In the general case of varying site abundance totals (\code{divide.by.site.totals = FALSE})
#' both the community-level test and the species-level test use
#' residualized predictor permutation (ter Braak 2022), so as to ensure that
#' \code{\link{anova.dcca}} is robust against differences
#' in species and site total abundance in the \code{response} (ter Braak & te Beest, 2022).
#' The community-level test uses \code{\link{anova_sites}}.
#' For the species-level test, \code{\link{anova_species}} is used.
#'
#' With equal site weights, obtained with \code{divide.by.site.total = TRUE},
#' the community-level test is obtained by applying \code{anova} to \code{object$RDAonEnv}
#' using \code{\link[vegan]{anova.cca}}.
#' This performs residualized response permutation which performs about equal
#' to residualized predictor permutation in the equi-weight case.
#' The function \code{\link[vegan]{anova.cca}} is implemented in C and therefore
#' a factor of 20 or so quicker than the native R-code used in \code{\link{anova_sites}}.
#' @return
#' A list of 3 of structures as from \code{\link[vegan]{anova.cca}}. The elements are \code{c("species","sites","max")}
#'
#' @references
#' ter Braak, C.J.F. & te Beest, D.E. 2022. Testing environmental effects
#' on taxonomic composition with canonical correspondence analysis:
#' alternative permutation tests are not equal.
#' Environmental and Ecological Statistics. 29 (4), 849-868.
#' https://doi.org/10.1007/s10651-022-00545-4
#'
#' ter Braak, C.J.F. (2022) Predictor versus response permutation
#' for significance testing in weighted regression and redundancy analysis.
#' Journal of statistical computation and simulation, 92, 2041-2059.
#'  https://doi.org/10.1080/00949655.2021.2019256
#' @example demo/dune_test.R
#' @export



anova.dcca <- function(object, permutations = 999, by = NULL){
# object dcca object; permat  a matrix of permutations. if set overrules permuations.
#
  if (!"dcca" %in% class(object)) stop("The first argument must be of class 'dccav', the result of the function dc_CA.")

  if (is.null(by)) by <- "omnibus"
  if (is.na(pmatch(by, c("axis","omnibus")) )) stop(" set argument 'by' to 'axis' or 'NULL'")

  if (length(permutations)==1)  permutations <- as.list(rep(permutations,2)) else
    if (!("list"%in% class(permutations)) ) stop("permutations must be a list of 2 elements specifying species and site permutations")

  for (k in 1:2){
    if (!(is.numeric(permutations[[k]])||"how" %in% class(permutations[[k]])))stop("permutations must be an integer of a list of 2 elements (numbers or of class how)")
  }

    f_species0 <- anova_species(object, by = by, permutations = permutations[[1]]) #
    object1 <- paste("Model:", c(object$call),"\n")
    if (!is.null(f_species0$table)){
    table <- f_species0$table
      #data.frame(f_species0$df, f_species0$ss, c(f_species0$F0, NA), c(f_species0$pval, NA))
    #colnames(table) <- c("Df", varname, "F", "Pr(>F)")

    head <- paste0("Species-level permutation test using dc-CA\n",
                   object1,
                   "Residualized predictor permutation\n",
                   howHead(attr(f_species0,"control") ))

    f_species <-structure(table, heading = head, #Random.seed = seed,
                          control = f_species0$table$how,
                          Random.seed =   attr(f_species0$table,"seed"),
                          control = attr(f_species0$table,"control"),
                          F.perm = attr(f_species0$table,"F.perm"),
                          class = c("anova.cca", "anova", "data.frame"))
    } else {f_species <- NULL}

    if (by=="axis") by1 <- by else by1 <-NULL

    if ("dccav" == class(object)[1]){

      f_sites <- stats::anova(object$RDAonEnv, by = by1, permutations = permutations[[2]])
      rownames(f_sites) <- paste("dcCA", seq_len(nrow(f_sites)), sep ="")
      attr(f_sites,"heading") <- paste0("Community-level equi-weigthed permutation test using vegan::rda\n",
                                        object1,
                                        howHead(attr(f_sites,"control")))
      names(f_sites)[2]<- "ChiSquare"
    } else {
      f_sites0 <- anova_sites(object, by = by1, permutations = permutations[[2]])
      table <- f_sites0$table
     # object1 <- paste("Model:", c(object$call),"\n")
      head <- paste0("sites-level permutation test using dc-CA\n",
                     object1,
                     "Residualized predictor permutation\n",
                     howHead(attr(f_sites0,"control") ))

      f_sites <-structure(table, heading = head, #Random.seed = seed,
                          control = f_sites0$table$how,
                          Random.seed =   attr(f_sites0$table,"seed"),
                          control = attr(f_sites0$table,"control"),
                          F.perm = attr(f_sites0$table,"F.perm"),
                          class = c("anova.cca", "anova", "data.frame"))
    }

    if (!is.null(f_species)){

    if (all(f_sites$`Pr(>F)` > f_species$`Pr(>F)`, na.rm=TRUE)) {
      f_max <- f_sites
      attr(f_max,"heading") <- paste0("Max test combining the community- and species- level tests \n", object1,
                                      "\nTaken from the community-level test:\n",
                                      "Residualized response permutation using vegan::rda\n",
                                      "which performs well in this equi-weight case.\n",
                                      howHead(attr(f_sites,"control")))
    } else if (all(f_sites$`Pr(>F)` <= f_species$`Pr(>F)`, na.rm=TRUE)) {
      f_max <- f_species
      attr(f_max,"heading") <- paste0("Max test combining the community- and species- level tests \n", object1,
                                      "\nTaken from the species-level test:\n",
                                      "Residualized predictor permutation\n",
                                      howHead(attr(f_species,"control")))
    } else { # mix
      id <- f_sites$`Pr(>F)` > f_species$`Pr(>F)`
      id <- id[-length(id)]
      f_max <- f_species
      f_max$R2 <- NULL
      f_max$F[id] <- f_sites$F[id]
      a <- cbind( traitP = f_species$`Pr(>F)`, envP=f_sites$`Pr(>F)` )
      a <- cbind(a, maxP = cummax(apply(a, 1, max)))
      f_max <- cbind(f_max[,c("df","ChiSquare", "F")], a )
      names(f_max)[ncol(f_max)] <- "Pr(>F)"
      head <- paste0("Max test combining the community- and species- level tests \n", object1,
                                      "\na mix the species- (traits) and community- (environment) levels:\n"
                                      )

      fmax<- structure(f_max, heading = head) #Random.seed = seed,

      attr(f_max,"heading") <- head
      class(f_max) <- c("anova.cca", "anova", "data.frame")
    }

      } else f_max <- NULL
    result <- list(species =f_species,
      sites= f_sites,
    max = f_max)

return(result)
}

