# Performs (weighted) double constrained correspondence analysis (dc-CA)

Double constrained correspondence analysis (dc-CA) for analyzing
(multi-)trait (multi-)environment ecological data using library `vegan`
and native R code. It has a `formula` interface which allows to assess,
for example, the importance of trait interactions in shaping ecological
communities. The function `dc_CA` has an option to divide the abundance
data of a site by the site total, giving equal site weights. This
division has the advantage that the multivariate analysis corresponds
with an unweighted (multi-trait) community-level analysis, instead of
being weighted (Kleyer et al. 2012, ter Braak and van Rossum, 2025).

## Usage

``` r
dc_CA(
  formulaEnv = NULL,
  formulaTraits = NULL,
  response = NULL,
  dataEnv = NULL,
  dataTraits = NULL,
  divideBySiteTotals = NULL,
  dc_CA_object = NULL,
  env_explain = TRUE,
  use_vegan_cca = FALSE,
  verbose = TRUE
)
```

## Arguments

- formulaEnv:

  two-sided or one-sided formula for the rows (samples) with row
  predictors in `dataEnv`. The left hand side of the formula is ignored
  if it is specified in the `response` argument. Specify row covariates
  (if any) by adding `+ Condition(covariate-formula)` to `formulaEnv` as
  in [`rda`](https://vegandevs.github.io/vegan/reference/cca.html). The
  `covariate-formula` should not contain a `~` (tilde). Default: `NULL`
  for `~.`, i.e. all variables in `dataEnv` are predictor variables.

- formulaTraits:

  formula or one-sided formula for the columns (species) with column
  predictors in `dataTraits`. When two-sided, the left hand side of the
  formula is not used. Specify column covariates (if any ) by adding
  `+ Condition(covariate-formula)` to `formulaTraits` as in
  [`cca`](https://vegandevs.github.io/vegan/reference/cca.html). The
  `covariate-formula` should not contain a `~` (tilde). Default: `NULL`
  for `~.`, i.e. all variables in `dataTraits` are predictor traits.

- response:

  matrix, data frame of the abundance data (dimension *n* x *m*) or list
  with community weighted means (CWMs) from [`fCWM_SNC`](fCWM_SNC.md),
  NULL. If `NULL`, the response should be at the left-hand side of
  `formulaEnv`. See Details for analyses starting from community
  weighted means. Rownames of `response`, if any, are carried through.

- dataEnv:

  matrix or data frame of the row predictors, with rows corresponding to
  those in `response`. (dimension *n* x *p*).

- dataTraits:

  matrix or data frame of the column predictors, with rows corresponding
  to the columns in `response`. (dimension *m* x *q*).

- divideBySiteTotals:

  logical; default `TRUE` for closing the data by dividing the rows in
  the `response` by their total. However, the default is `FALSE`, when
  the species totals are proportional to `N2*(N-N2)` with `N2` the Hill
  numbers of order 2 of the species and `N` the number of sites, as
  indicator that the response data have been pre-processed to `N2`-based
  marginals using [`ipf2N2`](ipf2N2.md).

- dc_CA_object:

  optional object from an earlier run of this function. Useful if the
  same formula for the columns (`formulaTraits`), `dataTraits` and
  `response` are used with a new formula for the rows. If set, the data
  of the previous run is used and the result of its first step is taken
  for the new analysis and `env_explain` is set to `FALSE`.

- env_explain:

  logical (default `TRUE`) for calculation of the inertia explained by
  the environmental variable (based on a CCA of abundance (with
  `divideBySiteTotals`, if true) on the environmental formula).

- use_vegan_cca:

  default `TRUE`.

- verbose:

  logical for printing a simple summary (default: TRUE)

## Value

A list of `class` `dcca`; that is a list with elements

- CCAonTraits:

  a
  [`cca.object`](https://vegandevs.github.io/vegan/reference/cca.object.html)
  from the [`cca`](https://vegandevs.github.io/vegan/reference/cca.html)
  analysis of the transpose of the closed `response` using formula
  `formulaTraits`.

- formulaTraits:

  the argument `formulaTraits`. If the formula was `~.`, it was changed
  to explicit trait names.

- data:

  a list of `Y`, `dataEnv` and `dataTraits`, after removing empty rows
  and columns in `response` and after closure if
  `divideBySiteTotals = TRUE` and with the corresponding rows in
  `dataEnv` and `dataTraits` removed.

- weights:

  a list of unit-sum weights of columns and rows. The names of the list
  are `c("columns", "rows")`, in that order.

- Nobs:

  number of sites (rows).

- CWMs_orthonormal_traits:

  Community weighted means w.r.t. orthonormalized traits.

- RDAonEnv:

  a [`wrda`](wrda.md) object or
  [`cca.object`](https://vegandevs.github.io/vegan/reference/cca.object.html)
  from the [`wrda`](wrda.md) or, if with equal row weights,
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.html)
  analysis, respectively of the column scores of the `cca`, which are
  the CWMs of orthonormalized traits, using formula `formulaEnv`.

- formulaEnv:

  the argument `formulaEnv`. If the formula was `~.`, it was changed to
  explicit environmental variable names.

- eigenvalues:

  the dc-CA eigenvalues (same as those of the
  [`rda`](https://vegandevs.github.io/vegan/reference/cca.html)
  analysis).

- c_traits_normed0:

  mean, sd, VIF and (regression) coefficients of the traits that define
  the dc-CA axes in terms of the traits with t-ratios missing indicated
  by `NA`s for 'tval1'.

- inertia:

  a one-column matrix with, at most, six inertias (weighted variances):

  - total: the total inertia.

  - conditionT: the inertia explained by the condition in
    `formulaTraits` if present (neglecting row constraints).

  - traits_explain: the trait-structured variation, *i.e.* the inertia
    explained by the traits (without constaints on the rows and
    conditional on the Condition in `formulaTraits`). This is the
    maximum that the row predictors could explain in dc-CA (the sum of
    the last two items is thus less than this value).

  - env_explain: the environmentally structured variation, *i.e.* the
    inertia explained by the environment (without constraints on the
    columns but conditional on the Condition `formulaEnv`). This is the
    maximum that the column predictors could explain in dc-CA (the item
    `constraintsTE` is thus less than this value). The value is `NA`, if
    there is collinearity in the environmental data.

  - conditionTE: the trait-constrained variation explained by the
    condition in `formulaEnv`.

  - constraintsTE: the trait-constrained variation explained by the
    predictors (without the row covariates).

If `verbose` is `TRUE` (or after `out <- print(out)` is invoked) there
are three more items.

- `c_traits_normed`: mean, sd, VIF and (regression) coefficients of the
  traits that define the dc-CA trait axes (composite traits), and their
  optimistic t-ratio.

- `c_env_normed`: mean, sd, VIF and (regression) coefficients of the
  environmental variables that define the dc-CA axes in terms of the
  environmental variables (composite gradients), and their optimistic
  t-ratio.

- `species_axes`: a list with four items

  - `species_scores`: a list with names
    `c("species_scores_unconstrained", "lc_traits_scores")` with the
    matrix with species niche centroids along the dc-CA axes (composite
    gradients) and the matrix with linear combinations of traits.

  - `correlation`: a matrix with inter-set correlations of the traits
    with their SNCs.

  - `b_se`: a matrix with (unstandardized) regression coefficients for
    traits and their optimistic standard errors.

  - `R2_traits`: a vector with coefficient of determination (R2) of the
    SNCs on to the traits. The square-root thereof could be called the
    species-trait correlation in analogy with the species-environment
    correlation in CCA.

- `sites_axes`: a list with four items

  - `site_scores`: a list with names
    `c("site_scores_unconstrained", "lc_env_scores")` with the matrix
    with community weighted means (CWMs) along the dc-CA axes (composite
    gradients) and the matrix with linear combinations of environmental
    variables.

  - `correlation`: a matrix with inter-set correlations of the
    environmental variables with their CWMs.

  - `b_se`: a matrix with (unstandardized) regression coefficients for
    environmental variables and their optimistic standard errors.

  - `R2_env`: a vector with coefficient of determination (R2) of the
    CWMs on to the environmental variables. The square-root thereof has
    been called the species-environmental correlation in CCA.

All scores in the `dcca` object are in scaling `"sites"` (1): the
scaling with *Focus on Case distances* .

## Details

Empty (all zero) rows and columns in `response` are removed from the
`response` and the corresponding rows from `dataEnv` and `dataTraits`.
Subsequently, any columns with missing values are removed from `dataEnv`
and `dataTraits`. It gives an error ('name_of_variable' not found), if
variables with missing entries are specified in `formulaEnv` and
`formulaTraits`.

Computationally, dc-CA can be carried out by a single singular value
decomposition (ter Braak et al. 2018), but it is here computed in two
steps. In the first step, the transpose of the `response` is regressed
on to the traits (the column predictors) using
[`cca`](https://vegandevs.github.io/vegan/reference/cca.html) with
`formulaTraits`. The column scores of this analysis (in scaling 1) are
community weighted means (CWM) of the orthonormalized traits. These are
then regressed on the environmental (row) predictors using
[`wrda`](wrda.md) with `formulaEnv` or using
[`rda`](https://vegandevs.github.io/vegan/reference/cca.html), if site
weights are equal.

A dc-CA can be carried out on, what statisticians call, the sufficient
statistics of the method. This is useful, when the abundance data are
not available or could not be made public in a paper attempting
reproducible research. In this case, `response` should be a list with as
first element community weighted means (e.g. `list(CWM = CWMs)`) with
respect to the traits, and the trait data, and, optionally, further list
elements, for functions related to `dc_CA`. The minimum is a
`list(CWM = CWMs, weight = list(columns = species_weights))` with CWM a
matrix or data.frame, but then `formulaEnv`, `formulaTraits`, `dataEnv`,
`dataTraits` must be specified in the call to `dc_CA`. The function
[`fCWM_SNC`](fCWM_SNC.md) and its example show how to set the `response`
for this and helps to create the `response` from abundance data in these
non-standard applications of dc-CA. Species and site weights, if not set
in `response$weights` can be set by a variable `weight` in the data
frames `dataTraits` and `dataEnv`, respectively, but formulas should
then not be `~.`.

The statistics and scores in the example `dune_dcCA.r`, have been
checked against the results in Canoco 5.15 (ter Braak & Šmilauer, 2018).

## References

Kleyer, M., Dray, S., Bello, F., Lepš, J., Pakeman, R.J., Strauss, B.,
Thuiller, W. & Lavorel, S. (2012) Assessing species and community
functional responses to environmental gradients: which multivariate
methods? Journal of Vegetation Science, 23, 805-821.
[doi:10.1111/j.1654-1103.2012.01402.x](https://doi.org/10.1111/j.1654-1103.2012.01402.x)

ter Braak, CJF, Šmilauer P, and Dray S. (2018). Algorithms and biplots
for double constrained correspondence analysis. Environmental and
Ecological Statistics, 25(2), 171-197.
[doi:10.1007/s10651-017-0395-x](https://doi.org/10.1007/s10651-017-0395-x)

ter Braak C.J.F. and P. Šmilauer (2018). Canoco reference manual and
user's guide: software for ordination (version 5.1x). Microcomputer
Power, Ithaca, USA, 536 pp.

ter Braak, C.J.F. and van Rossum, B. (2025). Linking Multivariate Trait
Variation to the Environment: Advantages of Double Constrained
Correspondence Analysis with the R Package Douconca. Ecological
Informatics, 88.
[doi:10.1016/j.ecoinf.2025.103143](https://doi.org/10.1016/j.ecoinf.2025.103143)

Oksanen, J., et al. (2024). vegan: Community Ecology Package. R package
version 2.6-6.1. <https://CRAN.R-project.org/package=vegan>.

## See also

[`plot.dcca`](plot.dcca.md), [`scores.dcca`](scores.dcca.md),
[`print.dcca`](print.dcca.md) and [`anova.dcca`](anova.dcca.md)

## Examples

``` r
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
abun <- dune_trait_env$comm[, -1]  # must delete "Sites"
mod <- dc_CA(formulaEnv = abun ~ A1 + Moist + Mag + Use + Manure,
             formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits,
       verbose = FALSE)

print(mod) # same output as with verbose = TRUE (the default of verbose).                                     
#> Step 1: the CCA ordination of the transposed matrix with trait constraints,
#>         useful in itself and also yielding CWMs of the orthonormalized traits for step 2.
#> 
#> Call: cca0(formula = formulaTraits, response = tY, data =
#> out0$data$dataTraits)
#> 
#>               Inertia Proportion Rank
#> Total          2.3490     1.0000     
#> Constrained    0.6776     0.2885    5
#> Unconstrained  1.6714     0.7115   19
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>    CCA1    CCA2    CCA3    CCA4    CCA5 
#> 0.26839 0.19597 0.12356 0.07003 0.01967 
#> 
#> Eigenvalues for unconstrained axes:
#> CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 
#>  NA  NA  NA  NA  NA  NA  NA  NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>                        Avg     SDS    VIF   Regr1   tval1
#> SLA                24.6468  6.3438 1.1888  0.4401  4.1794
#> Height             25.1272 15.6848 1.3033  0.1283  1.1640
#> LDMC              244.5084 70.9729 1.1791  0.0278  0.2647
#> Seedmass            0.6543  0.6688 1.0784  0.3708  3.6980
#> Lifespanperennial   0.9607  0.1944 1.0964 -0.0578 -0.5718
#> 
#> Step 2: the RDA ordination of CWMs of the orthonormalized traits 
#>         of step 1 with environmental constraints:
#> 
#> Call: rda(formula = out1$CWMs_orthonormal_traits ~ A1 + Moist + Mag + Use +
#> Manure, data = out1$data$dataEnv)
#> 
#>               Inertia Proportion Rank
#> Total          0.6776     1.0000     
#> Constrained    0.4454     0.6573    5
#> Unconstrained  0.2322     0.3427    5
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>    RDA1    RDA2    RDA3    RDA4    RDA5 
#> 0.23680 0.10903 0.05934 0.03792 0.00233 
#> 
#> Eigenvalues for unconstrained axes:
#>     PC1     PC2     PC3     PC4     PC5 
#> 0.12017 0.05836 0.03455 0.01337 0.00575 
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>         Avg    SDS    VIF   Regr1   tval1
#> A1     4.85 2.1240 1.6967 -0.0605 -0.9172
#> Moist  2.90 1.7578 1.7658  0.3250  4.8293
#> MagBF  0.15 0.3571 4.5016 -0.0475 -0.4421
#> MagHF  0.25 0.4330 2.6715 -0.0156 -0.1890
#> MagNM  0.30 0.4583 9.5666  0.1622  1.0352
#> Use    1.90 0.7681 1.7825  0.0219  0.3239
#> Manure 1.75 1.4448 9.3847 -0.1444 -0.9306
#>                        Avg     SDS    VIF   Regr1   tval1
#> SLA                24.6468  6.3438 1.1888 -0.8196 -3.6933
#> Height             25.1272 15.6848 1.3033 -0.1598 -0.6877
#> LDMC              244.5084 70.9729 1.1791 -0.0562 -0.2542
#> Seedmass            0.6543  0.6688 1.0784 -0.7586 -3.5896
#> Lifespanperennial   0.9607  0.1944 1.0964  0.1006  0.4722
#> 
#>                weighted variance
#> total                      2.349
#> traits_explain             0.678
#> env_explain                1.323
#> constraintsTE              0.445
#> attr(,"meaning")
#>                meaning                                                                
#> total          "total inertia (= weighted variation)"                                 
#> traits_explain "trait-constrained variation"                                          
#> env_explain    "environment-constrained variation"                                    
#> constraintsTE  "trait-constrained variation explained by the predictors in formulaEnv"
anova(mod, by = "axis")
#> $species
#> Species-level permutation test using dc-CA
#> Model: dc_CA(formulaEnv = abun ~ A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> Residualized predictor permutation
#> 
#>          df ChiSquare       R2      F Pr(>F)
#> dcCA1     1   0.23680 0.179002 5.9370  0.118
#> dcCA2     1   0.10903 0.082419 2.7336  0.488
#> dcCA3     1   0.05934 0.044853 1.4877  0.778
#> dcCA4     1   0.03792 0.028663 0.9507  0.799
#> dcCA5     1   0.00233 0.001760 0.0584  0.999
#> Residual 22   0.87749                       
#> 
#> $sites
#> sites-level permutation test using dc-CA
#> Model: dc_CA(formulaEnv = abun ~ A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> Residualized predictor permutation
#> 
#>          df ChiSquare      R2       F Pr(>F)    
#> dcCA1     1  0.236804 0.34946 12.2380  0.001 ***
#> dcCA2     1  0.109032 0.16090  5.6348  0.136    
#> dcCA3     1  0.059336 0.08757  3.0665  0.427    
#> dcCA4     1  0.037919 0.05596  1.9597  0.609    
#> dcCA5     1  0.002329 0.00344  0.1203  1.000    
#> Residual 12  0.232199                           
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $maxP
#> Max test combining the community- and species- level tests 
#> Model: dc_CA(formulaEnv = abun ~ A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> mixing the species- (traits) and community- (environment) levels:
#> 
#>       df ChiSquare traitP  envP Pr(>F)
#> dcCA1  1  0.236804  0.118 0.001  0.118
#> dcCA2  1  0.109032  0.488 0.136  0.488
#> dcCA3  1  0.059336  0.778 0.427  0.778
#> dcCA4  1  0.037919  0.799 0.609  0.799
#> dcCA5  1  0.002329  0.999 1.000  1.000
#> 
# For more demo on testing, see demo dune_test.r

mod_scores <- scores(mod)
# correlation of axes with a variable that is not in the model
scores(mod, display = "cor", scaling = "sym", which_cor = list(NULL, "X_lot"))
#>           dcCA1      dcCA2
#> X_lot 0.1425494 0.06695038
#> attr(,"meaning")
#> [1] "inter set correlation, correlation between environmental variables and the sites scores (CWMs)"

cat("head of unconstrained site scores, with meaning\n")
#> head of unconstrained site scores, with meaning
print(head(mod_scores$sites))
#>             dcCA1       dcCA2
#> Site_1 -1.2640026  0.19395232
#> Site_2 -0.7328822 -0.24929614
#> Site_3 -0.8250139 -0.44760064
#> Site_4 -0.4453407 -0.33741144
#> Site_5 -0.5783587  0.01524487
#> Site_6 -0.5722181  0.24340264

mod_scores_tidy <- scores(mod, tidy = TRUE)
print("names of the tidy scores")
#> [1] "names of the tidy scores"
print(names(mod_scores_tidy))
#> [1] "dcCA1"  "dcCA2"  "score"  "label"  "weight"
cat("\nThe levels of the tidy scores\n")
#> 
#> The levels of the tidy scores
print(levels(mod_scores_tidy$score))
#>  [1] "biplot"                       "biplot_traits"               
#>  [3] "centroids"                    "centroids_traits"            
#>  [5] "constraints_sites"            "constraints_species"         
#>  [7] "correlation"                  "correlation_traits"          
#>  [9] "factorbiplot"                 "factorbiplot_traits"         
#> [11] "intra_set_correlation"        "intra_set_correlation_traits"
#> [13] "regression"                   "regression_traits"           
#> [15] "sites"                        "species"                     
#> [17] "t_values"                     "t_values_traits"             

cat("\nFor illustration: a dc-CA model with a trait covariate\n")
#> 
#> For illustration: a dc-CA model with a trait covariate
mod2 <- dc_CA(formulaEnv = abun ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ SLA + Height + LDMC + Lifespan + Condition(Seedmass),
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits)
#> Step 1: the CCA ordination of the transposed matrix with trait constraints,
#>         useful in itself and also yielding CWMs of the orthonormalized traits for step 2.
#> 
#> Call: cca0(formula = formulaTraits, response = tY, data =
#> out0$data$dataTraits)
#> 
#>               Inertia Proportion Rank
#> Total         2.34899    1.00000     
#> Conditional   0.16671    0.07097    2
#> Constrained   0.51091    0.21750    4
#> Unconstrained 1.67137    0.71153   19
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>    CCA1    CCA2    CCA3    CCA4 
#> 0.21499 0.18622 0.07354 0.03616 
#> 
#> Eigenvalues for unconstrained axes:
#> CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 
#>  NA  NA  NA  NA  NA  NA  NA  NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>                        Avg     SDS    VIF   Regr1   tval1
#> Seedmass            0.6543  0.6688 1.0784  0.0223  0.2225
#> SLA                24.6468  6.3438 1.1888  0.4456  4.2398
#> Height             25.1272 15.6848 1.3033  0.3532  3.2091
#> LDMC              244.5084 70.9729 1.1791 -0.1043 -0.9967
#> Lifespanperennial   0.9607  0.1944 1.0964 -0.0439 -0.4349
#> 
#> Step 2: the RDA ordination of CWMs of the orthonormalized traits 
#>         of step 1 with environmental constraints:
#> 
#> Call: rda(formula = out1$CWMs_orthonormal_traits ~ A1 + Moist + Mag + Use +
#> Manure, data = out1$data$dataEnv)
#> 
#>               Inertia Proportion Rank
#> Total          0.5109     1.0000     
#> Constrained    0.3171     0.6207    4
#> Unconstrained  0.1938     0.3793    4
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>    RDA1    RDA2    RDA3    RDA4 
#> 0.16632 0.08100 0.04766 0.02215 
#> 
#> Eigenvalues for unconstrained axes:
#>     PC1     PC2     PC3     PC4 
#> 0.11386 0.04942 0.01888 0.01162 
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>         Avg    SDS    VIF   Regr1   tval1
#> A1     4.85 2.1240 1.6967 -0.0920 -1.1424
#> Moist  2.90 1.7578 1.7658  0.1002  1.2204
#> MagBF  0.15 0.3571 4.5016  0.0125  0.0951
#> MagHF  0.25 0.4330 2.6715  0.0405  0.4006
#> MagNM  0.30 0.4583 9.5666  0.3052  1.5969
#> Use    1.90 0.7681 1.7825  0.0927  1.1234
#> Manure 1.75 1.4448 9.3847 -0.1567 -0.8277
#>                        Avg     SDS    VIF   Regr1   tval1
#> Seedmass            0.6543  0.6688 1.0784 -0.0929 -0.3874
#> SLA                24.6468  6.3438 1.1888 -1.0644 -4.2276
#> Height             25.1272 15.6848 1.3033 -0.5285 -2.0049
#> LDMC              244.5084 70.9729 1.1791  0.0433  0.1728
#> Lifespanperennial   0.9607  0.1944 1.0964  0.1744  0.7211
#> 
#>                weighted variance
#> total                      2.349
#> conditionT                 0.167
#> traits_explain             0.511
#> env_explain                1.323
#> constraintsTE              0.317
#> attr(,"meaning")
#>                meaning                                                                
#> total          "total inertia (= weighted variation)"                                 
#> conditionT     "variation fitted by the trait condition"                              
#> traits_explain "trait-constrained variation"                                          
#> env_explain    "environment-constrained variation"                                    
#> constraintsTE  "trait-constrained variation explained by the predictors in formulaEnv"

cat("\nFor illustration: a dc-CA model with both environmental and trait covariates\n")
#> 
#> For illustration: a dc-CA model with both environmental and trait covariates
mod3 <- dc_CA(formulaEnv = abun ~ A1 + Moist + Use + Manure + Condition(Mag),
              formulaTraits = ~ SLA + Height + LDMC + Lifespan + Condition(Seedmass),
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits, 
        verbose = FALSE)

cat("\nFor illustration: same model but using dc_CA_object = mod2 for speed, ", 
    "as the trait model and data did not change\n")
#> 
#> For illustration: same model but using dc_CA_object = mod2 for speed,  as the trait model and data did not change
mod3B <- dc_CA(formulaEnv = abun ~ A1 + Moist + Use + Manure + Condition(Mag),
               dataEnv = dune_trait_env$envir,
               dc_CA_object = mod2, 
         verbose= FALSE)
cat("\ncheck on equality of mod3 (from data) and mod3B (from a dc_CA_object)\n",
    "the expected difference is in the component 'call'\n ")
#> 
#> check on equality of mod3 (from data) and mod3B (from a dc_CA_object)
#>  the expected difference is in the component 'call'
#>  

print(all.equal(mod3[-c(5,12)], mod3B[-c(5,12)])) #  only the component call differs
#> [1] TRUE
print(mod3$inertia[-c(3,5),]/mod3B$inertia) #        and mod3 has two more inertia items
#>                weighted variance
#> total                          1
#> conditionT                     1
#> traits_explain                 1
#> conditionTE                    1
#> constraintsTE                  1
#> attr(,"meaning")
#>                meaning                                                                
#> total          "total inertia (= weighted variation)"                                 
#> conditionT     "variation fitted by the trait condition"                              
#> traits_explain "trait-constrained variation"                                          
#> conditionTE    "trait-constrained variation explained by the condition in formulaEnv" 
#> constraintsTE  "trait-constrained variation explained by the predictors in formulaEnv"
```
