# Extract results of a double constrained correspondence analysis (dc-CA)

This function works very much like the `vegan`
[`scores`](https://vegandevs.github.io/vegan/reference/scores.html)
function, in particular
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html),
with the additional results such as regression coefficients and linear
combinations of traits `('reg_traits', 'lc_traits')`. All scores from CA
obey the so called transition formulas and so do the scores of CCA and
dc-CA. The differences are, for CCA, that the linear combinations of
environmental variables (the *constrained* site scores) replace the
usual (*unconstrained*) site scores, and for dc-CA, that the linear
combinations of traits (the *constrained* species scores) also replace
the usual (*unconstrained*) species scores in the transition formulas.

## Usage

``` r
# S3 method for class 'dcca'
scores(
  x,
  ...,
  choices = 1:2,
  display = "all",
  scaling = "sym",
  which_cor = "in model",
  normed = TRUE,
  tidy = FALSE
)
```

## Arguments

- x:

  object of class `"dcca"`, *i.e.* result of [`dc_CA`](dc_CA.md).

- ...:

  Other arguments passed to the function (currently ignored).

- choices:

  integer vector of which axes to obtain. Default: all dc-CA axes.

- display:

  a character vector, one or more of
  `c("all", "species", "sites", "sp", "wa", "lc", "bp", "cor", "ic", "reg", "tval", "cn", "lc_traits", "reg_traits", "tval_traits", "cor_traits", "ic_traits", "bp_traits", "cn_traits")`.
  The most items are as in
  [`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html),
  except `"cor"` and `"ic"`, for inter-set and intra-set correlations,
  respectively, and `"tval"` for the (over-optimistic) t-values of the
  regression coefficients. The remaining scores are analogous scores for
  species and traits.

- scaling:

  numeric (1,2 or 3) or character `"sites", "species" or "symmetric"`.
  Default: "symmetric". Either site- (1) or species- (2) related scores
  are scaled by eigenvalues, and the other set of scores have unit
  weighted mean square or with 3 both are scaled symmetrically to
  weighted mean squares equal to the square root of eigenvalues.
  Negative values are treated as the corresponding positive ones by
  `abs(scaling)`.

- which_cor:

  character or list of trait and environmental variables names (in this
  order) in the data frames for which inter-set correlations must
  calculated. Default: a character ("in_model") for all traits and
  variables in the model, including collinear variables and levels.

- normed:

  logical (default `TRUE`) giving standardized regression coefficients
  and biplot scores. When `FALSE`, (regular) regression coefficients and
  (unstandardized) biplot scores.

- tidy:

  Return scores that are compatible with `ggplot2`: all scores are in a
  single data.frame, score type is identified by factor variable
  `score`, the names by variable `label`, and species weights (in
  [`dc_CA`](dc_CA.md) are in variable `weight`. See
  [`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html).

## Value

A data frame if `tidy = TRUE`. Otherwise, a matrix if a single item is
asked for and a named list of matrices if more than one item is asked
for. The following names can be included:
`c("sites", "constraints_sites", "centroids", "regression", "t_values", "correlation", "intra_set_correlation", "biplot", "species", "constraints_species", "regression_traits", "t_values_traits", "correlation_traits", "intra_set_correlation_traits", "biplot_traits", "centroids_traits")`.
Each matrix has an attribute `"meaning"` explaining its meaning. With
`tidy = TRUE`, the resulting data frame has attributes `"scaling"` and
`"meaning"`; the latter has two columns: (1) name of score type and (2)
its meaning, usage and interpretation.

An example of the meaning of scores in scaling `"symmetric"` with
`display ="all"`:

- sites:

  CMWs of the trait axes (constraints species) in scaling 'symmetric'
  optimal for biplots and, almost so, for inter-site distances.

- constraints_sites:

  linear combination of the environmental predictors and the covariates
  (making the ordination axes orthogonal to the covariates) in scaling
  'symmetric' optimal for biplots and, almost so, for inter-site
  distances.

- regression:

  mean, sd, VIF, standardized regression coefficients and their
  optimistic t-ratio in scaling 'symmetric'.

- t_values:

  t-values of the coefficients of the regression of the CWMs of the
  trait composite on to the environmental variables

- correlation:

  inter set correlation, correlation between environmental variables and
  the sites scores (CWMs)

- intra_set_correlation:

  intra set correlation, correlation between environmental variables and
  the dc-ca axis (constrained sites scores)

- biplot:

  biplot scores of environmental variables for display with
  biplot-traits for fourth-corner correlations in scaling 'symmetric'.

- centroids:

  environmental category means of the site scores in scaling 'symmetric'
  optimal for biplots and, almost so, for inter-environmental category
  distances.

- species:

  SNC on the environmental axes (constraints sites) in scaling
  'symmetric' optimal for biplots and, almost so, for inter-species
  distances.

- constraints_species:

  linear combination of the traits and the trait covariates (making the
  ordination axes orthogonal to the covariates) in scaling 'symmetric'
  optimal for biplots and, almost so, for inter-species distances.

- regression_traits:

  mean, sd, VIF, standardized regression coefficients and their
  optimistic t-ratio in scaling 'symmetric'.

- t_values_traits:

  t-values of the coefficients of the regression of the SNCs along a
  dc-CA axis on to the traits

- correlation_traits:

  inter set correlation, correlation between traits and the species
  scores (SNCs)

- intra_set_correlation_traits:

  intra set correlation, correlation between traits and the dc-ca axis
  (constrained species scores)

- biplot_traits:

  biplot scores of traits for display with biplot scores for
  fourth-corner correlation in scaling 'symmetric'.

- centroids_traits:

  trait category means of the species scores in scaling 'symmetric'
  optimal for biplots and, almost so, for inter-trait category
  distances.

The statements on optimality for distance interpretations are based on
the `scaling` and the relative magnitude of the dc-CA eigenvalues of the
chosen axes.

## Details

The function is modeled after
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html).

The t-ratios are taken from a multiple regression of the unconstrained
species (or site) scores on to the traits (or environmental variables).

An example of `which_cor` is:
`which_cor = list(traits = "SLA", env = c("acidity", "humidity"))`.

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
#> dcCA1     1   0.23680 0.179002 5.9370  0.105
#> dcCA2     1   0.10903 0.082419 2.7336  0.519
#> dcCA3     1   0.05934 0.044853 1.4877  0.768
#> dcCA4     1   0.03792 0.028663 0.9507  0.790
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
#> dcCA2     1  0.109032 0.16090  5.6348  0.152    
#> dcCA3     1  0.059336 0.08757  3.0665  0.404    
#> dcCA4     1  0.037919 0.05596  1.9597  0.632    
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
#> dcCA1  1  0.236804  0.105 0.001  0.105
#> dcCA2  1  0.109032  0.519 0.152  0.519
#> dcCA3  1  0.059336  0.768 0.404  0.768
#> dcCA4  1  0.037919  0.790 0.632  0.790
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
