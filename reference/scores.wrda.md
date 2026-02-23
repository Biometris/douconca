# Extract results of a weighted redundancy analysis (wrda) or a cca0 object.

This function works very much like the `vegan`
[`scores`](https://vegandevs.github.io/vegan/reference/scores.html)
function, in particular
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html),
but with regression coefficients for predictors.

## Usage

``` r
# S3 method for class 'wrda'
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

  object of class `"wrda"`, *i.e.* result of [`wrda`](wrda.md) or
  [`cca0`](cca0.md).

- ...:

  Other arguments passed to the function (currently ignored).

- choices:

  integer vector of which axes to obtain. Default: all wrda axes.

- display:

  a character vector, one or more of
  `c("all", "species", "sites", "sp", "wa", "lc", "bp", "cor", "ic", "reg", "tval", "cn")`.
  The most items are as in
  [`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html),
  except `"cor"` and `"ic"`, for inter-set and intra-set correlations,
  respectively, and `"tval"` for the (over-optimistic) t-values of the
  regression coefficients.

- scaling:

  numeric (1,2 or 3) or character `"sites", "species" or "symmetric"`.
  Default: "symmetric". Either site- (1) or species- (2) related scores
  are scaled by eigenvalues, and the other set of scores have unit
  weighted mean square or with 3 both are scaled symmetrically to
  weighted mean squares equal to the square root of eigenvalues.
  Negative values are treated as the corresponding positive ones by
  `abs(scaling)`.

- which_cor:

  character vector environmental variables names in the data frames for
  which inter-set correlations must calculated. Default: a character
  ("in_model") for all predictors in the model, including collinear
  variables and levels.

- normed:

  logical (default `TRUE`) giving standardized regression coefficients
  and biplot scores. When `FALSE`, (regular) regression coefficients and
  (unstandardized) biplot scores.

- tidy:

  Return scores that are compatible with `ggplot2`: all variable
  `score`, the names by variable `label`. See weights (in
  [`dc_CA`](dc_CA.md) are in variable `weight`. See
  [`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html).

## Value

A data frame if `tidy = TRUE`. Otherwise, a matrix if a single item is
asked for and a named list of matrices if more than one item is asked
for. The following names can be included:
`c("sites", "constraints_sites", "centroids", "regression", "t_values", "correlation", "intra_set_correlation", "biplot", "species")`.
Each matrix has an attribute `"meaning"` explaining its meaning. With
`tidy = TRUE`, the resulting data frame has attributes `"scaling"` and
`"meaning"`; the latter has two columns: (1) name of score type and (2)
its meaning, usage and interpretation.

An example of the meaning of scores in scaling `"symmetric"` with
`display = "all"`:

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

The statements on optimality for distance interpretations are based on
the `scaling` and the relative magnitude of the dc-CA eigenvalues of the
chosen axes.

## Details

The function is modeled after
[`scores.cca`](https://vegandevs.github.io/vegan/reference/plot.cca.html).

An example of which_cor is: `which_cor = c("acidity", "humidity")`

## Examples

``` r
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
response <- dune_trait_env$comm[, -1]  # must delete "Sites"

w <- rep(1, 20) 
w[1:10] <- 8 
w[17:20] <- 0.5

object <- wrda(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure),
               data = dune_trait_env$envir, 
               weights = w)
object # Proportions equal to those Canoco 5.15
#> 
#> Call: wrda(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure),
#> data = dune_trait_env$envir, weights = w)
#> 
#>               Inertia Proportion Rank
#> Total         65.7007     1.0000     
#> Conditional    7.3839     0.1124    2
#> Constrained   36.4952     0.5555    6
#> Unconstrained 21.8217     0.3321   19
#> 
#> Inertia is weighted variance
#> 
#> Eigenvalues for constrained axes:
#>  wRDA1  wRDA2  wRDA3  wRDA4  wRDA5  wRDA6 
#> 18.133  8.233  4.047  3.155  2.238  0.689 
#> 
#> Eigenvalues for unconstrained axes:
#> wPCA1 wPCA2 wPCA3 wPCA4 wPCA5 wPCA6 wPCA7 wPCA8 
#>    NA    NA    NA    NA    NA    NA    NA    NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>           Avg    SDS     VIF   Regr1   tval1
#> Manure 2.4659 1.2056 10.5656 -0.2010 -0.1344
#> A1     4.1466 1.3786  1.4841  0.1013  0.1807
#> Moist  2.1761 1.4762  1.4883 -4.1077 -7.3195
#> MagBF  0.1932 0.3948  6.8216  0.7873  0.6553
#> MagHF  0.4545 0.4979  5.8054  1.8392  1.6594
#> MagNM  0.0455 0.2083  4.8878  0.1327  0.1304
#> Use    1.9205 0.7107  2.6787  0.2101  0.2791
#> 

mod_scores <- scores(object, display = "all")
scores(object, which_cor = c("A1", "X_lot"), display = "cor")
#>            wRDA1      wRDA2
#> A1    -0.2365172 -0.3510007
#> X_lot -0.4987793  0.4250562
#> attr(,"meaning")
#> [1] "inter set correlation, correlation between environmental variables and the sites scores (CWMs)"
anova(object)
#> $table
#> Permutation test for weighted reduncancy analysis
#> Model: wrda(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure), data = dune_trait_env$envir, weights = w) 
#> Residualized predictor permutation
#> Permutation: free
#> Number of permutations: 999
#> 
#>          df Variance      R2      F Pr(>F)  
#> wRDA      6   36.495 0.62581 3.3449  0.017 *
#> Residual 12   21.822                        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $eigenvalues
#>      wRDA1      wRDA2      wRDA3      wRDA4      wRDA5      wRDA6 
#> 18.1328543  8.2332883  4.0467247  3.1546579  2.2383295  0.6893052 
#> 
```
