# Calculate community weighted means and species niche centroids for double constrained correspondence analysis

Double constrained correspondence analysis (dc-CA) can be calculated
directly from community weighted means (CWMs), with the trait data from
which the CWMs are calculated, and the environmental data and weights
for species and sites (the abundance totals for species and sites).
Statistical testing at the species level requires also the species niche
centroids (SNCs). The function `fCWM_SNC` calculates the CWMs and SNCs
from the trait and environmental data, respectively, using a formula
interface, so as to allow categorical traits and environmental
variables. The resulting object can be set as the `response` argument in
[`dc_CA`](dc_CA.md) so as to give the same output as a call to
[`dc_CA`](dc_CA.md) with the abundance data as `response`, at least up
to sign changes of the axes.

## Usage

``` r
fCWM_SNC(
  response = NULL,
  dataEnv = NULL,
  dataTraits = NULL,
  formulaEnv = NULL,
  formulaTraits = NULL,
  divideBySiteTotals = NULL
)
```

## Arguments

- response:

  matrix, data frame of the abundance data (dimension *n* x *m*) or list
  with community weighted means (CWMs) from `fCWM_SNC`, NULL. If `NULL`,
  the response should be at the left-hand side of `formulaEnv`. See
  Details for analyses starting from community weighted means. Rownames
  of `response`, if any, are carried through.

- dataEnv:

  matrix or data frame of the row predictors, with rows corresponding to
  those in `response`. (dimension *n* x *p*).

- dataTraits:

  matrix or data frame of the column predictors, with rows corresponding
  to the columns in `response`. (dimension *m* x *q*).

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

- divideBySiteTotals:

  logical; default `TRUE` for closing the data by dividing the rows in
  the `response` by their total. However, the default is `FALSE`, when
  the species totals are proportional to `N2*(N-N2)` with `N2` the Hill
  numbers of order 2 of the species and `N` the number of sites, as
  indicator that the response data have been pre-processed to `N2`-based
  marginals using [`ipf2N2`](ipf2N2.md).

## Value

The default returns a list of CWM, SNC, weights, `formulaTraits`,
inertia (weighted variance explained by the traits and by the
environmental variables) and a list of data with elements `dataEnv` and
`dataTraits`.

## Details

The argument `formulaTraits` determines which CWMs are calculated. The
CWMs are calculated from the trait data (non-centered,
non-standardized). With trait covariates, the other predictor traits are
adjusted for the trait covariates by weighted regression, after which
the overall weighted mean trait is added. This has the advantage that
each CWM has the scale of the original trait.

The SNCs are calculated analogously from environmental data.

Empty (all zero) rows and columns in `response` are removed from the
`response` and the corresponding rows from `dataEnv` and `dataTraits`.
Subsequently, any columns with missing values are removed from `dataEnv`
and `dataTraits`. It gives an error (object 'name_of_variable' not
found), if variables with missing entries are specified in `formulaEnv`
and `formulaTraits`.

In the current implementation, `formulaEnv` and `formulaTraits` should
contain variable names as is, *i.e.* transformations of variables in the
formulas gives an error ('undefined columns selected') when the
[`scores`](https://vegandevs.github.io/vegan/reference/scores.html)
function is applied.

## References

Kleyer, M., Dray, S., Bello, F., Lepš, J., Pakeman, R.J., Strauss, B.,
Thuiller, W. & Lavorel, S. (2012) Assessing species and community
functional responses to environmental gradients: which multivariate
methods? Journal of Vegetation Science, 23, 805-821.
[doi:10.1111/j.1654-1103.2012.01402.x](https://doi.org/10.1111/j.1654-1103.2012.01402.x)

ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
double constrained correspondence analysis. Environmental and Ecological
Statistics, 25(2), 171-197.
[doi:10.1007/s10651-017-0395-x](https://doi.org/10.1007/s10651-017-0395-x)

ter Braak C.J.F. and P. Šmilauer (2018). Canoco reference manual and
user's guide: software for ordination (version 5.1x). Microcomputer
Power, Ithaca, USA, 536 pp.

Oksanen, J., et al. (2022) vegan: Community Ecology Package. R package
version 2.6-4. <https://CRAN.R-project.org/package=vegan>.

## See also

[`dc_CA`](dc_CA.md), [`plot.dcca`](plot.dcca.md),
[`scores.dcca`](scores.dcca.md), [`print.dcca`](print.dcca.md) and
[`anova.dcca`](anova.dcca.md)

## Examples

``` r
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

CWMSNC <- fCWM_SNC(formulaEnv = ~ A1 + Moist + Manure + Use + Condition(Mag),
                   formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits)
names(CWMSNC)
#> [1] "CWM"           "SNC"           "formulaEnv"    "formulaTraits"
#> [5] "inertia"       "weights"       "call"          "data"         
#CWMSNC$SNC <- NULL # would give correct dc-CA but no species-level t-values or test
mod <- dc_CA(response = CWMSNC) # formulas and data are in CWMSNC!
#> Step 1: the CCA ordination of the transposed matrix with trait constraints,
#>         useful in itself and also yielding CWMs of the orthonormalized traits for step 2.
#> NULL
#> Step 2: the RDA ordination of CWMs of the orthonormalized traits 
#>         of step 1 with environmental constraints:
#> 
#> Call: rda(formula = out1$CWMs_orthonormal_traits ~ A1 + Moist + Manure +
#> Use + Condition(Mag), data = out1$data$dataEnv)
#> 
#>               Inertia Proportion Rank
#> Total          0.5109     1.0000     
#> Conditional    0.1794     0.3512    3
#> Constrained    0.1377     0.2695    4
#> Unconstrained  0.1938     0.3793    4
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>    RDA1    RDA2    RDA3    RDA4 
#> 0.07351 0.03911 0.02064 0.00442 
#> 
#> Eigenvalues for unconstrained axes:
#>     PC1     PC2     PC3     PC4 
#> 0.11386 0.04942 0.01888 0.01162 
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>         Avg    SDS    VIF  Regr1  tval1
#> MagBF  0.15 0.3571 4.5016 0.1152 0.7224
#> MagHF  0.25 0.4330 2.6715 0.0658 0.5359
#> MagNM  0.30 0.4583 9.5666 0.0046 0.0199
#> A1     4.85 2.1240 1.6967 0.1562 1.5958
#> Moist  2.90 1.7578 1.7658 0.1492 1.4946
#> Manure 1.75 1.4448 9.3847 0.0320 0.1388
#> Use    1.90 0.7681 1.7825 0.1020 1.0172
#>                        Avg     SDS    VIF   Regr1   tval1
#> Seedmass            0.6543  0.6688 1.0784 -0.0724 -0.1662
#> SLA                24.6468  6.3438 1.1888 -0.4217 -0.9226
#> Height             25.1272 15.6848 1.3033  0.5546  1.1589
#> LDMC              244.5084 70.9729 1.1791 -0.8979 -1.9726
#> Lifespanperennial   0.9607  0.1944 1.0964  0.3589  0.8176
#> 
#>                weighted variance
#> traits_explain             0.511
#> env_explain                0.691
#> conditionTE                0.179
#> constraintsTE              0.138
#> attr(,"meaning")
#>                meaning                                                                
#> traits_explain "trait-constrained variation"                                          
#> env_explain    "environment-constrained variation"                                    
#> conditionTE    "trait-constrained variation explained by the condition in formulaEnv" 
#> constraintsTE  "trait-constrained variation explained by the predictors in formulaEnv"
# note that output also gives the environment-constrained inertia,
# which is a bonus compare to the usual way to carry out a dcCA.
anova(mod)
#> $species
#> Species-level permutation test using dc-CA
#> Model: dc_CA(response = CWMSNC) 
#> Residualized predictor permutation
#> 
#>          df ChiSquare      R2      F Pr(>F)
#> dcCA      4   0.13769 0.22615 1.6073  0.185
#> Residual 22   0.47115                      
#> 
#> $sites
#> sites-level permutation test using dc-CA
#> Model: dc_CA(response = CWMSNC) 
#> Residualized predictor permutation
#> 
#>          df ChiSquare      R2      F Pr(>F)  
#> dcCA      4   0.13769 0.41538 2.1316  0.041 *
#> Residual 12   0.19379                        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $maxP
#> Max test combining the community- and species- level tests 
#> Model: dc_CA(response = CWMSNC) 
#> 
#> Taken from the species-level test:
#> Residualized predictor permutation
#> Permutation: free
#> Number of permutations: 999
#> 
#>          df ChiSquare      R2      F Pr(>F)
#> dcCA      4   0.13769 0.22615 1.6073  0.185
#> Residual 22   0.47115                      
#> 
```
