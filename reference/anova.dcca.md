# Community- and Species-Level Permutation Test in Double Constrained Correspondence Analysis (dc-CA)

`anova.dcca` performs the community- and species-level permutation tests
of dc-CA and combines these in the 'max test', which takes the maximum
of the *P*-values. The function arguments are similar to (but more
restrictive than) those of
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html).

## Usage

``` r
# S3 method for class 'dcca'
anova(
  object,
  ...,
  rpp = TRUE,
  permutations = 999,
  max_axis = 10,
  by = c("omnibus", "axis"),
  n_axes = "all"
)
```

## Arguments

- object:

  an object from [`dc_CA`](dc_CA.md).

- ...:

  unused.

- rpp:

  Logical indicating residual predictor permutation (default `TRUE`).
  When `FALSE`, residual response permutation is used.

- permutations:

  a list of control values for the permutations for species and sites
  (species first, sites second, for traits and environment) as returned
  by the function [`how`](https://rdrr.io/pkg/permute/man/how.html), or
  the number of permutations required (default 999, or a two-vector with
  the number for the species-level test first and that for the
  sites-level second), or a list of two permutation matrices (again,
  species first, sites second) where each row gives the permuted
  indices.

- max_axis:

  maximum number of axes to test if `by = "axis"`.

- by:

  character `"axis"` which performs a series of tests, one for each
  axis, with the eigenvalue of the axis as test statistic. Default:
  `NULL` which sets the test statistic to the inertia (sum of all double
  constrained eigenvalues; named `constraintsTE` in the inertia element
  of [`dc_CA`](dc_CA.md).

  The interpretation of this inertia is, at the species-level, the
  environmentally constrained inertia explained by the traits (without
  trait covariates) and, at the community-level, the trait-constrained
  inertia explained by the environmental predictors (without
  covariates). The default (`NULL`) is computationally quicker as it
  avoids computation of an svd of permuted data sets.

- n_axes:

  number of axes used in the test statistic (default: `"all"`). Example,
  the test statistic is the sum of the first two eigenvalues, if
  `n_axes=2`. With a numeric `n_axes` and model `~X + Condition(Z)`, the
  residuals of `X` with respect to `Z` are permuted with a test
  statistic equal to the sum of the first `n_axes` eigenvalues of the
  fitted `Y` in the model `Y ~ X + Z`, with `Y` the response in the
  model. In the default `"all"`, the test statistic is all eigenvalues
  of the model `Y ~ X|Z`, *i.e.* the effects of `X` after adjustment for
  the effects on `Y` of `Z`. If `by = "axis"`, the setting of `n_axes`
  is ignored.

## Value

A list of 3 of structures as from
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html).
The elements are `c("species", "sites", "maxP")`

## Details

In the general case of varying site abundance totals
(`divideBySiteTotals = FALSE`) both the community-level test and the
species-level test use residualized predictor permutation (ter Braak
2022), so as to ensure that `anova.dcca` is robust against differences
in species and site total abundance in the `response` (ter Braak & te
Beest, 2022). The community-level test uses
[`anova_sites`](anova_sites.md). For the species-level test,
[`anova_species`](anova_species.md) is used.

With equal site weights, obtained with `divide.by.site.total = TRUE`,
the community-level test is obtained by applying `anova` to
`object$RDAonEnv` using
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html).
This performs residualized response permutation which performs about
equal to residualized predictor permutation in the equi-weight case. The
function
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html)
is implemented in C and therefore a factor of 20 or so quicker than the
native R-code used in [`anova_sites`](anova_sites.md).

The `n_axes` argument is new, and can be considered experimental.

## References

ter Braak, C.J.F. & te Beest, D.E. 2022. Testing environmental effects
on taxonomic composition with canonical correspondence analysis:
alternative permutation tests are not equal. Environmental and
Ecological Statistics. 29 (4), 849-868.
[doi:10.1007/s10651-022-00545-4](https://doi.org/10.1007/s10651-022-00545-4)

ter Braak, C.J.F. (2022) Predictor versus response permutation for
significance testing in weighted regression and redundancy analysis.
Journal of statistical computation and simulation, 92, 2041-2059.
[doi:10.1080/00949655.2021.2019256](https://doi.org/10.1080/00949655.2021.2019256)

## Examples

``` r
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

mod <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
             formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits,
             verbose = FALSE)
anova(mod)
#> $species
#> Species-level permutation test using dc-CA
#> Model: dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, response = dune_trait_env$comm[, -1], dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> Residualized predictor permutation
#> 
#>          df ChiSquare     R2      F Pr(>F)  
#> dcCA      5   0.44542 0.3367 2.2335  0.019 *
#> Residual 22   0.87749                       
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $sites
#> sites-level permutation test using dc-CA
#> Model: dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, response = dune_trait_env$comm[, -1], dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> Residualized predictor permutation
#> 
#>          df ChiSquare      R2      F Pr(>F)    
#> dcCA      7   0.44542 0.65733 3.2885  0.001 ***
#> Residual 12   0.23220                          
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $maxP
#> Max test combining the community- and species- level tests 
#> Model: dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, response = dune_trait_env$comm[, -1], dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> 
#> Taken from the species-level test:
#> Residualized predictor permutation
#> Permutation: free
#> Number of permutations: 999
#> 
#>          df ChiSquare     R2      F Pr(>F)  
#> dcCA      5   0.44542 0.3367 2.2335  0.019 *
#> Residual 22   0.87749                       
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

a_species <- anova_species(mod)
a_species
#> $table
#> Species-level permutation test using dc-CA
#> Model: dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, response = dune_trait_env$comm[, -1], dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> Residualized predictor permutation
#> Permutation: free
#> Number of permutations: 999
#> 
#>          df ChiSquare     R2      F Pr(>F)  
#> dcCA      5   0.44542 0.3367 2.2335  0.012 *
#> Residual 22   0.87749                       
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $eigenvalues
#>      dcCA1      dcCA2      dcCA3      dcCA4      dcCA5 
#> 0.23680387 0.10903220 0.05933626 0.03791909 0.00232876 
#> 
# anova_species can be used for manual forward selection of
# trait variables, as done for environmental variables in the demo 
# dune_FS_dcCA.r, based on the first eigenvalue and its significance
# and adding the first axis of the constrained species scores from mod to 
# the Condition of a new mod.
(eig1_and_pval <- c(eig = a_species$eig[1], pval = a_species$table$`Pr(>F)`[1]))
#> eig.dcCA1      pval 
#> 0.2368039 0.0120000 
a_species$eig
#>      dcCA1      dcCA2      dcCA3      dcCA4      dcCA5 
#> 0.23680387 0.10903220 0.05933626 0.03791909 0.00232876 
anova_sites(mod)
#> $table
#> Community-level permutation test using dc-CA
#> Model: dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure, formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan, response = dune_trait_env$comm[, -1], dataEnv = dune_trait_env$envir, dataTraits = dune_trait_env$traits, verbose = FALSE) 
#> Residualized predictor permutation
#> Permutation: free
#> Number of permutations: 999
#> 
#>          df ChiSquare      R2      F Pr(>F)    
#> dcCA      7   0.44542 0.65733 3.2885  0.001 ***
#> Residual 12   0.23220                          
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $eigenvalues
#>      dcCA1      dcCA2      dcCA3      dcCA4      dcCA5 
#> 0.23680387 0.10903220 0.05933626 0.03791909 0.00232876 
#> 
```
