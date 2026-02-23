# Utility function: community-level permutation test in Double Constrained Correspondence Analysis (dc-CA)

`anova_sites` performs the community-level permutation test of dc-CA
when site weights vary. The test uses residual predictor permutation
(ter Braak 2022), which is robust against differences in site total
abundance in the `response` in [`dc_CA`](dc_CA.md) (ter Braak & te
Beest, 2022). The arguments of the function are similar to those of
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html),
but more restricted. With equal site-totals as in [`dc_CA`](dc_CA.md),
`anova(object$RDAonEnv)` is much faster.

## Usage

``` r
anova_sites(
  object,
  permutations = 999,
  rpp = TRUE,
  n_axes = "all",
  max_axis = 10,
  by = NULL
)
```

## Arguments

- object:

  an object from [`dc_CA`](dc_CA.md).

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required (default 999), or a permutation matrix
  where each row gives the permuted indices.

- rpp:

  Logical indicating residual predictor permutation (default `TRUE`).
  When `FALSE`, residual response permutation is used.

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

- max_axis:

  maximum number of axes to test if `by = "axis"`.

- by:

  character `"axis"` which sets the test statistic to the first
  eigenvalue of the dc-CA model. Default: `NULL` which sets the test
  statistic to the inertia (sum of all double constrained eigenvalues;
  named `constraintsTE` in the inertia element of [`dc_CA`](dc_CA.md)).
  This is the trait constrained inertia explained by the environmental
  predictors (without covariates), which is equal to the
  environmentally-constrained inertia explained by the traits (without
  trait covariates). The default is quicker computationally as it avoids
  computation of an svd of permuted data sets.

## Value

A list with two elements with names `table` and `eigenvalues`. The
`table` is as from
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html)
and `eigenvalues` gives the dc-CA eigenvalues.

## Details

The algorithm is analogous to that of [`anova.wrda`](anova.wrda.md). The
function is used in [`anova.dcca`](anova.dcca.md).

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
#> dcCA      5   0.44542 0.3367 2.2335  0.012 *
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
#> dcCA      5   0.44542 0.3367 2.2335  0.012 *
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
#> dcCA      5   0.44542 0.3367 2.2335  0.011 *
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
#> 0.2368039 0.0110000 
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
