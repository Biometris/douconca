# Forward selection of predictor variables using wrda or cca0

Forward selection of predictor variables using wrda or cca0

## Usage

``` r
# S3 method for class 'wrda'
FS(
  mod,
  ...,
  consider = NULL,
  permutations = 999,
  n_axes = "all",
  initial_model = "1",
  factor2categories = TRUE,
  test = TRUE,
  threshold_P = 0.1,
  PvalAdjustMethod = "holm",
  max_step = 10,
  verbose = FALSE
)
```

## Arguments

- mod:

  initial wrda or cca0 model with at least on predictor variable,

- ...:

  unused.

- consider:

  character vector of names in `mod$data` to consider for addition.

- permutations:

  a list of control values for the permutations as returned by the
  function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
  number of permutations required (default 999), or a permutation matrix
  where each row gives the permuted indices.

- n_axes:

  number of eigenvalues to select upon. The sum of `n_axes` eigenvalues
  is taken as criterion. Default `"full"` for selection without
  dimension reduction to `n_axes`. If `n_axes =1`, selection is on the
  first eigenvalue for selection of variables that form an optimal
  one-dimensional model.

- initial_model:

  character specifying what should be inside `Condition()`. Default:
  `"1"` (nothing, the intercept only). Examples: `"region"` for a
  within-region analysis or `"A*B"` for a within analysis specified by
  the interaction of factors `A` and `B`, with `region, A, B` in the
  data.

- factor2categories:

  logical, default `TRUE`, to convert factors to their categories
  (set(s) of indicator values). If `FALSE`, the selection uses, the fit
  of a factor divided by its number of categories minus 1.

- test:

  logical; default: `TRUE`.

- threshold_P:

  significance level, after adjustment for testing multiplicity, for
  addition of a variable to the model.

- PvalAdjustMethod:

  method for correction for multiple testing in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html), default `"holm"`,
  which is an improved version Bonferroni.

- max_step:

  maximal number of variables selected.

- verbose:

  show progress, default: `TRUE`.

## Value

list with three elements: `final...` with selected variables and
`model_final` and `process` with account of the selection process If
`is.numeric(n_axes)`, then the variance in the returned table is the sum
of the n_axes eigenvalues of the current model (all variables so far
included).

## Details

The selection is on the basis of the additional fit (inertia) of a
variable given the variables already in the model.

The names in `consider` may include transformations of predictor
variables, such as `log(.)`, if `consider` does not include factors or
if `factor2categories=FALSE`. If `consider` does include factors, such
transformations give in a error in the default setting
(`factor2categories=TRUE`).

## See also

[`cca0`](cca0.md), [`wrda`](wrda.md) and [`FS.dcca`](FS.dcca.md)

## Examples

``` r
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
abun <- dune_trait_env$comm[, -1] # must delete "Sites"

mod <- dc_CA(formulaEnv = abun ~ Moist + Mag,
             formulaTraits = ~ F + R + N + L,
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits,
             verbose = FALSE)

# selection of traits with environmental model of mod (~ Moist+Mag)
out1 <- FS(mod, consider = c("F", "R", "N", "L"), 
           select = "traits", verbose = FALSE) 

names(out1)
#> [1] "finalWithOneExtra" "model_final"       "formula"          
#> [4] "process"          
out1$finalWithOneExtra
#>   df   variance  pval pvaladj
#> F  1 0.28633635 0.001   0.004
#> N  1 0.24243384 0.001   0.003
#> L  1 0.03114819 0.137   0.274
out1$model_final
#> Step 1: the CCA ordination of the transposed matrix with trait constraints,
#>         useful in itself and also yielding CWMs of the orthonormalized traits for step 2.
#> 
#> Call: cca0(formula = formulaTraits, response = tY, data =
#> out0$data$dataTraits)
#> 
#>               Inertia Proportion Rank
#> Total          2.3490     1.0000     
#> Constrained    0.7657     0.3260    2
#> Unconstrained  1.5833     0.6740   19
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>   CCA1   CCA2 
#> 0.4576 0.3081 
#> 
#> Eigenvalues for unconstrained axes:
#> CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 
#>  NA  NA  NA  NA  NA  NA  NA  NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>      Avg    SDS    VIF   Regr1   tval1
#> F 6.1126 1.6685 1.0575  0.6413  8.7440
#> N 5.1674 1.6303 1.0575 -0.1125 -1.5334
#> 
#> Step 2: the RDA ordination of CWMs of the orthonormalized traits 
#>         of step 1 with environmental constraints:
#> 
#> Call: rda(formula = out1$CWMs_orthonormal_traits ~ Moist + Mag, data =
#> out1$data$dataEnv)
#> 
#>               Inertia Proportion Rank
#> Total          0.7657     1.0000     
#> Constrained    0.5288     0.6906    2
#> Unconstrained  0.2370     0.3094    2
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>    RDA1    RDA2 
#> 0.31310 0.21567 
#> 
#> Eigenvalues for unconstrained axes:
#>     PC1     PC2 
#> 0.17114 0.06581 
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>        Avg    SDS    VIF  Regr1  tval1
#> Moist 2.90 1.7578 1.3001 0.4769 4.6192
#> MagBF 0.15 0.3571 1.4553 0.0075 0.0684
#> MagHF 0.25 0.4330 1.4214 0.0527 0.4879
#> MagNM 0.30 0.4583 1.4393 0.1948 1.7933
#>      Avg    SDS    VIF  Regr1   tval1
#> F 6.1126 1.6685 1.0575  0.726  4.9779
#> N 5.1674 1.6303 1.0575 -0.539 -3.6956
#> 
#>                weighted variance
#> total                      2.349
#> traits_explain             0.766
#> env_explain                0.967
#> constraintsTE              0.529
#> attr(,"meaning")
#>                meaning                                                                
#> total          "total inertia (= weighted variation)"                                 
#> traits_explain "trait-constrained variation"                                          
#> env_explain    "environment-constrained variation"                                    
#> constraintsTE  "trait-constrained variation explained by the predictors in formulaEnv"

# selection of environmental variables with trait model of mod (~ F + R + N + L)
out2 <- FS(mod, consider =  c("A1", "Moist", "Mag", "Use", "Manure"), 
           select= "env", verbose = FALSE) 

names(out2)
#> [1] "finalWithOneExtra" "model_final"       "formula"          
#> [4] "process"          
out2$finalWithOneExtra
#>        df   variance  pval pvaladj
#> Moist   1 0.33061402 0.001   0.008
#> Manure  1 0.20050625 0.001   0.007
#> A1      1 0.04503044 0.085   0.510
out2$model_final
#> Step 1: the CCA ordination of the transposed matrix with trait constraints,
#>         useful in itself and also yielding CWMs of the orthonormalized traits for step 2.
#> 
#> Call: cca0(formula = formulaTraits, response = tY, data =
#> out0$data$dataTraits)
#> 
#>               Inertia Proportion Rank
#> Total          2.3490     1.0000     
#> Constrained    0.8701     0.3704    4
#> Unconstrained  1.4789     0.6296   19
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>   CCA1   CCA2   CCA3   CCA4 
#> 0.4771 0.3122 0.0580 0.0228 
#> 
#> Eigenvalues for unconstrained axes:
#> CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 
#>  NA  NA  NA  NA  NA  NA  NA  NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>      Avg    SDS    VIF   Regr1   tval1
#> F 6.1126 1.6685 1.1323  0.6048  8.2493
#> R 5.8061 1.1894 1.7564  0.0170  0.1865
#> N 5.1674 1.6303 1.7822 -0.1358 -1.4764
#> L 7.3839 0.8774 1.1378  0.1424  1.9382
#> 
#> Step 2: the RDA ordination of CWMs of the orthonormalized traits 
#>         of step 1 with environmental constraints:
#> 
#> Call: rda(formula = out1$CWMs_orthonormal_traits ~ Moist + Manure +
#> Condition(1), data = out1$data$dataEnv)
#> 
#>               Inertia Proportion Rank
#> Total          0.8701     1.0000     
#> Constrained    0.5311     0.6104    2
#> Unconstrained  0.3389     0.3896    4
#> 
#> Inertia is variance
#> 
#> Eigenvalues for constrained axes:
#>   RDA1   RDA2 
#> 0.3331 0.1980 
#> 
#> Eigenvalues for unconstrained axes:
#>     PC1     PC2     PC3     PC4 
#> 0.18449 0.10693 0.02887 0.01864 
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>         Avg    SDS    VIF   Regr1   tval1
#> Moist  2.90 1.7578 1.0645  0.5519  6.2344
#> Manure 1.75 1.4448 1.0645 -0.0809 -0.9142
#>      Avg    SDS    VIF   Regr1   tval1
#> F 6.1126 1.6685 1.1323  0.7661  5.4523
#> R 5.8061 1.1894 1.7564 -0.1975 -1.1287
#> N 5.1674 1.6303 1.7822 -0.2553 -1.4483
#> L 7.3839 0.8774 1.1378  0.3167  2.2482
#> 
#>                weighted variance
#> total                      2.349
#> traits_explain             0.870
#> env_explain                0.751
#> constraintsTE              0.531
#> attr(,"meaning")
#>                meaning                                                                
#> total          "total inertia (= weighted variation)"                                 
#> traits_explain "trait-constrained variation"                                          
#> env_explain    "environment-constrained variation"                                    
#> constraintsTE  "trait-constrained variation explained by the predictors in formulaEnv"

# selection of environmental variables without a trait model 
# i.e. with a single constraint
mod3 <- cca0(mod$data$Y ~ Moist, data = mod$data$dataEnv)
out3 <- FS(mod3, consider = c("A1", "Moist", "Mag", "Use", "Manure"), 
           threshold_P = 0.05)

out3$finalWithOneExtra
#>       df  variance  pval pvaladj
#> Moist  1 0.4672020 0.001   0.008
#> MagNM  1 0.3088722 0.001   0.007
#> A1     1 0.1841344 0.037   0.222
out3$model_final
#> 
#> Call: FUN(formula = formula2, data = mydata, traceonly = FALSE, cca_object
#> = mod)
#> 
#>               Inertia Proportion Rank
#> Total          2.3490     1.0000     
#> Constrained    0.7761     0.3304    2
#> Unconstrained  1.5729     0.6696   19
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>   CCA1   CCA2 
#> 0.4803 0.2958 
#> 
#> Eigenvalues for unconstrained axes:
#> CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 
#>  NA  NA  NA  NA  NA  NA  NA  NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>       Avg    SDS    VIF   Regr1   tval1
#> Moist 2.9 0.3931 1.1374 -0.5995 -9.0496
#> MagNM 0.3 0.1025 1.1374 -0.1970 -2.9746
#> 

# selection of traits without an environmental model 
#                         i.e. with a single constraint
tY <- t(mod$data$Y)

mod4 <- cca0(tY ~ L, data = mod$data$dataTraits)

names(mod$data$dataTraits)
#>  [1] "Species"      "Species_abbr" "SLA"          "Height"       "LDMC"        
#>  [6] "Seedmass"     "Lifespan"     "F"            "R"            "N"           
#> [11] "L"           
out4 <- FS(mod4, 
           consider =  c("SLA", "Height", "LDMC", "Seedmass", "Lifespan", 
                         "F", "R", "N", "L"))

out4$finalWithOneExtra
#>          df  variance  pval pvaladj
#> F         1 0.4536777 0.001   0.010
#> N         1 0.3120442 0.001   0.009
#> Seedmass  1 0.1060672 0.047   0.376
out4$model_final
#> 
#> Call: FUN(formula = formula2, data = mydata, traceonly = FALSE, cca_object
#> = mod)
#> 
#>               Inertia Proportion Rank
#> Total          2.3490     1.0000     
#> Constrained    0.7657     0.3260    2
#> Unconstrained  1.5833     0.6740   19
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>   CCA1   CCA2 
#> 0.4576 0.3081 
#> 
#> Eigenvalues for unconstrained axes:
#> CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 
#>  NA  NA  NA  NA  NA  NA  NA  NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>      Avg    SDS    VIF   Regr1   tval1
#> F 6.1126 1.6685 1.0575  0.6413  8.7440
#> N 5.1674 1.6303 1.0575 -0.1125 -1.5334
#> 
```
