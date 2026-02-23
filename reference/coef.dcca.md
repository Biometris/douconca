# Coefficients of double-constrained correspondence analysis (dc-CA)

Fourth-corner coefficients and regression coefficients (of full or
reduced rank) to predict traits from environment, environment from
traits and response from trait and environment data.

## Usage

``` r
# S3 method for class 'dcca'
coef(
  object,
  ...,
  type = c("fourth_corner", "all_reg", "env2traits_reg", "traits2env_reg"),
  rank = "full",
  normed = TRUE
)
```

## Arguments

- object:

  return value of [`dc_CA`](dc_CA.md).

- ...:

  Other arguments passed to the function (currently ignored).

- type:

  type of coefficients,
  `c("fourth_corner", "all_reg", "env2traits_reg", "traits2env_reg")`
  for fourth-corner coefficients and regression coefficients for all
  trait x environmental predictors, environmental predictors only and
  trait predictors only for prediction of the (transformed) response,
  traits and environmental values, respectively.

- rank:

  rank (number of axes to use). Default "full" for all axes (no
  rank-reduction).

- normed:

  logical (default `TRUE`) giving standardized regression coefficients
  and biplot scores. When `FALSE`, (regular) regression coefficients and
  (unstandardized) biplot scores.

## Value

a matrix with coefficients. The exact content of the matrix depends on
the `type` of coefficient that is asked for.

Regression coefficients for a response variable are usually
column-vectors. With **X** the matrix of units-by-predictors and **B**
the matrix of predictors-by-response-variables, predictions or fits are
of the form **Y = XB**. Analogously, `type = "trait2env"` gives a
trait-by-environment matrix and `type = "env2traits"` gives an
environment-by-trait matrix.

## Details

Regression coefficients are for standardized traits and environmental
variables.

With covariates, [`coef()`](https://rdrr.io/r/stats/coef.html) gives
partialfourth-corner correlations. With `rank = 2`,
[`coef()`](https://rdrr.io/r/stats/coef.html) gives the two-dimensional
approximation of the full-rank fourth-corner correlations in the biplot
that displays the traits and environmental variables at arrow heads or
points at `scores(mod, display = c("bp", "bp_traits"))`.

## Examples

``` r
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

mod <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Condition(Manure),
             formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits,
             verbose = FALSE)

# regression coefficients
coef(mod, type = "env2traits")
#>                SLA     Height        LDMC Lifespanannual Lifespanperennial
#> Manure -0.23423414  0.2540139  0.32912396    -0.37130684        0.37130684
#> A1      0.05226663  0.1025887 -0.08560552    -0.12390838        0.12390838
#> Moist  -0.18017363  0.1836804  0.04051036    -0.02158823        0.02158823
#> MagBF  -0.01504873  0.1005833  0.10645047    -0.25809310        0.25809310
#> MagHF  -0.05899194  0.1044629  0.03221358    -0.18057689        0.18057689
#> MagNM  -0.29468992  0.1739784  0.33095579    -0.39679946        0.39679946
#> Use    -0.04563582 -0.1100948 -0.17509044     0.03719925       -0.03719925
#> attr(,"meaning")
#> [1] "standardized regression coefficients to predict traits from environment."
coef(mod, type = "traits2env")
#>                            A1       Moist       MagSF        MagBF        MagHF
#> Seedmass          -0.02240045 -0.04978638 -0.01023869  0.041848273  0.009631558
#> SLA               -0.06845546 -0.16774065 -0.02649254  0.146106902  0.034893995
#> Height             0.15457078  0.15961290  0.02806485 -0.073452767  0.041610702
#> LDMC              -0.17828444 -0.10419687  0.01403797  0.005732252 -0.079309518
#> Lifespanperennial  0.08878173  0.02537622 -0.05136037 -0.005425793  0.007615053
#>                         MagNM          Use
#> Seedmass          -0.03147019  0.004500904
#> SLA               -0.12032477 -0.090142614
#> Height            -0.01014934  0.018558815
#> LDMC               0.05643595 -0.180870869
#> Lifespanperennial  0.04839256  0.080107192
#> attr(,"meaning")
#> [1] "standardized regression coefficients to predict environment from traits."
coef(mod, type = "fourth")
#>                            A1       Moist       MagSF       MagBF        MagHF
#> SLA               -0.07275100 -0.18945586 -0.04438289  0.15899470  0.031696233
#> Height             0.11655218  0.16432389  0.03216548 -0.10690577  0.005095068
#> LDMC              -0.10383747 -0.02955614  0.01837840 -0.03300820 -0.065595295
#> Lifespanannual    -0.06904553  0.00293052  0.05112542 -0.01535635 -0.008559657
#> Lifespanperennial  0.06904553 -0.00293052 -0.05112542  0.01535635  0.008559657
#>                         MagNM          Use
#> SLA               -0.10945491 -0.059594125
#> Height             0.04632044 -0.009109448
#> LDMC               0.06932311 -0.150609601
#> Lifespanannual    -0.03107173 -0.039539841
#> Lifespanperennial  0.03107173  0.039539841
#> attr(,"meaning")
#> [1] "fourth-corner correlation"
coef(mod, type = "all_reg")
#>                        Manure           A1       Moist       MagBF       MagHF
#> Seedmass          -0.12493711 -0.006613744 -0.05931607 -0.04200761 -0.03770288
#> SLA               -0.28292594  0.060443140 -0.15449656 -0.05168545 -0.08321524
#> Height             0.05781765  0.163735690  0.15681800  0.04023132  0.07189288
#> LDMC               0.23683091 -0.153975995 -0.03174838  0.05461067 -0.02638899
#> Lifespanperennial  0.38895095  0.113943476  0.03774092  0.25684954  0.19263442
#>                         MagNM         Use
#> Seedmass          -0.12349465  0.02839218
#> SLA               -0.38998539 -0.08819784
#> Height            -0.07279848 -0.08547934
#> LDMC               0.26477639 -0.16089710
#> Lifespanperennial  0.44747415  0.01244177
#> attr(,"meaning")
#> [1] "standardized regression coefficients for traits and environment to predict the response."
```
