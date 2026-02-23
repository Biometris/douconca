# Performs a canonical correspondence analysis

`cca0` is formula-based implementation of canonical correspondence
analysis.

## Usage

``` r
cca0(
  formula,
  response = NULL,
  data,
  traceonly = FALSE,
  cca_object = NULL,
  object4QR = NULL
)
```

## Arguments

- formula:

  one or two-sided formula for the rows (samples) with row predictors in
  `data`. The left hand side of the formula is ignored if it is
  specified in the next argument (`response`). Specify row covariates
  (if any ) by adding `+ Condition(covariate-formula)` to `formula` as
  in [`rda`](https://vegandevs.github.io/vegan/reference/cca.html). The
  `covariate-formula` should not contain a `~` (tilde).

- response:

  matrix or data frame of the abundance data (dimension *n* x *m*).
  Rownames of `response`, if any, are carried through. BEWARE: all rows
  and columns should have positive sums! Can be `NULL` if `cca_object`
  is supplied or if the response is `formula` is two-sided.

- data:

  matrix or data frame of the row predictors, with rows corresponding to
  those in `response` (dimension *n* x *p*).

- traceonly:

  logical, default `FALSE`. If `TRUE`, only the explained variance of
  the predictors and the `Condition()` are returned, *i.e* without
  performing a singular value decompostion.

- cca_object:

  a vegan-type cca-object of *transposed* `response`, from which
  chisq_residuals and row and column weights can be obtained.

- object4QR:

  a vegan-type cca-object with weighted QR's for `formula`, i.e. `qr(Z)`
  and `qr(XZ)` obtainable via `get_QR(object4QR, model = "pCCA")` and
  `get_QR(object4QR, model = "CCA")`, respectively.

## Value

All scores in the `cca0` object are in scaling `"sites"` (1): the
scaling with *Focus on Case distances*.

The returned object has class `c("cca0" "wrda")` so that the methods
`print`, `predict` and `scores` can use the `wrda` variant.

## Details

The algorithm is a wrda on the abundance data after transformation to
chi-square residuals.

It is much slower than
[`cca`](https://vegandevs.github.io/vegan/reference/cca.html). The only
reason to use it, is that [`anova.cca0`](anova.cca0.md) does
residualized predictor permutation. It is unknown to the authors of
`douconca` which method
[`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html)
implements. See [`anova.cca0`](anova.cca0.md).

Compared to
[`cca`](https://vegandevs.github.io/vegan/reference/cca.html), `cca0`
does not have residual axes, *i.e.* no CA of the residuals is performed.

## References

ter Braak C.J.F. and P. Šmilauer (2018). Canoco reference manual and
user's guide: software for ordination (version 5.1x). Microcomputer
Power, Ithaca, USA, 536 pp.

Oksanen, J., et al. (2022) vegan: Community Ecology Package. R package
version 2.6-8. <https://CRAN.R-project.org/package=vegan>.

## See also

[`scores.wrda`](scores.wrda.md), [`anova.cca0`](anova.cca0.md),
[`print.wrda`](print.wrda.md) and [`predict.wrda`](predict.wrda.md)

## Examples

``` r
data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
abun <- dune_trait_env$comm[, -1]  # must delete "Sites"

mod <- cca0(formula = abun ~ A1 + Moist + Mag + Use + Manure,
            data = dune_trait_env$envir)

mod # Proportions equal to those Canoco 5.15
#> 
#> Call: cca0(formula = abun ~ A1 + Moist + Mag + Use + Manure, data =
#> dune_trait_env$envir)
#> 
#>               Inertia Proportion Rank
#> Total          2.1866     1.0000     
#> Constrained    1.2511     0.5722    7
#> Unconstrained  0.9355     0.4278   19
#> 
#> Inertia is scaled Chi-square
#> 
#> Eigenvalues for constrained axes:
#>   CCA1   CCA2   CCA3   CCA4   CCA5   CCA6   CCA7 
#> 0.4759 0.3241 0.1764 0.1313 0.0690 0.0412 0.0333 
#> 
#> Eigenvalues for unconstrained axes:
#> CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 
#>  NA  NA  NA  NA  NA  NA  NA  NA 
#> (Showing 8 of 19 unconstrained eigenvalues)
#> 
#> mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:
#>           Avg    SDS    VIF   Regr1   tval1
#> A1     4.6363 1.8096 1.5022 -0.0962 -1.2237
#> Moist  2.7604 1.7163 1.7279 -0.5275 -6.2574
#> MagBF  0.1773 0.3819 4.7947  0.0107  0.0758
#> MagHF  0.3179 0.4657 2.9965  0.1034  0.9314
#> MagNM  0.2029 0.4021 7.1714 -0.1599 -0.9308
#> Use    1.9042 0.7523 1.6739 -0.1310 -1.5794
#> Manure 1.9521 1.3534 7.8029 -0.0354 -0.1975
#> 

scores(mod, which_cor = c("A1", "X_lot"), display = "cor")
#>             CCA1        CCA2
#> A1    -0.5433670  0.15530408
#> X_lot -0.1808011 -0.08064616
#> attr(,"meaning")
#> [1] "inter set correlation, correlation between environmental variables and the sites scores (CWMs)"

set.seed(123)
anova(mod)
#> $table
#> Permutation test for canonical correspondence analysis
#> Model: cca0(formula = abun ~ A1 + Moist + Mag + Use + Manure, data = dune_trait_env$envir) 
#> Residualized predictor permutation
#> Permutation: free
#> Number of permutations: 999
#> 
#>          df ChiSquare      R2      F Pr(>F)    
#> cca       7   1.25113 0.57218 2.2927  0.001 ***
#> Residual 12   0.93548                          
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $eigenvalues
#>       cca1       cca2       cca3       cca4       cca5       cca6       cca7 
#> 0.47593373 0.32414453 0.17637629 0.13125103 0.06899011 0.04116257 0.03326681 
#> 
anova(mod, by = "axis")
#> $table
#> Permutation test for canonical correspondence analysis
#> Model: cca0(formula = abun ~ A1 + Moist + Mag + Use + Manure, data = dune_trait_env$envir) 
#> Residualized predictor permutation
#> Permutation: free
#> Number of permutations: 999
#> 
#>          df ChiSquare       R2      F Pr(>F)   
#> cca1      1   0.47593 0.217659 6.1051  0.002 **
#> cca2      1   0.32414 0.148241 4.1580  0.018 * 
#> cca3      1   0.17638 0.080662 2.2625  0.236   
#> cca4      1   0.13125 0.060025 1.6836  0.404   
#> cca5      1   0.06899 0.031551 0.8850  0.951   
#> cca6      1   0.04116 0.018825 0.5280  0.996   
#> cca7      1   0.03327 0.015214 0.4267  0.996   
#> Residual 12   0.93548                          
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $eigenvalues
#>       cca1       cca2       cca3       cca4       cca5       cca6       cca7 
#> 0.47593373 0.32414453 0.17637629 0.13125103 0.06899011 0.04116257 0.03326681 
#> 

mod2 <- vegan::cca(abun ~ A1 + Moist + Mag + Use + Manure,
                   data = dune_trait_env$envir)
anova(mod2, by = "axis")
#> Permutation test for cca under reduced model
#> Forward tests for axes
#> Permutation: free
#> Number of permutations: 999
#> 
#> Model: cca(formula = abun ~ A1 + Moist + Mag + Use + Manure, data = dune_trait_env$envir)
#>          Df ChiSquare      F Pr(>F)    
#> CCA1      1   0.47593 6.1051  0.001 ***
#> CCA2      1   0.32414 4.5163  0.002 ** 
#> CCA3      1   0.17638 2.6289  0.080 .  
#> CCA4      1   0.13125 2.0812  0.217    
#> CCA5      1   0.06899 1.1811  0.765    
#> CCA6      1   0.04116 0.7570  0.903    
#> CCA7      1   0.03327 0.6478  0.903    
#> Residual 12   0.93548                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dat <- dune_trait_env$envir
dat$Mag <- "SF"
predict(mod, type = "lc", newdata = dat)
#>               CCA1      CCA2        CCA3        CCA4        CCA5        CCA6
#> Site_1   0.6904890 0.2738131  0.85706483  0.21494213  0.04961136  0.72100459
#> Site_2   0.7086329 0.2734695  0.33481865 -1.07270103 -1.35928061  0.53160107
#> Site_3   0.2244798 0.5945085  0.60570803  0.44083874  0.37245365  0.46727185
#> Site_4   0.2308792 0.5789014  0.63318079  0.42798545  0.35582702  0.47363638
#> Site_5   0.7391771 0.6853586 -0.03558872 -0.62800472 -0.21511036 -0.01408255
#> Site_6   0.6574380 0.3983263  0.11503656 -0.96987465 -1.22626765  0.48068489
#> Site_7   0.5122323 0.2441261  0.29326645 -0.55866996 -1.39165267  1.01605565
#> Site_8  -1.0574335 0.8089811  0.55158630 -0.24633521 -0.86510811  0.29389272
#> Site_9  -0.1730258 0.4845439  0.99593856 -1.55170714 -1.18971277 -0.39582562
#> Site_10  0.5926090 0.2489378  0.78436035 -1.66931462 -1.40310518 -0.05383772
#> Site_11  0.5303763 0.2437824 -0.22897973 -1.84631313 -2.80054465  0.82665213
#> Site_12 -0.5486085 0.8921994  0.18514902 -0.67778374 -0.75653943 -0.08957766
#> Site_13 -0.9628957 1.0647992  0.45590655  0.06982812  0.11279594 -0.18814568
#> Site_14 -1.2893919 1.4405532 -1.34442983 -1.65724140 -2.30506793 -0.24797572
#> Site_15 -1.2204515 1.7587997 -1.55000063 -1.28966488 -1.26065740 -0.75547220
#> Site_16 -1.1534240 1.0430876  0.13949487 -0.05353574 -0.61570879  0.19842488
#> Site_17  0.5792832 0.3033908  0.42708259 -2.26814966 -2.04935800 -0.17081531
#> Site_18  0.9109057 0.3104445  0.10151139 -2.22412702 -2.02304124 -0.05073754
#> Site_19 -0.5115748 0.5163361  0.99170477 -2.20741811 -1.87890897 -0.62651646
#> Site_20 -0.4987760 0.4851219  1.04665030 -2.23312471 -1.91216221 -0.61378742
#>                CCA7
#> Site_1  -0.21001080
#> Site_2  -0.96564817
#> Site_3  -0.19205750
#> Site_4  -0.18164028
#> Site_5  -1.37467944
#> Site_6  -1.04898590
#> Site_7  -0.43402001
#> Site_8   0.11698515
#> Site_9  -0.92255561
#> Site_10 -1.22930984
#> Site_11 -1.18965737
#> Site_12 -0.68260950
#> Site_13 -0.18787396
#> Site_14 -1.43836815
#> Site_15 -1.78489613
#> Site_16 -0.03927309
#> Site_17 -1.64358878
#> Site_18 -1.88030362
#> Site_19 -1.08970250
#> Site_20 -1.06886806
```
