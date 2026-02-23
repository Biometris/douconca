# Print a summary of a wrda or cca0 object

Print a summary of a wrda or cca0 object

## Usage

``` r
# S3 method for class 'wrda'
print(x, ...)
```

## Arguments

- x:

  an object from [`wrda`](wrda.md) or [`cca0`](cca0.md)

- ...:

  Other arguments passed to the function (currently ignored).

## Value

No return value, results are printed to console.

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
#> wRDA      6   36.495 0.62581 3.3449  0.011 *
#> Residual 12   21.822                        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $eigenvalues
#>      wRDA1      wRDA2      wRDA3      wRDA4      wRDA5      wRDA6 
#> 18.1328543  8.2332883  4.0467247  3.1546579  2.2383295  0.6893052 
#> 
```
