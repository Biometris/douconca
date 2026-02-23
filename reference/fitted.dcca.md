# Fitted values of double-constrained correspondence analysis (dc-CA)

Community weighted means (CWM) and species-niche centroids (SNC), as
fitted (in full or reduced rank) from the environmental data and trait
data, respectively, and the fitted response from trait and environment
data.

## Usage

``` r
# S3 method for class 'dcca'
fitted(object, ..., type = c("CWM", "SNC", "response"), rank = "full")
```

## Arguments

- object:

  return value of [`dc_CA`](dc_CA.md).

- ...:

  Other arguments passed to the function (currently ignored).

- type:

  type of prediction, `c( "CWM","SNC", "response")` for environmental
  values, values of traits, response (expected abundance).

- rank:

  rank (number of axes to use). Default "full" for all axes (no
  rank-reduction).

## Value

a matrix with fitted value. The exact content of the matrix depends on
the `type` of fits that are asked for.

## Details

If `type="response"` the rowsums of `object$data$Y` are used to scale
the fit to these sums. Many of the predicted response values may be
negative, indicating expected absences (0) or small expected response
values.

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
# fit the mean traits at each site (20x6),
# that is CWM at each site
CWM.traits <- fitted(mod, type = "CWM")
head(CWM.traits)
#>             SLA   Height     LDMC Lifespanannual Lifespanperennial
#> Site_1 24.69054 23.04477 262.1525    0.047500047         0.9525000
#> Site_2 26.58938 22.47821 248.9740   -0.001033559         1.0010336
#> Site_3 24.27447 25.82005 259.4974    0.028101542         0.9718985
#> Site_4 24.25886 25.74429 259.7835    0.029235611         0.9707644
#> Site_5 26.80648 26.21316 241.2643    0.017241993         0.9827580
#> Site_6 26.11737 22.44992 230.8071    0.049338100         0.9506619

# fit the mean environment for each species (28x8)
# that is SNC of each species
SNC.env <- fitted(mod, type = "SNC")
head(SNC.env)
#>                              A1    Moist     MagSF     MagBF     MagHF
#> Achillea millefolium   5.342999 3.326274 0.2987008 0.1012023 0.2658685
#> Agrostis stolonifera   4.529007 2.576826 0.2877323 0.1864249 0.2361958
#> Aira praecox           3.642275 2.771151 0.4244468 0.1056802 0.1759198
#> Alopecurus geniculatus 4.774155 2.531686 0.2733250 0.2228758 0.2796069
#> Anthoxanthum odoratum  3.917107 2.475008 0.3079576 0.1444949 0.1555197
#> Bellis perennis        4.651703 2.424855 0.2610354 0.2306304 0.2703286
#>                            MagNM      Use
#> Achillea millefolium   0.3342285 2.102427
#> Agrostis stolonifera   0.2896470 1.794929
#> Aira praecox           0.2939532 1.624523
#> Alopecurus geniculatus 0.2241924 1.874519
#> Anthoxanthum odoratum  0.3920278 1.618153
#> Bellis perennis        0.2380056 1.944774

fit.resp <- fitted(mod, type = "response")
# fitted often gives negative values and dc_CA cannot have negatives in the 
# response so, modify fit.resp,
# which gives about similar eigenvalues as the original data
fit.resp[fit.resp < 0] <- 0
mod3 <- dc_CA(formulaEnv = mod$formulaEnv,
              formulaTraits = mod$formulaTraits,
              response = fit.resp, 
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits, 
              verbose = FALSE)
mod3$eigenvalues / mod$eigenvalues
#>     dcCA1     dcCA2     dcCA3     dcCA4 
#> 0.8803851 0.8934459 0.8390431 0.7971850 
```
