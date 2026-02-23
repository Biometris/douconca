# Prediction for double-constrained correspondence analysis (dc-CA)

Prediction of traits from environment, environment from traits and
response from trait and environment data.

## Usage

``` r
# S3 method for class 'dcca'
predict(
  object,
  ...,
  type = c("envFromTraits", "traitsFromEnv", "response", "SNC", "CWM", "lc", "lc_traits"),
  rank = "full",
  newdata = NULL,
  weights = NULL,
  scaling = "symmetric"
)
```

## Arguments

- object:

  return value of [`dc_CA`](dc_CA.md).

- ...:

  Other arguments passed to the function (currently ignored).

- type:

  type of prediction,
  `c("envFromTraits", "traitsFromEnv", "response", "lc", "lc_traits")`
  for environmental values, values of traits, response (expected
  abundance) and constrained scores for sites and species. `"SNC"` is
  equivalent with `"envFromTraits"`. `"CWM"` is equivalent with
  `"traitsFromEnv"`.

- rank:

  rank (number of axes to use). Default "full" for all axes (no
  rank-reduction).

- newdata:

  Data in which to look for variables with which to predict. For
  `type = "envFromTraits" or "traitsFromEnv"` or
  `type = "lc_traits" or "lc"`, `newdata` is a data frame of trait and
  environmental values, respectively, which are used for the prediction
  or the calculation of scores. If omitted, fitted values are generated
  (use [`fitted.dcca`](fitted.dcca.md) instead). For
  `type = "response"`, newdata is a list of two data frames with trait
  and environmental values in this order, *e.g.*
  `list(traits = dataTraits, env = dataEnv)`.

- weights:

  list of weights of species and of sites in `newdata` when
  `type = "response"`, else ignored (default NULL yielding equal species
  and site weights, both summing to 1). Example: weights = list(species
  = c(100, 1, 1), sites = c(1, 1, 1, 1)), in that order, with traits of
  three new species in newdata\[\[1\]\] and environmental values (and
  levels of factors) of four new sites in newdata\[\[2\]\]. Species
  weights are scaled to a sum of one.

- scaling:

  numeric (1,2 or 3) or character `"sites", "species" or "symmetric"`.
  Default: "symmetric". Either site- (1) or species- (2) related scores
  are scaled by eigenvalues, and the other set of scores have unit
  weighted mean square or with 3 both are scaled symmetrically to
  weighted mean squares equal to the square root of eigenvalues.
  Negative values are treated as the corresponding positive ones by
  `abs(scaling)`.

## Value

a matrix with the predictions. The exact content of the matrix depends
on the `type` of predictions that are being made.

## Details

Variables that are in the model but not in `newdata` are set to their
weighted means in the training data. Predictions are thus at the
(weighted) mean of the quantitative variables not included. Predictions
with not-included factors are at the reference level (the first level of
the factor).

For `type = "response"`, many of the predicted values may be negative,
indicating expected absences (0) or small expected response values.

With `type = "traitsFromEnv"` and `newdata = NULL`, predict gives the
fitted mean traits, *i.e.* the fitted community weighted means. With
`type = "envFromTraits"` and `newdata = NULL`, predict gives the fitted
mean environment, *i.e.* the fitted species niche centroids (see
[`fitted.dcca`](fitted.dcca.md)). See `fitted.dcca`.

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

# Ten 'new' sites with a subset of the variables in mod 
# X_lot will be ignored as it is not part of the model
newEnv <- dune_trait_env$envir[1:10, c("A1", "Mag", "Manure", "X_lot")]
newEnv[2,"A1"] <- 3.0
rownames(newEnv) <- paste0("NewSite", 1:10)

pred.traits <- predict(mod, type = "traitsFromEnv", newdata = newEnv)
#> Warning: newdata does not contain the predictor variables
#>  Moist,Use
#> These are set at their mean values and,
#> for factors, at the reference level
#> The current formula is
#>  ~ Manure + A1 + Moist + Mag + Use
head(pred.traits)
#>               SLA   Height     LDMC Lifespanannual Lifespanperennial
#> NewSite1 23.49280 26.38356 266.8780   0.0420223885         0.9579776
#> NewSite2 25.31359 25.43822 255.1297  -0.0008408709         1.0008409
#> NewSite3 23.72696 27.51990 262.5873   0.0250113500         0.9749886
#> NewSite4 23.71135 27.44414 262.8733   0.0261454193         0.9738546
#> NewSite5 25.23184 27.30383 229.8117   0.0211790573         0.9788209
#> NewSite6 24.91963 25.78871 235.5326   0.0438604419         0.9561396


# Eight 'new' species with a subset of traits that are included in the model 
# Variable "L" will be ignored as it is not in the model 
newTraits <- dune_trait_env$traits[1:8, c("Species", "SLA", "LDMC", "L")]
newTraits[3,"SLA"]<- 18 
rownames(newTraits) <- paste("Species",LETTERS[1:8] )# or another meaningful name.

pred.env <- predict(mod, type = "envFromTraits", newdata = newTraits)
#> Warning: newdata does not contain the predictor variables
#>  Height,Lifespan,Seedmass
#> These are set at their mean values and,
#> for factors, at the reference level
#> The current formula is
#>  ~ Seedmass + SLA + Height + LDMC + Lifespan
head(pred.env)
#>                 A1    Moist     MagSF     MagBF     MagHF     MagNM      Use
#> Species A 4.418925 3.099357 0.4193545 0.1162310 0.2567447 0.2076698 1.792234
#> Species B 3.655381 2.386016 0.4099394 0.1990894 0.2308547 0.1601165 1.487752
#> Species C 3.932918 2.921976 0.4313648 0.1056551 0.2053991 0.2575811 1.618013
#> Species D 3.891226 2.356745 0.3966054 0.2306193 0.2703075 0.1024678 1.564131
#> Species E 3.133064 2.374767 0.4345196 0.1466711 0.1530670 0.2657424 1.313215
#> Species F 4.117969 2.529417 0.3968530 0.2149759 0.2829311 0.1052401 1.651777


pred.resp <- predict(mod, type = "response", newdata = list(newTraits, newEnv),
                     weights = list(species = rep(1:2, 4), sites = rep(1, 10)))
#> Warning: newdata does not contain the predictor variables
#>  Height,Lifespan,Seedmass
#> These are set at their mean values and,
#> for factors, at the reference level
#> The current formula is
#>  ~ Seedmass + SLA + Height + LDMC + Lifespan
#> Warning: newdata does not contain the predictor variables
#>  Moist,Use
#> These are set at their mean values and,
#> for factors, at the reference level
#> The current formula is
#>  ~ Manure + A1 + Moist + Mag + Use
colSums(pred.resp) # about alternating 0.8 and 1.6 (reflecting the new species weights)
#> Species A Species B Species C Species D Species E Species F Species G Species H 
#> 0.6754447 1.4335014 0.7063500 1.3996597 0.7529561 1.3738211 0.7113281 1.3499565 
rowSums(pred.resp) # about equal rowsums
#>   NewSite1   NewSite2   NewSite3   NewSite4   NewSite5   NewSite6   NewSite7 
#>  0.9528552  0.1295043  0.5920024  0.6160593  0.7202928  1.2014298  0.1325513 
#>   NewSite8   NewSite9  NewSite10 
#> -0.2042446  2.7755022  1.4870650 
```
