# The package douconca performs double constrained correspondence analysis for trait-environment analysis in ecology

Double constrained correspondence analysis (dc-CA) for analyzing
(multi-)trait (multi-)environment ecological data using library `vegan`
and native R code. It has a `formula` interface which allows to assess,
for example, the importance of trait interactions in shaping ecological
communities. The function `dc_CA` has an option to divide the abundance
data of a site by the site total, giving equal site weights. This
division has the advantage that the multivariate analysis corresponds
with an unweighted (multi-trait) community-level analysis, instead of
being weighted.

Throughout the two step algorithm of ter Braak et al. (2018) is used.
This algorithm combines and extends community- (sample-) and
species-level analyses, *i.e.* (1) the usual community weighted means
(CWM) regression analysis and (2) the species-level analysis of
species-niche centroids (SNC) regression analysis. The SNC is the center
of the realized niche of the species along an environmental variable or,
in the case of dc-CA, an environmental gradient, *i.e.* the dc-CA
ordination axis. Computationally, dc-CA can be carried out by a single
singular value decomposition (ter Braak et al. 2018), but it is here
computed in two steps.

The first step uses canonical correspondence analysis
([`cca`](https://vegandevs.github.io/vegan/reference/cca.html)) to
regress the (transposed) abundance data on to the traits and the second
step uses weighed redundancy analysis ([`wrda`](wrda.md) or, with equal
site weights,
[`rda`](https://vegandevs.github.io/vegan/reference/cca.html)) to
regress the CWMs of the orthonormalized traits, obtained from the first
step, on to the environmental predictors. The second step is thus a
community-level analysis.

If `divideBySiteTotals = FALSE`, the second step uses [`wrda`](wrda.md)
and performs a weighted redundancy analysis of the CWMs on to the
environmental variables.

Division of the abundance data by the site totals has the advantage that
the resulting analysis (without dimension reduction, *i.e.* retaining
all dc-CA axes) corresponds with a series of unweighted community-level
analyses, instead of the analyses being weighted.

Warning: The `dcCA` package was built from `vegan` version 2.6-4 and
uses some of the internal structure of the `vegan`
[`cca.object`](https://vegandevs.github.io/vegan/reference/cca.object.html)
in the not-exported functions `f_inertia` and `get_QR` in the source
code file `functions_using_cca_object_internals.r`.

The main user-functions are [`dc_CA`](dc_CA.md),
[`plot.dcca`](plot.dcca.md), [`scores.dcca`](scores.dcca.md),
[`print.dcca`](print.dcca.md) and [`anova.dcca`](anova.dcca.md).

## References

ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
double constrained correspondence analysis. Environmental and Ecological
Statistics, 25(2), 171-197.
[doi:10.1007/s10651-017-0395-x](https://doi.org/10.1007/s10651-017-0395-x)

Oksanen, J., et al. (2022) vegan: Community Ecology Package. R package
version 2.6-4. <https://CRAN.R-project.org/package=vegan>.

## See also

[`cca`](https://vegandevs.github.io/vegan/reference/cca.html) and
[`rda`](https://vegandevs.github.io/vegan/reference/cca.html)

## Author

**Maintainer**: Bart-Jan van Rossum <bart-jan.vanrossum@wur.nl>
([ORCID](https://orcid.org/0000-0002-8673-2514))

Authors:

- Cajo J.F ter Braak <cajo.terbraak@wur.nl>
  ([ORCID](https://orcid.org/0000-0002-0414-8745))
