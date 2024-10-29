This repository contains code and datasets used in the analyses for the manuscript titled: 
"Winner-loser plant trait replacements in human-modified tropical forests"- Pinho et al. 2024 Nature Ecology & Evolution (accepted)

These data is the result of the work of several people who applied for grants, sampled the data and kept running long-term plots at great expenses. As such, it would be appreciated if data owners were consulted and invited for any papers using this dataset. Here is the list of data owners for each study region:

- Paragominas and Santarém: Erika Berenguer (erikaberenguer@gmail.com) and Jos Barlow (josbarlow@gmail.com)
- Balbina: Maíra Benchimol (mairabs02@gmail.com) and Carlos Peres (c.peres@uea.ac.uk)
- Una: Deborah Faria (deborahuesc@gmail.com) and Eduardo Mariano (marianon@gmail.com)
- Serra Grande: Marcelo Tabarelli (marcelotabarelli9@gmail.com) and Braulio Santos (braulio@dse.ufpb.br)
- Paraíba: Braulio Santos (braulio@dse.ufpb.br), Felipe Melo (fplmelo@gmail.com) and Bruno X. Pinho (bxpinho@hotmail.com)

The data are from 1,207 species occurring across 271 forest plots distributed in six regions of the Amazonian and Atlantic forests in Brazil.
All plots were located in lowland (up to 800 m a.s.l.), old-growth forests within a variable matrix of land uses.

The following datasets are included: species abundance by plot, species functional traits, and environmental (disturbance and climatie) variables.
Landscape-scale variables (forest loss, number of patches, and edge density) were measured at three scales (500, 1,000, 2,000 m), using concentric buffers from the center of each forest plot.

The R scripts in the main folder contain the code related to each component of the data analysis section (see Methods) - i.e. community-level, species-level, and multivariate analyses.
An additional script ('plots.R') reproduces all the graphics in the main and supplementary materials.

The 'comm_analysis.R' file includes (1) the calculation of functional composition - CWM and and upper 5% percentiles of community trait distributions, and
(2) linear models for each region-trait-driver-scale-response combinations, following our causal framework (see Fig. 1) and controlling for spatial autocorrelation in models residuals.

The 'sp_analysis.R' file contains (1) the calculation of observed and null species niche centroids (SNC) and definition of winners and losers from this, and 
(2) linear models of trait-SNC relationships for each trait-region combination, at the largest analyzed landscape scale.
Species-level analyses were performed considering only species occuring in at least five plots within regions, leading to 484 species and 765 species-region combinations.

The 'multivariate_analysis.R' runs the main steps of the dcCA analysis.

Species traits are averages based on multiple records from multiple data sources. Missing values in the species trait dataset were filled by genus-level or imputed data (see Methods)  

Information of species abundance by plot are not public, it must be requested to the authors. 

All R code was run using R version 4.3.2.