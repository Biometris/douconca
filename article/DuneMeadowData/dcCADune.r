rm(list=ls())
library(douconca)
data("dune_trait_env")

abun <-dune_trait_env$comm[, -1]  # must delete column/variable "Sites"
envData <- dune_trait_env$envir
traitData <- dune_trait_env$traits
names(envData)[c(3,4)] <- c("Moisture","Management")

mod_funct_traits <- dc_CA(
  formulaEnv = ~ Moisture + Management,
  formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
  response = abun, envData, traitData)
mod_ecol_traits <- dc_CA(
  formulaEnv = ~ Moisture + Management,
  formulaTraits = ~F+R+N+L,
  response = abun, envData, traitData, verbose = FALSE)
Y <- abun / rowSums(abun)
mod_cca <-cca(Y~Moisture+Management, data = envData)
mod_CCA <- dc_CA(formulaEnv = ~ Moisture+Management,
                 formulaTraits = ~Species,
                 response = abun, envData, traitData, verbose = FALSE)
sum(eigenvals(mod_cca, model = "constrained")) # 0.967
sum(mod_CCA$eigenvalues)                       # 0.967
sum(mod_ecol_traits$eigenvalues)               # 0.567
sum(mod_funct_traits$eigenvalues)              # 0.348
# 
set.seed(1457)
anova(mod_ecol_traits, by = "axis") # axis 1 P =  0.001 axis 2 0.003
anova(mod_funct_traits, by = "axis")# axis 1 P =  0.048 axis 2 0.427
