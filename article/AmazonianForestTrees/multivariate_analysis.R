rm(list=ls(all=TRUE)) 

#Pinho, Bruno; ter Braak, Cajo; P. L. Melo, Felipe;
#Bauman, David; Barlow, Jos (2024). 
#Data and code from Pinho et al. - Winner-Loser plant trait replacements in
#human-modified tropical forests. 
#figshare. Dataset. https://doi.org/10.6084/m9.figshare.25565169

# see README.txt for further use of this data

# Libraries----
library(data.table)
library(tidyr)
library(dplyr)
library(vegan)
library(tibble)

#remotes::install_github("CajoterBraak/douconca") # or:
#install.packages("douconca") # preferred 
library(douconca)#dc-CA

# Data----
comm <- fread("Data/comm.csv")
comm <- dcast.data.table(data = comm, Region_Plot ~ Binomial_correct)
comm <- column_to_rownames(comm, "Region_Plot")

env  <- fread("Data/env.csv") 
env$invBA_ha <- 1/ env$BA_ha

traits <- fread("Data/traits.csv")
traits <- column_to_rownames(traits, "Binomial_correct")

#Create community matrix

# Analysis----

# an initial analysis required 
# all traits are being used in this selection as formulaTraits =~.
names(traits)

out1 <- dc_CA(formulaEnv = ~lnFL2000 + Condition(Region), formulaTraits =~., 
                    response = comm, dataEnv = env, dataTraits = traits)

out1$eigenvalues

## Fitting models according to the DAG 
## (as in Fig.3), but now using dc-CA 

### step 1 Forest loss------------------------------------------------------------------

considered <- NULL
consider <- c("lnFL500","lnFL1000", "lnFL2000")
names(consider) <- consider
consider

fit_measures <- matrix(NA, nrow = length(consider), ncol = 2)
colnames(fit_measures) <- c("variance","pval1")
rownames(fit_measures) <- consider

cntr <- how(within = Within(type = "free"), plots = Plots(strata=  env$Region, type = "none"), nper= 1999)
set.seed(123)

test <- TRUE

for (k in seq_along(consider)){
  
  formulaE_FS <- as.formula(paste("~", consider[k], "+Condition(Region)" ))
  
  out_FS <- dc_CA(formulaE_FS,
                  dc_CA_object = out1, verbose = FALSE)
  
  if (test) {
    an <- anova(out_FS$RDAonEnv, permutations = cntr)
    pval <- an$`Pr(>F)`[1]} else pval <- NA
  fit_measures[k,] <- c(out_FS$inertia["constraintsTE","weighted variance"], pval)
}


pvaladj3 <- c(p.adjust(fit_measures[1:3,"pval1"], method = "holm"), rep(NA, nrow(fit_measures)-3))
#pvaladj <- p.adjust(fit_measures[,"pval1"], method = "holm")
fit_measures <- cbind(fit_measures,pvaladj3)
round(fit_measures, 5) 
# best = lnFL2000 ; all significant with and without correction for multiple testing

fit_measuresL <- list()
fit_measuresL[[1]] <- fit_measures

considered <- c(considered, consider["lnFL2000"])
considered

### step 2 Fragmentation------------------------------------------------------------------
stepk <-2

considerk <-  c("lnNfrag500","lnNfrag1000", "lnNfrag2000")

fit_measures <- matrix(NA, nrow = length(considerk), ncol = 2)
colnames(fit_measures) <- c("variance","pval1")
rownames(fit_measures) <- considerk

for (k in seq_along(considerk)){
  
  formulaE_FS <- as.formula(paste("~", considerk[k], "+Condition(Region +", paste(considered, collapse = "+") ,")", sep ="" ))
  
  out_FS <- dc_CA(formulaE_FS,
                  dc_CA_object = out1, verbose = FALSE)
  
  if (test) {
    an <- anova(out_FS$RDAonEnv, permutations = cntr)
    pval <- an$`Pr(>F)`[1]} else pval <- NA
  fit_measures[k,] <- c(out_FS$inertia["constraintsTE","weighted variance"], pval)
}

pvaladj3 <- c(p.adjust(fit_measures[1:3,"pval1"], method = "holm"), rep(NA, nrow(fit_measures)-3))
#pvaladj <- p.adjust(fit_measures[,"pval1"], method = "holm")
fit_measures <- cbind(fit_measures,pvaladj3 )
round(fit_measures, 5) 
# best = lnNfrag1000 ; not significant, neither with nor without multiple testing correction

fit_measuresL[[stepk]] <- fit_measures

considered <- c(considered, "lnNfrag1000")
considered

### step 3 Edge density------------------------------------------------------------------

stepk <-3
considerk <-  c("ln_edg500" , "ln_edg1000" ,"ln_edg2000")

fit_measures <- matrix(NA, nrow = length(considerk), ncol = 2)
colnames(fit_measures) <- c("variance","pval1")
rownames(fit_measures) <- considerk

for (k in seq_along(considerk)){
  
  formulaE_FS <- as.formula(paste("~", considerk[k], "+Condition(Region +", paste(considered, collapse = "+") ,")", sep ="" ))
  
  out_FS <- dc_CA(formulaE_FS,
                  dc_CA_object = out1, verbose = FALSE)
  
  if (test) {
    an <- anova(out_FS$RDAonEnv, permutations = cntr)
    pval <- an$`Pr(>F)`[1]} else pval <- NA
  fit_measures[k,] <- c(out_FS$inertia["constraintsTE","weighted variance"], pval)
}

pvaladj3 <- c(p.adjust(fit_measures[1:3,"pval1"], method = "holm"), rep(NA, nrow(fit_measures)-3))
fit_measures <- cbind(fit_measures,pvaladj3)
round(fit_measures, 5) #
# best ln_edg500

fit_measuresL[[stepk]] <- fit_measures #
fit_measuresL
considered <- c(considered, "ln_edg500")
considered
### step 4 Local degradation------------------------------------------------------------------

stepk <-4
considerk <-  c("invBA_ha")

fit_measures <- matrix(NA, nrow = length(considerk), ncol = 2)
colnames(fit_measures) <- c("variance","pval1")
rownames(fit_measures) <- considerk

for (k in seq_along(considerk)){
  
  formulaE_FS <- as.formula(paste("~", considerk[k], "+Condition(Region +", paste(considered, collapse = "+") ,")", sep ="" ))
  
  out_FS <- dc_CA(formulaE_FS,
                  dc_CA_object = out1, verbose = FALSE)
  
  if (test) {
    an <- anova(out_FS$RDAonEnv, permutations = cntr)
    pval <- an$`Pr(>F)`[1]} else pval <- NA
  fit_measures[k,] <- c(out_FS$inertia["constraintsTE","weighted variance"], pval)
}

round(fit_measures, 5)

fit_measuresL[[stepk]] <- fit_measures
fit_measuresL


# Final model with lnFL2000 and local degradation -------------------------


outF <- dc_CA(formulaEnv = ~lnFL2000 + invBA_ha + Condition(Region), formulaTraits =~.,
              response = comm, dataEnv = env, dataTraits = traits)

anova(outF, by= "axis", permutations = list(999,cntr))

plot(outF, widths = c(3,1,1.2), remove_centroids = TRUE)


# Extra code showing log(invBA_ha) ----------------------------------------

# avoids an outlier

outF2 <- dc_CA(formulaEnv = ~lnFL2000 + log(invBA_ha) + Condition(Region), formulaTraits =~.,
              response = comm, dataEnv = env, dataTraits = traits)

anova(outF2, by= "axis", permutations = list(999,cntr))


plot(outF2, widths = c(3,1,1.2), remove_centroids = TRUE, flip_axis = TRUE)

