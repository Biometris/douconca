data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
Y <- dune_trait_env$comm[, -1] # must delete "Sites"
Y_N2 <- ipf2N2(Y, updateN2 = FALSE, N2N_N2_species = FALSE)
attr(Y_N2, "iter") # 6

# show that column margins of the transform matrix are
# equal to the Hill N2 values
diff(range(colSums(Y_N2) / apply(X = Y, MARGIN = 2, FUN = fN2))) #  8.881784e-16
diff(range(rowSums(Y_N2) / apply(X = Y, MARGIN = 1, FUN = fN2))) #  0.0201

Y_N2i <- ipf2N2(Y, updateN2 = TRUE, N2N_N2_species = FALSE)
attr(Y_N2i, "iter") # 9
diff(range(colSums(Y_N2i) / apply(X = Y_N2i, MARGIN = 2, FUN = fN2))) #  1.776357e-15
diff(range(rowSums(Y_N2i) / apply(X = Y_N2i, MARGIN = 1, FUN = fN2))) #  0.0229

# the default version:
Y_N2N_N2i <- ipf2N2(Y)
# ie. 
# Y_N2N_N2i <- ipf2N2(Y, updateN2 = TRUE, N2N_N2_species = TRUE)
attr(Y_N2N_N2i, "iter") # 16
N2 <- apply(X = Y_N2N_N2i, MARGIN = 2, FUN = fN2)
N <- ncol(Y)
diff(range(colSums(Y_N2N_N2i) / (N2 * (N - N2)))) #  8.326673e-17

N2 <- apply(X = Y_N2N_N2i, MARGIN = 1, FUN = fN2)
N <- nrow(Y)
diff(range(rowSums(Y_N2N_N2i) / (N2 * (N - N2)))) #  0.207

sum(Y_N2N_N2i) - sum(Y)

mod0 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
              response = Y,  
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits, 
              divide = FALSE,
              verbose = FALSE)

mod1 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
              response = Y_N2N_N2i,  
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits, 
              divide = FALSE,
              verbose = FALSE)

mod1$eigenvalues / mod0$eigenvalues
# as the ratios of eigenvalues are greater than 1,
# all axes have higher (squared) fourth-corner correlation

Y_PA <- 1 * (Y > 0)
Y_PA_N2 <- ipf2N2(Y_PA, N2N_N2_species = FALSE)
attr(Y_PA_N2, "iter") # 1
diff(range(Y_PA - Y_PA_N2)) # 4.440892e-16

Y_PA_N2i <- ipf2N2(Y_PA, N2N_N2_species = TRUE)
attr(Y_PA_N2i, "iter") # 9
N_occ <- colSums(Y_PA) # number of occurrences of species
N <- nrow(Y_PA)
plot(N_occ, colSums(Y_PA_N2i))
cor(colSums(Y_PA_N2i), N_occ * (N - N_occ)) # 0.79
