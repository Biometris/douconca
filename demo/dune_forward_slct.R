data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

mod <- dc_CA(formulaEnv = ~ A1,
             formulaTraits = ~ F + R + N + L,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir[, -c(1, 7:10)], # without Sites, X, Y etc.
             dataTraits = dune_trait_env$traits,
             verbose = FALSE)

out1 <- dc_CAForwardSelectionEnv(mod, verbose = FALSE) 
names(out1)
out1$process
out1$finalWithOneExtra
out1$model_final

mod2 <- dc_CA(formulaEnv = ~ Moist + Mag,
              formulaTraits = ~ F,
              response = dune_trait_env$comm[, -1],  # must delete "Sites"
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits[, -c(1, 2)],
              verbose = FALSE)

out2 <- dc_CAForwardSelectionTraits(mod2, verbose = FALSE) 
names(out2)
out2$process
out2$finalWithOneExtra
out2$model_final
