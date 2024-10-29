data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# use vegan::rda in step 2
divide <- TRUE # divide by site.totals if TRUE

Y <- dune_trait_env$comm[, -1]  # must delete "Sites"
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
traits <- dune_trait_env$traits
envir <- dune_trait_env$envir

# test a normal model
mod1a <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
               formulaTraits = ~ SLA + Height + LDMC + Seedmass +Lifespan,
               response = Y,
               dataEnv = envir,
               dataTraits = traits,
               divideBySiteTotals = divide, 
               verbose = FALSE)

# five types of full rank prediction:
predDivT1a_reg <- coef(mod1a, type = "env2traits")
predDivT1a_reg_traits <- coef(mod1a, type = "traits2env")
predDivT1a_env <- predict(mod1a, type = "env")
predDivT1a_traits <- predict(mod1a, type = "traits")
predDivT1a_response <- predict(mod1a, type = "response")

expect_equal_to_reference(predDivT1a_reg, "predDivT1a_reg") 
expect_equal_to_reference(predDivT1a_reg_traits, "predDivT1a_reg_traits") 
expect_equal_to_reference(predDivT1a_env, "predDivT1a_env") 
expect_equal_to_reference(predDivT1a_traits, "predDivT1a_traits") 
expect_equal_to_reference(predDivT1a_response, "predDivT1a_response") 

# five type of rank 1 predDivTiction
predDivT1a_reg1 <- coef(mod1a, type = "env2traits", rank = 1)
predDivT1a_reg_traits1 <- coef(mod1a, type = "traits2env", rank = 1)
predDivT1a_env1 <- predict(mod1a, type = "env", rank = 1)
predDivT1a_traits1 <- predict(mod1a, type = "traits", rank = 1)
predDivT1a_response1 <- predict(mod1a, type = "response", rank = 1)

expect_equal_to_reference(predDivT1a_reg1, "predDivT1a_reg1") 
expect_equal_to_reference(predDivT1a_reg_traits1, "predDivT1a_reg_traits1") 
expect_equal_to_reference(predDivT1a_env1, "predDivT1a_env1") 
expect_equal_to_reference(predDivT1a_traits1, "predDivT1a_traits1") 
expect_equal_to_reference(predDivT1a_response1, "predDivT1a_response1") 

expect_equal(predict(mod1a, type = "env", rank = 1),
             fitted(mod1a, type = "SNC", rank = 1))
expect_equal(predict(mod1a, type = "traits", rank = 1),
             fitted(mod1a, type = "CWM", rank = 1))
expect_equal(predict(mod1a, type = "response", rank = 1, weights = mod1a$weights) *
               sum(mod1a$data$Y), 
             fitted(mod1a, type = "response", rank = 1))

# test factors only
envir$fUse <- factor(envir$Use)
traits$fSLA <- cut(traits$SLA, 3)
mod_fact<- dc_CA(formulaEnv = ~ fUse + Mag,
                 formulaTraits = ~ fSLA + Lifespan ,
                 response = Y, 
                 dataEnv = envir,
                 dataTraits = traits,
                 divideBySiteTotals = divide,
                 verbose = FALSE)

expect_equal(colnames(predict(mod_fact, type = "env")), 
             c("fUse1", "fUse2", "fUse3", "MagSF", "MagBF", "MagHF", "MagNM"))

expect_equal(colnames(predict(mod_fact, type = "traits")), 
             c( "fSLA(5.77,15]", "fSLA(15,24.2]", "fSLA(24.2,33.4]", 
                "Lifespanannual", "Lifespanperennial"))

traits$lifequant <- 1 * traits$Lifespan %in% "perennial"
mod_factb <- dc_CA(formulaEnv = ~ fUse + Mag,
                   formulaTraits = ~ fSLA + lifequant,
                   response = Y, 
                   dataEnv = envir,
                   dataTraits = traits,
                   divideBySiteTotals = divide,
                   verbose = FALSE)

expect_silent(mod_fact)
expect_equal(mod_fact$eigenvalues, mod_factb$eigenvalues)
expect_equivalent(mod_fact$c_traits_normed0, mod_factb$c_traits_normed0)

# test quant variable only, 1 quant trait

modq11<- dc_CA(formulaEnv = ~ Manure,
               formulaTraits = ~ SLA,
               response = Y, 
               dataEnv = envir, 
               dataTraits = traits,
               divideBySiteTotals = divide,
               verbose = FALSE)

modq11_predDivTenv <- predict(modq11, type = "env")
expect_equal_to_reference(modq11_predDivTenv, "modq11_predDivTenv")

modq11_predDivTtraits <- predict(modq11, type = "traits")
expect_equal_to_reference(modq11_predDivTtraits, "modq11_predDivTtraits")

# check how douconca manages collinear models A
envir$Sites <- factor(dune_trait_env$envir$Sites)
mod_dccaFF <- dc_CA(formulaEnv = ~ Sites,
                    formulaTraits = ~ Species,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv = envir,
                    dataTraits = traits,
                    divideBySiteTotals = divide,
                    verbose = FALSE)

Yfit <- fitted(mod_dccaFF, type = "response")
expect_equivalent(Yfit, as.matrix(mod_dccaFF$data$Y)) # full rank fit!

# check how douconca manages collinear models

#  A11 collinear variable
envir$Sites <- factor(dune_trait_env$envir$Sites)
envir$A11 <- envir$A1

mod_DivT_dccaA1 <- dc_CA(formulaEnv = ~ A1 + Manure + Moist,
                         formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                         response = Y, 
                         dataEnv = envir,
                         dataTraits = traits,
                         divideBySiteTotals = divide,
                         verbose = FALSE)
expect_stdout( mod_DivT_dccaA11 <- 
                 dc_CA(formulaEnv = ~ A1 + A11 + Manure + Moist,
                       formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                       response = Y, 
                       dataEnv = envir,
                       dataTraits = traits,
                       divideBySiteTotals = divide,
                       verbose = FALSE))

expect_equal(predict(mod_DivT_dccaA1, type = "response"),
             predict(mod_DivT_dccaA11, type = "response"))

# SLA11 collinear
traits$SLA11 <- traits$SLA
mod_DivT_dccaSLA <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                          response = Y, 
                          dataEnv = envir, 
                          dataTraits = traits,
                          divideBySiteTotals = divide,
                          verbose = FALSE)

expect_message(
  mod_DivT_dccaSLA11 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                              formulaTraits = ~ SLA + SLA11 + Height + LDMC + Seedmass + Lifespan,
                              response = Y, 
                              dataEnv = envir,
                              dataTraits = traits,
                              divideBySiteTotals = divide,
                              verbose = FALSE))

expect_equal(predict(mod_DivT_dccaSLA, type = "response"),
             predict(mod_DivT_dccaSLA11, type = "response"))
