rm(list=ls(all=TRUE))  # remove all existing items from the workspace

resp<- read.csv("data/Wisconsin_response.csv")
names(resp)[1:10]
Y <- resp[,-1]; rownames(resp) <- resp[,"X"]
envir<- read.csv("data/Wisconsin_envir.csv")[,-1]
envir$site<- paste0("site", round(envir$site,0))
traits<- read.csv("data/Wisconsin_traits.csv")
names(envir)
names(traits)
envir$BA <- sqrt(envir$BA)
envir$X.N <- log(envir$X.N)
summary(traits)
nams<- c("SLA","VH","LS","CN","LT") #log transformed
for (k in  nams) traits[[k]] <- log(traits[[k]])
names(traits)[3] <- "LMA"; traits[,3] <- - traits[,3]
names(traits)
keytraits <- c("VH","LS","LMA","LCC")
keyenv <- c("MAT","MAP","TSD","X.Sand","X.N","BA")

formulaEnvkey <- as.formula(paste0("~", paste(keyenv,collapse = "+")))
formulaTraitskey <- as.formula(paste0("~", paste(keytraits,collapse = "+")))


library(douconca)
divide = TRUE
mod1 <- dc_CA(formulaEnv = formulaEnvkey,formulaTraits = formulaTraitskey,
              response = Y, dataEnv = envir, dataTraits = traits, divideBySiteTotals = divide)
set.seed(1457)
anova(mod1,by="axis")
# the default plot
plot(mod1) # gradient_description = "correlation"

plot(mod1, gradient_description = "tvalues")

# nitty-gritty, adapt names, change widths,  and flip the axis -----------------------------------------------------------
# plot with gradient_description = tvalues 

# adapt names
newnams <- getPlotdata(mod1)$newNameList$newnames
newnams$env[4:5] <- c("%Sand", "%N")
plot(mod1, gradient_description = "tvalues", widths = c(3,1,1), newnames = newnams,
      flip_axis  = TRUE, expand =1.5)
