rm(list=ls(all=TRUE))  # remove all existing items from the workspace
# download the spreadsheet from https://doi.org/10.5061/dryad.4qrfj6qb0
dat <- as.data.frame(readxl::read_xlsx("Rolhauser__Waller___Tucker_J_Ecol_-_data.xlsx"))
#View(dat)
names(dat)
envlong <- dat[,c(1,2,13:21)]
names(envlong)
env <- unique(envlong)
dim(env)
# clearly no intra-site variation
names(env)
#View(env)
write.csv(env, "Wisconsin_envir.csv")
spplong<-dat[,c(1:12)]
names(spplong)
#View(spplong)
traitlong <- dat[,c(3, 5:12)]
names(traitlong)
traits<- unique(traitlong)
dim(traits)
# clearly no intra-species variation
write.csv(traits,"Wisconsin_traits.csv")
abunlong <- dat[,c(1,2,3,4)]
names(abunlong)
m <- nrow(traits)
n <- nrow(env)
Y <- matrix(abunlong$freq,nrow = n, ncol = m,byrow = TRUE)
rownames(Y) <- as.character(env$site)
colnames(Y) <- traits$sp
dim(Y)
#View(Y)
write.csv(Y,"Wisconsin_response.csv")
