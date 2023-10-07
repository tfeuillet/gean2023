##################################################
############## GWR-ML hybrid model ###############
################## October 2023 ##################
##### Feuillet et al., Geographical analysis #####
########## thierry.feuillet@unicaen.fr ###########
##################################################

## Warning: since our analyses involved personal health data, the code
## cannot be reproducible. You have to put your own data instead of "data"
## in the following script


## Load packages
library(sf)
library(lme4)
library(GWmodel)
library(spdep)
library(dplyr)
library(ClustGeo)
library(GenSA)

## Import your spatial data
data <- "Your own level-1 data" # Not reproducible
dataSp <- data %>% as_Spatial() # Needed for GWR

## Fix parameters for subsequent analyses
formula_gwr <- "Your GWR formula"
formula_ml <- "Your multilevel model formula" # Your formula needs including '|clust' somewhere (i.e. random intercept and/or slopes by clust, clust = GWR beta-based clusters further defined)
matDist <- gw.dist(dp.locat = coordinates(dataSp)) # Distance matrix for GWR

## GWR-ML hybrid model AIC function to be further minimized
func <-function(par){
  
  # GWR model
  GWR <- gwr.basic(data=dataSp, bw=par[1], kernel="bisquare", 
                    adaptive=T, dMat=matDistTmp, formula = formula_gwr) # Choose your own spatial weighting scheme
  GWR1betas <- GWR$SDF[,c(1:X)] # X = number of predictors including intercept
  
  # Regionalization of GWR betas with clustGeo package
  D0 <- dist(GWR1betas@data) # Dissimilarity matrix based on GWR betas
  D1 <- as.dist(matDist) # Distance matrix
  tree <- hclustgeo(D0, D1, alpha = 0.75) # Choose your own spatial constraint (alpha)
  Pk <- cutree(tree, par[2]) # cut the dendrogram to get the partition in k clusters
  data$clust <- Pk %>% as.factor()
  
  # Multilevel model
  lme <- lmer(formula_ml, data = data)
  
  # AIC outcome
  AIC(lme)
}

## Optimization for finding out the optimal GWR bandwith (bw) and number of clusters (k)
## based on the generalized simulated annealing algorithm

optim <- GenSA::GenSA(par = c(50,31), fn = func, lower = c(50,30), upper = c(200,80),
                      control = list(max.time=30000)) # Can take looooong time :)
optim$par # Optimal bw + k

## Final GWR
GWR_fin <- gwr.basic(data=dataSp, bw=optim$par[1], kernel="bisquare", adaptive=T, dMat=matDist, formula = formula_gwr)

## Test for spatial nonstationnarity with Leung tests (Leung et al., 2000)
GWR_leung <- gwr.basic(data = dataSp, bw = optim$par[1], kernel = "bisquare", 
                       adaptive = T, dMat = matDist, 
                       formula = formula_gwr, F123.test = TRUE)
GWR_leung$Ftests

## GWR betas regionalization

GWR1betas <- GWR_fin$SDF[,c(1:X)] # X being the number of predictors including intercept

## clustering with ClustGeo
D <- dist(GWR_fin$SDF[,1])
D0 <- dist(GWR1betas@data) # Beta dissimilarity matrix
D1 <- as.dist(matDist) # Distance matrix
tree <- hclustgeo(D, D1, alpha = 0.7)
k <- optim$par[2]
Pk <- cutree(tree, k) # Cut the dendrogram to get the partition in k clusters
data$clust <- Pk %>% as.factor()

## Multilevel model with GWR beta regions as random effect
lme <- lmer(formula_ml, data = data) # # Your formula needs including '|clust' somewhere (i.e. random intercept and/or slopes)

# Enjoy the results :)


