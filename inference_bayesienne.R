##################################################
##### PROJET HMM ET SMC
##### Inference bayesienne par PMCMC
##### (28/12/2017)
##################################################

rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
setwd(wd)
getwd()
source("vraisemblance.R")

##################################################
##### PARAMETRES DU MODELE

n <- 100
N <- 100

rho <- 0.5
delta <- 0.2
beta1 <- 0
beta2 <- 0
beta <- c(beta1, beta2)
d <- rep(12, n)

##################################################
##### FONCTIONS DU PROBLEME



##################################################
##### TESTS

### Generation du processus de comptage de photons

simulPhotoCount <- genPhotonCount(d, beta, delta, rho, n)
plot(simulPhotoCount)

### Estimation de la vraisemblance

likelihoodPhotonCount <- sapply(1:10, function(i) {
  return(likelihoodBootstrapParticleFilter(d, beta, delta, rho, simulPhotoCount, N))
})
likelihoodPhotonCount
