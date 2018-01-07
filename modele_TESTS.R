##################################################
##### PROJET HMM ET SMC
##### Generation de la serie
##### --- TESTS
##### (26/12/2017)
##################################################

rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
setwd(wd)
getwd()
source("modele.R")

##################################################
##### PARAMETRES DU MODELE

n <- 1000

rho <- 0.5
delta <- 0.2
beta1 <- 0
beta2 <- 0
beta <- c(beta1, beta2)
d <- rep(12, n)

##################################################
##### TESTS

### Generation du processus 'Xi'

simulXi <- genXi(delta, rho, n)

plot(simulXi,
     main = "Processus 'Xi' (latent) auto-régressif stationnaire",
     xlab = "Temps t",
     ylab = "Valeur de Xi observée")

### Generation du processus de comptage de photons

simulPhotoCount <- genPhotonCount(d, beta, delta, rho, n)

plot(simulPhotoCount)
