##################################################
##### PROJET HMM ET SMC
##### Calcul de la vraisemblance
##### (26/12/2017)
##################################################

rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
setwd(wd)
getwd()

##################################################
##### PARAMETRES DU MODELE

n <- 10000

rho <- 0.5
delta <- 0.2
beta1 <- 0
beta2 <- 0
beta <- c(beta1, beta2)
d <- rep(12, n)

##################################################
##### FONCTIONS DU MODELE
