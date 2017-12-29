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

genNormMultivariate <- function(theta) {
  require(mvtnorm)
  n <- length(theta)
  res <-
    rmvnorm(n = 1,
            mean = theta,
            sigma = diag(x = 1, nrow = n, ncol = n))
  return(res)
}

genThetaPosterior <- function(theta0, d, Y, N, nb) {
  require(mvtnorm)
  theta <- theta0
  beta1 <- c(theta$beta1)
  beta2 <- c(theta$beta2)
  delta <- c(theta$delta)
  rho <- c(theta$rho)
  objLkh <-
    likelihoodBootstrapParticleFilter(d, c(beta1[1], beta2[1]), delta[1], rho[1], Y, N)
  for (i in 1:nb) {
    newTheta <- genNormMultivariate(theta)
    propLkh <- dmvnorm(,sigma = diag(x = 1, nrow = 4, ncol = 4))
  }
}

##################################################
##### TESTS

### Generation du processus de comptage de photons

simulPhotoCount <- genPhotonCount(d, beta, delta, rho, n)
# plot(simulPhotoCount)

### Estimation de la vraisemblance

likelihoodPhotonCount <- sapply(1:10, function(i) {
  return(likelihoodBootstrapParticleFilter(d, beta, delta, rho, simulPhotoCount, N))
})
likelihoodPhotonCount

### Generation d'une normale multivariee
genNormMultivariate(rep(0, 5))
