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
nb <- 100

rho <- 0.5
delta <- 0.2
beta1 <- 0
beta2 <- 0
beta <- c(beta1, beta2)
d <- rep(12, n)

##################################################
##### FONCTIONS DU PROBLEME

getThetaNameList <- function() {
  return(c("beta1", "beta2", "delta", "rho"))
}

genNormMultivariate <- function(theta) {
  require(mvtnorm)
  thetaNameList <- getThetaNameList()
  thetaDim <- length(thetaNameList)
  thetaVect <- sapply(thetaNameList, function(name) {
    return(theta[[name]])
  })
  gen <-
    rmvnorm(
      n = 1,
      mean = thetaVect,
      sigma = diag(x = 1, nrow = thetaDim, ncol = thetaDim)
    )
  res <- list()
  for (i in 1:thetaDim) {
    res[[thetaNameList[i]]] <- gen[i]
  }
  return(res)
}

densityNormMultivariate <- function(theta, mean) {
  require(mvtnorm)
  thetaNameList <- getThetaNameList()
  thetaDim <- length(thetaNameList)
  thetaVect <- sapply(thetaNameList, function(name) {
    return(theta[[name]])
  })
  meanVect <- sapply(thetaNameList, function(name) {
    return(mean[[name]])
  })
  res <-
    dmvnorm(
      x = thetaVect,
      mean = meanVect,
      sigma = diag(x = 1, nrow = thetaDim, ncol = thetaDim)
    )
  return(res)
}

estimateLikelyhood <- function(d, theta, Y, N) {
  return(likelihoodBootstrapParticleFilter(d, c(theta$beta1, theta$beta2), theta$delta, theta$rho, Y, N))
}

genThetaPosterior <- function(theta0, d, Y, N, nb) {
  require(mvtnorm)
  res <- list()
  theta <- theta0
  thetaLkh <-
    estimateLikelyhood(d, theta, Y, N)
  for (i in 1:nb) {
    newTheta <- genNormMultivariate(theta)
    newThetaLkh <- estimateLikelyhood(d, newTheta, Y, N)
    lkhRatio <-
      (newThetaLkh * densityNormMultivariate(theta, newTheta)) / (thetaLkh *
                                                                    densityNormMultivariate(newTheta, theta))
    proba <- min(lkhRatio, 1)
    if (rbinom(1, 1, proba) == 1) {
      theta <- newTheta
      thetaLkh <- newThetaLkh
    }
    res[[i]] <- theta
  }
  return(res)
}

formatResThetaPosterior <- function(rawRes) {
  res <- list()
  n <- length(rawRes)
  thetaNameList <- getThetaNameList
  for (name in thetaNameList) {
    res[[name]] <- sapply(1:n, function(i) {
      return(rawRes[[i]][[name]])
    })
  }
  return(res)
}

##################################################
##### TESTS

### Generation du processus de comptage de photons

simulPhotoCount <- genPhotonCount(d, beta, delta, rho, n)

### Simulation de theta a posteriori

# genThetaPosterior <- function(theta0, d, Y, N, nb) {

theta0 <- list()
theta0$beta1 <- beta1
theta0$beta2 <- beta2
theta0$delta <- delta
theta0$rho <- rho
simulTheta <- genThetaPosterior(theta0, d, simulPhotoCount, N, nb)
