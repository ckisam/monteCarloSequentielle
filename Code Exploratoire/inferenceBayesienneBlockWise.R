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
library(truncnorm)

##################################################
##### FONCTIONS UTILITAIRES

getThetaNameList <- function() {
  return(c("beta1", "beta2", "delta", "rho"))
}

##################################################
##### INSTRUMENT : NORMALE MULTIVARIEE (MAUVAIS CHOIX)

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

genNewProposalNormMultivariate <- function(theta) {
  res <- list()
  res$value <- genNormMultivariate(theta)
  res$density <- densityNormMultivariate(res$value, theta)
  res$densityInv <- densityNormMultivariate(theta, res$value)
  return(res)
}

##################################################
##### INSTRUMENT : VARIABLES IID

genSimpleIid <- function(theta) {
  require(truncnorm)
  res <- list()
  res$beta1 <- rnorm(1, mean = theta$beta1)
  res$beta2 <- rnorm(1, mean = theta$beta2)
  res$delta <- rtruncnorm(1, a = 0, mean = theta$delta)
  res$rho <- rtruncnorm(1,
                        a = -.99,
                        b = .99,
                        mean = theta$rho)
  return(res)
}

densitySimpleIid <- function(theta, ancestor) {
  res <- 1
  res <- res * dnorm(theta$beta1, mean = ancestor$beta1)
  res <- res * dnorm(theta$beta2, mean = ancestor$beta2)
  res <- res * dtruncnorm(theta$delta, a = 0, mean = ancestor$delta)
  res <- res * dtruncnorm(theta$rho,
                          a = -.99,
                          b = .99,
                          mean = ancestor$rho)
  return(res)
}

densitySimpleIidByLog <- function(theta, ancestor) {
  res <- 0
  res <- res + dnorm(theta$beta1, mean = ancestor$beta1, log = TRUE)
  res <- res + dnorm(theta$beta2, mean = ancestor$beta2, log = TRUE)
  res <-
    res + log(dtruncnorm(theta$delta, a = 0, mean = ancestor$delta))
  res <- res + log(dtruncnorm(
    theta$rho,
    a = -.99,
    b = .99,
    mean = ancestor$rho
  ))
  return(res)
}

genNewProposalSimpleIid <- function(theta) {
  res <- list()
  res$value <- genSimpleIid(theta)
  res$density <- densitySimpleIid(res$value, theta)
  res$densityInv <- densitySimpleIid(theta, res$value)
  return(res)
}

genNewProposalSimpleIidByLog <- function(theta) {
  res <- list()
  res$value <- genSimpleIid(theta)
  res$density <- densitySimpleIidByLog(res$value, theta)
  res$densityInv <- densitySimpleIidByLog(theta, res$value)
  return(res)
}

##################################################
##### ECHANTILLONNAGE PMCMC DE THETA

estimateLikelyhood <- function(d, theta, Y, N) {
  return(likelihoodBootstrapParticleFilter(d, c(theta$beta1, theta$beta2), theta$delta, theta$rho, Y, N))
}

genThetaPosterior <- function(theta0, d, Y, N, nb, genNewProposal) {
  require(mvtnorm)
  res <- list()
  theta <- theta0
  thetaLkh <-
    estimateLikelyhood(d, theta, Y, N)
  for (i in 1:nb) {
    proposal <- genNewProposal(theta)
    newTheta <- proposal$value
    newThetaLkh <- estimateLikelyhood(d, newTheta, Y, N)
    lkhRatio <-
      (newThetaLkh * proposal$densityInv / (thetaLkh * proposal$density))
    proba <- min(lkhRatio, 1)
    print(proba)
    if (rbinom(1, 1, proba) == 1) {
      theta <- newTheta
      thetaLkh <- newThetaLkh
    }
    res[[i]] <- theta
  }
  return(res)
}

genThetaPosteriorComponentWise <-
  function(theta0, d, Y, N, nb, genNewProposal) {
    require(mvtnorm)
    res <- list()
    theta <- theta0
    thetaLkh <-
      estimateLikelyhood(d, theta, Y, N)
    for (i in 1:nb) {
      proposal <- genNewProposal(theta)
      newTheta <- proposal$value
      newThetaLkh <- estimateLikelyhood(d, newTheta, Y, N)
      lkhRatio <-
        (newThetaLkh * proposal$densityInv / (thetaLkh * proposal$density))
      proba <- min(lkhRatio, 1)
      print(proba)
      if (rbinom(1, 1, proba) == 1) {
        theta <- newTheta
        thetaLkh <- newThetaLkh
      }
      res[[i]] <- theta
    }
    return(res)
  }

estimateLogLikelyhood <- function(d, theta, Y, N) {
  return(logLikelihoodBootstrapParticleFilter(d, c(theta$beta1, theta$beta2), theta$delta, theta$rho, Y, N))
}

genThetaPosteriorByLog <-
  function(theta0, d, Y, N, nb, genNewProposal) {
    require(mvtnorm)
    res <- list()
    theta <- theta0
    thetaLkh <-
      estimateLogLikelyhood(d, theta, Y, N)
    for (i in 1:nb) {
      proposal <- genNewProposal(theta)
      newTheta <- proposal$value
      newThetaLkh <- estimateLogLikelyhood(d, newTheta, Y, N)
      lkhRatio <-
        newThetaLkh + proposal$densityInv - thetaLkh - proposal$density
      proba <- min(exp(lkhRatio), 1)
      print(proba)
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
  thetaNameList <- getThetaNameList()
  for (name in thetaNameList) {
    res[[name]] <- sapply(1:n, function(i) {
      return(rawRes[[i]][[name]])
    })
  }
  return(res)
}
