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

getStart <- function() {
  return(as.numeric(Sys.time()))
}

formatIntTime <- function(i, l) {
  res <- as.character(i)
  remain <- l - nchar(res)
  res <-
    paste(paste(rep("0", remain), collapse = ""), res, sep = "")
  return(res)
}

formatSpendTime <- function(spend) {
  seconds <- as.integer(spend)
  milliseconds <- formatIntTime(as.integer(1000 * (spend %% 1)), 3)
  minutes <- seconds %/% 60
  seconds <- formatIntTime(seconds %% 60, 2)
  hours <- formatIntTime(minutes %/% 60, 2)
  minutes <- formatIntTime(minutes %% 60, 2)
  res <-
    paste(hours, "h", minutes, "m", seconds, ".", milliseconds, "s", sep = "")
  return(res)
}

logWithTime <- function(start, text) {
  current <- getStart()
  spend <- current - start
  print(paste("(", formatSpendTime(spend), ") - ", text, sep = ""))
}

logTotalTime <- function(start) {
  end <- getStart()
  spend <- end - start
  print("##################################################")
  print(paste("#####   Temps écoulé : ", formatSpendTime(spend), sep = ""))
  print("##################################################")
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
  res <- list()
  res$beta1 <- dnorm(theta$beta1, mean = ancestor$beta1)
  res$beta2 <- dnorm(theta$beta2, mean = ancestor$beta2)
  res$delta <-
    dtruncnorm(theta$delta, a = 0, mean = ancestor$delta)
  res$rho <- dtruncnorm(theta$rho,
                        a = -.99,
                        b = .99,
                        mean = ancestor$rho)
  return(res)
}

densitySimpleIidByLog <- function(theta, ancestor) {
  res <- list()
  res$beta1 <-
    dnorm(theta$beta1, mean = ancestor$beta1, log = TRUE)
  res$beta2 <-
    dnorm(theta$beta2, mean = ancestor$beta2, log = TRUE)
  res$delta <-
    log(dtruncnorm(theta$delta, a = 0, mean = ancestor$delta))
  res$rho <- log(dtruncnorm(
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

genSimpleIidIndep <- function(theta) {
  require(truncnorm)
  res <- list()
  res$beta1 <- rnorm(1)
  res$beta2 <- rnorm(1)
  res$delta <- rtruncnorm(1, a = 0)
  res$rho <- rtruncnorm(1,
                        a = -.99,
                        b = .99)
  return(res)
}

densitySimpleIidIndep <- function(theta, ancestor) {
  res <- list()
  res$beta1 <- dnorm(theta$beta1)
  res$beta2 <- dnorm(theta$beta2)
  res$delta <-
    dtruncnorm(theta$delta, a = 0)
  res$rho <- dtruncnorm(theta$rho,
                        a = -.99,
                        b = .99)
  return(res)
}

genNewProposalSimpleIidIndep <- function(theta) {
  res <- list()
  res$value <- genSimpleIidIndep(theta)
  res$density <- densitySimpleIidIndep(res$value, theta)
  res$densityInv <- densitySimpleIidIndep(theta, res$value)
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
  start <- getStart()
  logWithTime(start, "Generation d'un echantillon PMCMC de theta")
  pct <- nb %/% 10
  thetaNameList <- getThetaNameList()
  theta <- theta0
  thetaLkh <-
    estimateLikelyhood(d, theta, Y, N)
  for (i in 1:nb) {
    proposal <- genNewProposal(theta)
    for (name in thetaNameList) {
      newTheta <- theta
      newTheta[[name]] <- proposal$value[[name]]
      newThetaLkh <- estimateLikelyhood(d, newTheta, Y, N)
      lkhRatio <-
        (newThetaLkh * proposal$densityInv[[name]] / (thetaLkh * proposal$density[[name]]))
      proba <- min(lkhRatio, 1)
      # print(proba)
      if (rbinom(1, 1, proba) == 1) {
        theta <- newTheta
        thetaLkh <- newThetaLkh
      }
    }
    res[[i]] <- theta
    if (i %% pct == 0) {
      logWithTime(start, paste((i / pct) * 10, "% généré", sep = ""))
    }
  }
  logTotalTime(start)
  return(res)
}

estimateLogLikelyhood <- function(d, theta, Y, N) {
  return(logLikelihoodBootstrapParticleFilter(d, c(theta$beta1, theta$beta2), theta$delta, theta$rho, Y, N))
}

genThetaPosteriorByLog <-
  function(theta0, d, Y, N, nb, genNewProposal) {
    require(mvtnorm)
    res <- list()
    start <- getStart()
    logWithTime(start, "Generation d'un echantillon PMCMC de theta")
    pct <- nb %/% 10
    thetaNameList <- getThetaNameList()
    theta <- theta0
    thetaLkh <-
      estimateLogLikelyhood(d, theta, Y, N)
    for (i in 1:nb) {
      proposal <- genNewProposal(theta)
      for (name in thetaNameList) {
        newTheta <- theta
        newTheta[[name]] <- proposal$value[[name]]
        newThetaLkh <- estimateLogLikelyhood(d, newTheta, Y, N)
        lkhRatio <-
          newThetaLkh + proposal$densityInv[[name]] - thetaLkh - proposal$density[[name]]
        proba <- min(exp(lkhRatio), 1)
        # print(proba)
        if (rbinom(1, 1, proba) == 1) {
          theta <- newTheta
          thetaLkh <- newThetaLkh
        }
      }
      res[[i]] <- theta
      if (i %% pct == 0) {
        logWithTime(start, paste((i / pct) * 10, "% généré", sep = ""))
      }
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
