##################################################
##### PROJET HMM ET SMC
##### PMCMC appliqué à l'astrophysique
##### (26/12/2017)
##################################################

rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
setwd(wd)
getwd()


library(truncnorm)
library(mvtnorm)

##################################################
##### FONCTIONS UTILITAIRES

### Calcule un ecart-type pondere
weighted.sd <- function(x, w) {
  return(sqrt(sum(w * (x - mean(
    x
  )) ^ 2) / (length(x) - 1)))
}

getThetaNameList <- function() {
  return(c("beta1", "beta2", "delta", "rho"))
}

### Fonctions de log

getStart <- function(text) {
  return(as.numeric(Sys.time()))
}

logStart <- function(text) {
  print("##################################################")
  print(paste("#####   ", text, sep = ""))
  print("##################################################")
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
##### SIMULATION DU MODELE

### Engendre le processus latent auto-regressif stationnaire Xi
### sur un temps de duree n, avec les parametres delta et rho
genXi <- function(delta, rho, n) {
  sd <- delta / sqrt(1 - (rho ^ 2))
  xi <- rnorm(n = 1, mean = 0, sd = sd)
  res <- c(xi)
  for (i in 2:n) {
    xi <- rnorm(n = 1,
                mean = (rho * xi),
                sd = delta)
    res <- append(res, xi)
  }
  return(res)
}

### Permet de simuler un processus de nombre de photons
### selon le modele propose par Yu et Meng (2011)
genPhotonCount <- function(d, theta, n) {
  res <- list()
  
  ## Extraction des parametres
  beta1 <- theta$beta1
  beta2 <- theta$beta2
  delta <- theta$delta
  rho <- theta$rho
  
  ## Generation du processus xi latent
  xi <- genXi(delta, rho, n)
  res$xi <- xi
  
  ## Generation du processus Y
  res$y <- sapply(1:n, function(i) {
    lambda <- d[i] * exp(beta1 + (beta2 * i / n) + xi[i])
    y <- rpois(n = 1, lambda = lambda)
    return(y)
  })
  
  return(res)
}

### Permet de calculer la vraisemblance des données Y
### lorsque le parametre theta est connu ainsi que le
### processus latent xi.
### - si log=TRUE, la log-vraisemblance est calculee
computeTrueLikelihood <- function(d, theta, xi, Y, log) {
  indLkh <- sapply(1:length(y), function(i) {
    lambda <- d[i] * exp(theta$beta1 + (theta$beta2 * i / n) + xi[i])
    y <- dpois(y, lambda = lambda, log = log)
    return(y)
  })
  if (log) {
    likelihood <- sum(indLkh)
  } else{
    likelihood <- prod(indLkh)
  }
  return(likelihood)
}

##################################################
##### ESTIMATION DE LA VRAISEMBLANCE (theta fixe)
##### PAR FILTRE PARTICULAIRE

### Fonctions de reechantillonnage

### Algorithme de reechantillonnage par defaut
### utilisable dans un filtre particulaire
noResampling <- function(xi, weights) {
  return(xi)
}

### Algorithme de reechantillonnage multinomial
### utilisable dans un filtre particulaire
multinomialResampling <- function(xi, weights, N = length(xi)) {
  res <- sample(
    x = xi,
    size = N,
    replace = TRUE,
    prob = weights
  )
  return(res)
}

### Algorithme de reechantillonnage residuel
### utilisable dans un filtre particulaire
residualResampling <- function(xi, weights) {
  if (sum(weights) == 0) {
    return(xi)
  }
  N <- length(xi)
  w <- weights
  if (sum(w) != 1) {
    w <- w / sum(w)
  }
  
  ## On affecte des copies des particules
  ## proportionnellement a leur poids
  nbCopy <- floor(N * w)
  copy <- do.call(c, lapply(1:N, function(i) {
    return(rep(xi[i], nbCopy[i]))
  }))
  
  ## On tire les particules residuelles par
  ## reechantillonnage multinomial
  w <- N * w - nbCopy
  if (sum(w) > 0) {
    w <- w / sum(w)
  }
  N <- N - length(copy)
  if (N > 0) {
    residual <- multinomialResampling(xi, w, N)
  } else{
    residual <- c()
  }
  
  ## On retourne les particules echanillonnees
  res <- c(copy, residual)
  # sample(res)
  return(res)
}

### Calcule la taille d'echantillonnage effective i.e
### un indicateur permettant de savoir s'il est
### necessaire de reechantillonner pour reduire
### la degenerescence du processus de filtrage
computeEffectiveSampleSize <- function(weights) {
  res <- sum(weights * weights)
  if (res == 0) {
    return(1)
  } else{
    return(1 / res)
  }
}

### Filtre particulaire 'bootstrap'

### Calcule la vraisemblance des donnees observees
### a theta donne, avec le filtre particulaire
### 'bootstrap'.
### - si log=TRUE, la log-vraisemblance est calculee
### - algoResample est une fonction de reechantillonnage
###   prenant en argument les particules et leur poids.
###   L'algorithme par defaut ne reechantillonne pas.
### - si exportXi=TRUE, les particules sont exportees
likelihoodBootstrapParticleFilter <-
  function(d,
           theta,
           Y,
           N,
           log = FALSE,
           algoResample = noResampling,
           exportXi = FALSE) {
    ## Extraction des parametres
    beta1 <- theta$beta1
    beta2 <- theta$beta2
    delta <- theta$delta
    rho <- theta$rho
    
    likelihood <- 1
    if (log) {
      likelihood <- 0
    }
    
    if (exportXi) {
      xiMean <- c()
      xiSd <- c()
      xiSampleSize <- c()
    }
    
    ## t = 1
    # -1- Initialisation
    xiGen <- rnorm(n = N,
                   mean = 0,
                   sd = (abs(delta) / sqrt(1 - rho ^ 2)))
    # -2- Ponderation
    w <- sapply(1:N, function(j) {
      lambda <- d[1] * exp(beta1 + beta2 / n + xiGen[j])
      lkh <- dpois(x = Y[1], lambda = lambda)
      return(lkh)
    })
    if (sum(w) > 0) {
      W <- w / sum(w)
    } else{
      W <- w
    }
    
    if (exportXi) {
      xiMean <- append(xiMean, weighted.mean(xiGen, W))
      xiSd <- append(xiSd, weighted.sd(xiGen, W))
      xiSampleSize <-
        append(xiSampleSize, computeEffectiveSampleSize(W))
    }
    
    if (log) {
      likelihood <- likelihood + log(mean(w))
    } else{
      likelihood <- likelihood * mean(w)
    }
    # -3- Selection
    xiSel <- algoResample(xiGen, W)
    
    ## t > 1
    for (i in 2:length(Y)) {
      # -1- Mutation
      xiGen <- sapply(1:N, function(j) {
        return(rnorm(
          n = 1,
          mean = rho * xiSel[j],
          sd = abs(delta)
        ))
      })
      # -2- Ponderation
      w <- sapply(1:N, function(j) {
        lambda <- d[i] * exp(beta1 + (beta2 * i / n) + xiGen[j])
        lkh <- dpois(x = Y[i], lambda = lambda)
        return(lkh)
      })
      if (sum(w) > 0) {
        W <- w / sum(w)
      } else{
        W <- w
      }
      
      if (exportXi) {
        xiMean <- append(xiMean, weighted.mean(xiGen, W))
        xiSd <- append(xiSd, weighted.sd(xiGen, W))
        xiSampleSize <-
          append(xiSampleSize, computeEffectiveSampleSize(W))
      }
      
      if (log) {
        likelihood <- likelihood + log(mean(w))
      } else{
        likelihood <- likelihood * mean(w)
      }
      # -3- Selection
      xiSel <- algoResample(xiGen, W)
    }
    if (exportXi) {
      res <- list()
      res$likelihood <- likelihood
      res$xi$mean <- xiMean
      res$xi$sd <- xiSd
      res$xi$sampleSize <- xiSampleSize
      return(res)
    } else{
      return(likelihood)
    }
  }

##################################################
##### ECHANTILLONNAGE PMCMC

genThetaPosterior <-
  function(theta0,
           d,
           Y,
           N,
           nb,
           estimateLikelihood = likelihoodBootstrapParticleFilter,
           log = FALSE,
           algoResample = noResampling,
           componentWise = FALSE,
           exportProba = FALSE,
           covariance = diag(rep(1, 4))) {
    res <- list()
    start <- getStart()
    logStart("Generation d'un echantillon PMCMC de theta")
    pct <- nb %/% 10
    thetaNameList <- getThetaNameList()
    theta <- theta0
    thetaLkh <-
      estimateLikelihood(d, theta, Y, N, log, algoResample)
    if (log) {
      thetaLkh <- exp(thetaLkh)
    }
    for (i in 1:nb) {
      if (componentWise) {
        for (j in 1:length(thetaNameList)) {
          param <- thetaNameList[j]
          noise <- rnorm(1, sd = sqrt(covariance[j, j]))
          newTheta <- theta
          newTheta[[param]] <- newTheta[[param]] + noise
          if (abs(newTheta$rho) > .99) {
            if (exportProba) {
              res[[i]]$proba$rho <- proba
            }
            next
          }
          newThetaLkh <-
            estimateLikelihood(d, newTheta, Y, N, log, algoResample)
          if (log) {
            newThetaLkh <- exp(newThetaLkh)
          }
          if (thetaLkh == 0) {
            if (newThetaLkh == 0) {
              proba <- .5
            } else{
              proba <- 1
            }
          } else{
            ratio <- newThetaLkh / thetaLkh
            proba <- min(1, ratio)
          }
          modif <- rbinom(1, 1, proba)
          if (modif == 1) {
            theta <- newTheta
            thetaLkh <- newThetaLkh
          }
          if (exportProba) {
            res[[i]]$proba[[param]] <- proba
          }
        }
      } else{
        noise <- as.vector(rmvnorm(1, rep(0, 4), sigma = covariance))
        newTheta <- theta
        newTheta$beta1 <- newTheta$beta1 + noise[1]
        newTheta$beta2 <- newTheta$beta2 + noise[2]
        newTheta$delta <- newTheta$delta + noise[3]
        newTheta$rho <- newTheta$rho + noise[4]
        if (abs(newTheta$rho) > .99) {
          proba <- 0
        } else{
          newThetaLkh <-
            estimateLikelihood(d, newTheta, Y, N, log, algoResample)
          if (log) {
            newThetaLkh <- exp(newThetaLkh)
          }
          if (thetaLkh == 0) {
            if (newThetaLkh == 0) {
              proba <- .5
            } else{
              proba <- 1
            }
          } else{
            ratio <- newThetaLkh / thetaLkh
            proba <- min(1, ratio)
          }
        }
        modif <- rbinom(1, 1, proba)
        if (modif == 1) {
          theta <- newTheta
          thetaLkh <- newThetaLkh
        }
        if (exportProba) {
          res[[i]]$proba <- proba
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

formatResThetaPosterior <- function(rawRes, exportProba = FALSE) {
  res <- list()
  n <- length(rawRes)
  thetaNameList <- getThetaNameList()
  if (exportProba) {
    thetaNameList <- append(thetaNameList, "proba")
  }
  for (name in thetaNameList) {
    res[[name]] <- do.call(c, lapply(1:n, function(i) {
      return(rawRes[[i]][[name]])
    }))
  }
  return(res)
}

plotSimulResult <-
  function(simulResult,
           burnin,
           pmcmcSize,
           format = TRUE,
           exportProba = FALSE) {
    resultFormat <- simulResult
    if (format) {
      resultFormat <-
        formatResThetaPosterior(resultFormat, exportProba = exportProba)
    }
    par(mfrow = c(2, 2))
    for (param in getThetaNameList()) {
      plot(
        resultFormat[[param]][burnin:pmcmcSize],
        type = "l",
        main = param,
        xlab = "Itération",
        ylab = param
      )
    }
    if (exportProba) {
      for (i in 1:length(getThetaNameList())) {
        param <- getThetaNameList()[i]
        plot(
          resultFormat$proba[i,][burnin:pmcmcSize],
          type = "l",
          main = param,
          xlab = "Itération",
          ylab = param
        )
      }
    }
    par(mfrow = c(1, 1))
  }

printMeanAndVariance <- function(estim, burnin, pmcmcSize) {
  require(mcmcse)
  estimFormat <- formatResThetaPosterior(estim)
  beta1 <- mcse(estimFormat$beta1[burnin:pmcmcSize], method = "bm")
  beta2 <- mcse(estimFormat$beta2[burnin:pmcmcSize], method = "bm")
  delta <- mcse(estimFormat$delta[burnin:pmcmcSize], method = "bm")
  rho <- mcse(abs(estimFormat$rho[burnin:pmcmcSize]), method = "bm")
  print(paste(
    "Beta1 = ",
    beta1,
    " - Beta2 = ",
    beta2,
    " - Delta = ",
    delta,
    " - rho = ",
    rho,
    sep = ""
  ))
}

##################################################
##### SOLUTION DE L'ARTICLE

genThetaGibbsSampler <- function(theta0, d, Y, nb) {
  res <- list()
  beta <- c(theta0$beta1, theta0$beta2)
  delta <-
    for (i in 1:nb) {
      # - 1 - Estimation de Xi
      newXi <- step1GibbsSampler()
      # - 2 - Estimation de Beta
      beta <- step1GibbsSampler()
      # - 3 - Estimation de Rho et Delta
      
    }
  return(res)
}
