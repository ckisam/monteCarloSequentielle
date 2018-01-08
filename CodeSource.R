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
  w <- w / sum(w)
  N <- N - length(copy)
  residual <- multinomialResampling(xi, w, N)
  
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
  res <- 1 / (weights * weights)
  res <- sum(res)
  return(res)
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
    }
    
    ## t = 1
    # -1- Initialisation
    xiGen <- rnorm(n = N,
                   mean = 0,
                   sd = (delta / sqrt(1 - rho ^ 2)))
    # -2- Ponderation
    w <- sapply(1:N, function(j) {
      lambda <- d[1] * exp(beta1 + (beta2 / n) + xiGen[j])
      lkh <- dpois(x = Y[1], lambda = lambda)
      return(lkh)
    })
    W <- w / sum(w)
    
    if (exportXi) {
      xiMean <- append(xiMean, weighted.mean(xiGen, W))
      xiSd <- append(xiSd, weighted.sd(xiGen, W))
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
      # print(computeEffectiveSampleSize(W) / N)
      
      # -1- Mutation
      xiGen <- sapply(1:N, function(j) {
        return(rnorm(
          n = 1,
          mean = rho * xiSel[j],
          sd = delta
        ))
      })
      # -2- Ponderation
      w <- sapply(1:N, function(j) {
        lambda <- d[1] * exp(beta1 + (beta2 / n) + xiGen[j])
        lkh <- dpois(x = Y[i], lambda = lambda)
        return(lkh)
      })
      W <- w / sum(w)
      
      if (exportXi) {
        xiMean <- append(xiMean, weighted.mean(xiGen, W))
        xiSd <- append(xiSd, weighted.sd(xiGen, W))
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
      return(res)
    } else{
      return(likelihood)
    }
  }

##################################################
##### ECHANTILLONNAGE PMCMC

### Proposal independante et autocorrelee

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

densitySimpleIid <- function(theta, ancestor, log = FALSE) {
  require(truncnorm)
  res <- list()
  res$beta1 <- dnorm(theta$beta1, mean = ancestor$beta1, log = log)
  res$beta2 <- dnorm(theta$beta2, mean = ancestor$beta2)
  res$delta <-
    dtruncnorm(theta$delta, a = 0, mean = ancestor$delta)
  if (log) {
    res$delta <- log(res$delta)
  }
  res$rho <- dtruncnorm(theta$rho,
                        a = -.99,
                        b = .99,
                        mean = ancestor$rho)
  if (log) {
    res$rho <- log(res$rho)
  }
  return(res)
}

genNewProposalSimpleIid <- function(theta) {
  res <- list()
  res$value <- genSimpleIid(theta)
  res$density <- densitySimpleIid(res$value, theta)
  res$densityInv <- densitySimpleIid(theta, res$value)
  return(res)
}

### Proposal independante non autocorrelee

genSimpleIidIndep <- function(theta) {
  require(truncnorm)
  res <- list()
  res$beta1 <- rnorm(1, mean = 0)
  res$beta2 <- rnorm(1, mean = 0)
  res$delta <- rtruncnorm(1, a = 0, mean = 0)
  res$rho <- rtruncnorm(1,
                        a = -.99,
                        b = .99,
                        mean = 0)
  return(res)
}

densitySimpleIidIndep <- function(theta, log = FALSE) {
  require(truncnorm)
  res <- list()
  res$beta1 <- dnorm(theta$beta1, mean = 0, log = log)
  res$beta2 <- dnorm(theta$beta2, mean = 0)
  res$delta <-
    dtruncnorm(theta$delta, a = 0, mean = 0)
  if (log) {
    res$delta <- log(res$delta)
  }
  res$rho <- dtruncnorm(theta$rho,
                        a = -.99,
                        b = .99,
                        mean = 0)
  if (log) {
    res$rho <- log(res$rho)
  }
  return(res)
}

genNewProposalSimpleIidIndep <- function(theta) {
  res <- list()
  res$value <- genSimpleIidIndep(theta)
  res$density <- densitySimpleIidIndep(res$value)
  res$densityInv <- densitySimpleIidIndep(theta)
  return(res)
}

genThetaPosterior <-
  function(theta0,
           d,
           Y,
           N,
           nb,
           estimateLikelihood = likelihoodBootstrapParticleFilter,
           log = FALSE,
           algoResample = noResampling,
           genNewProposal = genNewProposalSimpleIid) {
    res <- list()
    start <- getStart()
    logStart("Generation d'un echantillon PMCMC de theta")
    pct <- nb %/% 10
    thetaNameList <- getThetaNameList()
    theta <- theta0
    thetaLkh <-
      estimateLikelihood(d, theta, Y, N, log, algoResample)
    for (i in 1:nb) {
      proposal <- genNewProposal(theta)
      for (name in thetaNameList) {
        newTheta <- theta
        newTheta[[name]] <- proposal$value[[name]]
        newThetaLkh <-
          estimateLikelihood(d, newTheta, Y, N, log, algoResample)
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