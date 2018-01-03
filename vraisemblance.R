##################################################
##### PROJET HMM ET SMC
##### Calcul de la vraisemblance
##### (26/12/2017)
##################################################

rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
setwd(wd)
getwd()
source("modele.R")

##################################################
##### FONCTIONS DU PROBLEME

likelihoodBootstrapParticleFilter <-
  function(d, beta, delta, rho, Y, N) {
    ## Premiere etape
    #Initialisation
    xiGen <- rnorm(n = N,
                   mean = 0,
                   sd = (delta / sqrt(1 - rho ^ 2)))
    # Ponderation
    w <- sapply(1:N, function(i) {
      lambda <- d[1] * exp(beta[1] + (beta[2] / n) + xiGen[i])
      lkh <- dpois(x = Y[i], lambda = lambda)
      return(lkh)
    })
    W <- w / sum(w)
    likelihood <- mean(w)
    # Selection
    xiSel <- sample(
      x = xiGen,
      size = N,
      replace = TRUE,
      prob = W
    )
    ## Etapes suivantes
    for (i in 2:length(Y)) {
      # Mutation
      xiGen <- sapply(1:N, function(j) {
        return(rnorm(
          n = 1,
          mean = rho * xiSel[j],
          sd = delta
        ))
      })
      # Ponderation
      w <- sapply(1:N, function(i) {
        lambda <- d[1] * exp(beta[1] + (beta[2] / n) + xiGen[i])
        lkh <- dpois(x = Y[i], lambda = lambda)
        return(lkh)
      })
      W <- w / sum(w)
      likelihood <- likelihood * mean(w)
      # Selection
      xiSel <- sample(
        x = xiGen,
        size = N,
        replace = TRUE,
        prob = W
      )
    }
    return(likelihood)
  }
