##################################################
##### PROJET HMM ET SMC
##### PMCMC appliqué à l'astrophysique
#####
##### NB : le code source des fonctions se trouve
#####    dans le fichier R 'CodeSource'
##### (26/12/2017)
##################################################

rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
setwd(wd)
getwd()
source("CodeSource.R")

##################################################
##### SIMULATION DU MODELE

### Choix des paramètres

n <- 100

theta <- list()
theta$rho <- 0.5
theta$delta <- 0.2
theta$beta1 <- 0
theta$beta2 <- 0
d <- rep(12, n)

### Generation du processus de comptage de photons

simul <- genPhotonCount(d, theta, n)
simulXi <- simul$xi
plot(simulXi,
     main = "Processus 'Xi' (latent) auto-régressif stationnaire",
     xlab = "Temps t",
     ylab = "Valeur de Xi observée")
simulY <- simul$y
plot(simulY,
     main = "Nombre de photons Y au cours du temps",
     xlab = "Temps t",
     ylab = "Y")

##################################################
##### FILTRAGE PARTICULAIRE - CALCUL DE LA VRAISEMBLANCE

N <- 10

resLkh100 <-
  likelihoodBootstrapParticleFilter(d, theta, simulY, N, exportXi = TRUE)

plotPartAndCi <- function(resLkhFun, nbPart) {
  plot(resLkhFun$xi$mean,
       main = "Processus Xi estimé",
       xlab = "Temps t",
       ylab = "Xi")
  lines(simulXi)
  lines(resLkhFun$xi$mean - (qt(.975, df = nbPart) * resLkhFun$xi$sd / sqrt(nbPart)), col =
          "red")
  lines(resLkhFun$xi$mean + (qt(.975, df = nbPart) * resLkhFun$xi$sd / sqrt(nbPart)), col =
          "red")
}

plotPartAndCi(resLkh100, N)

# Avec reechantillonnage residuel
resLkhRes100 <-
  likelihoodBootstrapParticleFilter(d,
                                    theta,
                                    simulY,
                                    N,
                                    algoResample = residualResampling,
                                    exportXi = TRUE)

plotPartAndCi(resLkhRes100, N)

N <- 50

resLkh1000 <-
  likelihoodBootstrapParticleFilter(d, theta, simulY, N, exportXi = TRUE)

plotPartAndCi(resLkh1000, N)

# Avec reechantillonnage residuel
resLkhRes1000 <-
  likelihoodBootstrapParticleFilter(d,
                                    theta,
                                    simulY,
                                    N,
                                    algoResample = residualResampling,
                                    exportXi = TRUE)

plotPartAndCi(resLkhRes1000, N)

### Nombre optimal de particules

drawOptimal <- FALSE

if (drawOptimal) {
  essByParticles <-
    function(d, theta, n, testZone, algoResample = noResampling) {
      start <- getStart()
      logStart(paste("Optimisation de l'ESS pour n=", n, sep = ""))
      d <- rep(d, n)
      Y <- genPhotonCount(d, theta, n)$y
      effSize <- sapply(testZone, function(i) {
        res <- likelihoodBootstrapParticleFilter(d,
                                                 theta,
                                                 Y,
                                                 i,
                                                 algoResample = algoResample,
                                                 exportXi = TRUE)
        res <- res$xi$sampleSize
        return(mean(res))
      })
      logTotalTime(start)
      return(effSize)
    }
  testZone <- 2:1000
  essNoResmp100 <- essByParticles(12, theta, 100, testZone)
  essResmp100 <-
    essByParticles(12, theta, 100, testZone, algoResample = residualResampling)
  essNoResmp200 <- essByParticles(12, theta, 200, testZone)
  essResmp200 <-
    essByParticles(12, theta, 200, testZone, algoResample = residualResampling)
  
  par(mfrow = c(2, 2))
  plot(
    testZone,
    essNoResmp100 / testZone,
    main = "T=100 - Pas de rééchantillonnage",
    ylab = "ESS/N",
    xlab = "N",
    type = "l"
  )
  plot(
    testZone,
    essNoResmp200 / testZone,
    main = "T=200 - Pas de rééchantillonnage",
    ylab = "ESS/N",
    xlab = "N",
    type = "l"
  )
  plot(
    testZone,
    essResmp100 / testZone,
    main = "T=100 - Rééchantillonnage",
    ylab = "ESS/N",
    xlab = "N",
    type = "l"
  )
  plot(
    testZone,
    essResmp200 / testZone,
    main = "T=200 - Rééchantillonnage",
    ylab = "ESS/N",
    xlab = "N",
    type = "l"
  )
  par(mfrow = c(1, 1))
  
}

### Comportement du filtre

drawXiMean <- FALSE

if (drawXiMean) {
  par(mfrow = c(2, 2))
  for (i in c(10, 50, 100, 200)) {
    particleList <- lapply(1:10, function(i) {
      res <-
        likelihoodBootstrapParticleFilter(d,
                                          theta,
                                          simulY,
                                          100,
                                          algoResample = residualResampling,
                                          exportXi = TRUE)
    })
    allData <- c(simulXi, do.call(c, lapply(1:10, function(j) {
      return(particleList[[j]]$xi$mean)
    })))
    yLim <- c(min(allData), max(allData))
    print(yLim)
    plot(
      particleList[[1]]$xi$mean,
      ylim = yLim,
      main = i,
      xlab = "Temps t",
      ylab = "Xi",
      type = "l"
    )
    for (i in 2:10) {
      lines(particleList[[i]]$xi$mean)
    }
    points(simulXi)
  }
  par(mfrow = c(1, 1))
  
}

##################################################
##### ECHANTILLONNAGE DE THETA PAR PMCMC

burnin <- 100
pmcmcSize <- 200
N <- 20

simulThetaIid <-
  genThetaPosterior(theta,
                    d,
                    simulY,
                    N,
                    pmcmcSize,
                    componentWise = TRUE,
                    exportProba = TRUE)

plotSimulResult(simulThetaIid, exportProba = TRUE)

### Variance optimale

pmcmcSize <- 100
N <- 10

varianceList <- 10^(-6:6)
monitorChain <- list()
for (v in varianceList) {
  temp <-
    genThetaPosterior(
      theta,
      d,
      simulY,
      N,
      pmcmcSize,
      componentWise = TRUE,
      exportProba = TRUE,
      covariance = diag(c(v, rep(1, 3)))
    )
  temp <- formatResThetaPosterior(temp, exportProba = TRUE)
  monitorChain$corr <- append(monitorChain$corr,
                              acf(temp$beta1, lag.max = 1, plot = FALSE)[["acf"]][2, 1, 1])
  monitorChain$proba <-
    append(monitorChain$proba, mean(temp$proba[1,]))
}

best <-
  genThetaPosterior(
    theta,
    d,
    simulY,
    20,
    500,
    componentWise = TRUE,
    exportProba = TRUE,
    covariance = diag(c(.0001,.01,.001,.005))
  )
plotSimulResult(best, 0, 500)

####################################################################################
##### TOUTE LA SUITE DOIT ETRE MISE A JOUR

cov <- acf(simulThetaIid$beta1, lag.max = 1, plot = FALSE)

simulThetaIidRes <-
  genThetaPosterior(theta, d, simulY, N, pmcmcSize, algoResample = residualResampling)
simulThetaIidRes <- formatResThetaPosterior(simulThetaIidRes)
par(mfrow = c(2, 2))
for (param in getThetaNameList()) {
  plot(
    simulThetaIidRes[[param]][burnin:pmcmcSize],
    type = "l",
    main = param,
    xlab = "Itération",
    ylab = param
  )
}
par(mfrow = c(1, 1))

simulThetaIidIndep <-
  genThetaPosterior(theta, d, simulY, N, pmcmcSize, genNewProposal = genNewProposalSimpleIidIndep)
simulThetaIidIndep <- formatResThetaPosterior(simulThetaIidIndep)
par(mfrow = c(2, 2))
for (param in getThetaNameList()) {
  plot(
    simulThetaIidIndep[[param]][burnin:pmcmcSize],
    type = "l",
    main = param,
    xlab = "Itération",
    ylab = param
  )
}
par(mfrow = c(1, 1))

simulThetaIidResIndep <-
  genThetaPosterior(theta,
                    d,
                    simulY,
                    N,
                    pmcmcSize,
                    algoResample = residualResampling,
                    genNewProposal = genNewProposalSimpleIidIndep)
simulThetaIidResIndep <-
  formatResThetaPosterior(simulThetaIidResIndep)
par(mfrow = c(2, 2))
for (param in getThetaNameList()) {
  plot(
    simulThetaIidResIndep[[param]][burnin:pmcmcSize],
    type = "l",
    main = param,
    xlab = "Itération",
    ylab = param
  )
}
par(mfrow = c(1, 1))

for (param in getThetaNameList()) {
  print(param)
  print(mean(simulThetaIidResIndep[[param]]))
  print(mean(simulThetaIidResIndep[[param]][burnin:pmcmcSize]))
}

### Variance des MCMC
library(mcmc)
initseq(1:10)$var.con
