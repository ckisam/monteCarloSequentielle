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
N <- 100
covariance <- diag(c(.1, .1, .1, .1))

simulThetaIid <-
  genThetaPosterior(
    theta,
    d,
    simulY,
    N,
    pmcmcSize,
    componentWise = TRUE,
    exportProba = TRUE,
    covariance = covariance
  )

plotSimulResult(simulThetaIid, burnin, pmcmcSize, exportProba = TRUE)

### Variance optimale

# for (v in varianceList) {
#   temp <-
#     genThetaPosterior(
#       theta,
#       d,
#       simulY,
#       N,
#       pmcmcSize,
#       componentWise = TRUE,
#       exportProba = TRUE,
#       covariance = diag(c(v, rep(1, 3)))
#     )
#   temp <- formatResThetaPosterior(temp, exportProba = TRUE)
#   monitorChain$corr <- append(monitorChain$corr,
#                               acf(temp$beta1, lag.max = 1, plot = FALSE)[["acf"]][2, 1, 1])
#   monitorChain$proba <-
#     append(monitorChain$proba, mean(temp$proba[1,]))
# }

pmcmcSize <- 100
N <- 10

# varianceList <- 10 ^ (-6:6)
varianceList <- c(.001, 1, 100)
monitorChain <- list()

for (i in 1:length(varianceList)) {
  monitorChain[[i]] <- list()
  for (j in 1:length(varianceList)) {
    monitorChain[[i]][[j]] <- list()
    for (k in 1:length(varianceList)) {
      monitorChain[[i]][[j]][[k]] <- list()
      for (l in 1:length(varianceList)) {
        monitorChain[[i]][[j]][[k]][[l]] <- list()
        covariance <- diag(c(
          varianceList[i],
          varianceList[j],
          varianceList[k],
          varianceList[l]
        ))
        temp <-
          genThetaPosterior(
            theta,
            d,
            simulY,
            N,
            pmcmcSize,
            componentWise = TRUE,
            exportProba = TRUE,
            covariance = covariance
          )
        temp <- formatResThetaPosterior(temp, exportProba = TRUE)
        for (param in getThetaNameList()) {
          monitorChain[[i]][[j]][[k]][[l]]$corr[[param]] <-
            acf(temp[[param]], lag.max = 1, plot = FALSE)[["acf"]][2, 1, 1]
        }
        monitorChain[[i]][[j]][[k]][[l]]$proba <-
          append(monitorChain$proba, mean(temp$proba[1, ]))
      }
    }
  }
}

maxProba <- list()
for (i in 1:length(varianceList)) {
  for (j in 1:length(varianceList)) {
    for (k in 1:length(varianceList)) {
      for (l in 1:length(varianceList)) {
        if (monitorChain[[i]][[j]][[k]][[l]]$proba >= .2) {
          maxProba[[length(maxProba) + 1]] <- c(i, j, k, l)
        }
      }
    }
  }
}

for(l in maxProba) {
  print(l)
  print(monitorChain[[l[1]]][[l[2]]][[l[3]]][[l[4]]]$proba)
}


# alterTheta <- list(
#   beta1 = 10,
#   beta2 = -5,
#   delta = 3,
#   rho = 0.9
# )
alterTheta <- list(
  beta1 = 0,
  beta2 = 0,
  delta = 10,
  rho = 0.7
)
burnin <- 1000
pmcmcSize <- 2000
N <- 100
# covariance <- diag(c(.001, .01, .01, .001))
covariance <- diag(c(.1, .1, .1, .1))

testAVirer <-
  genThetaPosteriorNew(
    list(
      beta1 = 5,
      beta2 = 5,
      delta = 5,
      rho = .2
    ),
    # list(
    #   beta1 = 0,
    #   beta2 = 1,
    #   delta = .1,
    #   rho = .5
    # ),
    d,
    genPhotonCount(d, list(
      beta1 = 0,
      beta2 = 1,
      delta = .1,
      rho = .5
    ), n)$y,
    N,
    pmcmcSize,
    exportProba = TRUE,
    #algoResample = residualResampling,
    log = TRUE,
    covariance = covariance
  )
plotSimulResult(testAVirer, 1, pmcmcSize)

testAVirer <- formatResThetaPosterior(testAVirer)
print(mean(testAVirer$beta1[burnin:pmcmcSize]))
print(mean(testAVirer$beta2[burnin:pmcmcSize]))
print(mean(testAVirer$delta[burnin:pmcmcSize]))
print(mean(abs(testAVirer$rho[burnin:pmcmcSize])))

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
