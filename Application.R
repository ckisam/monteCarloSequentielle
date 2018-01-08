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

N <- 100

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

N <- 1000

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

##################################################
##### ECHANTILLONNAGE DE THETA PAR PMCMC

burnin <- 1000
pmcmcSize <- 2000
N <- 100

simulThetaIid <- genThetaPosterior(theta, d, simulY, N, pmcmcSize)
simulThetaIid <- formatResThetaPosterior(simulThetaIid)
par(mfrow = c(2, 2))
for (param in getThetaNameList()) {
  plot(
    simulThetaIid[[param]][burnin:pmcmcSize],
    type = "l",
    main = param,
    xlab = "Itération",
    ylab = param
  )
}
par(mfrow = c(1, 1))

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
