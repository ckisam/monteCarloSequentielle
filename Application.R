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
lkh <- likelihoodBootstrapParticleFilter(d, theta, simulY, N)

### Comportement des particules
resLkh100 <-
  likelihoodBootstrapParticleFilter(d, theta, simulY, 100, exportXi = 100)

plot(
  resLkh100$xi$mean,
  main = "Processus Xi estimé",
  xlab = "Temps t",
  ylab = "Xi",
  # ylim = c(
  #   min(resLkh100$xi$mean, simulXi),
  #   max(resLkh100$xi$mean, simulXi)
  # )
)
#lines(simulXi)
lines(resLkh100$xi$mean-(qt(.975,df=n)*resLkh100$xi$sd/sqrt(n)),col="red")
lines(resLkh100$xi$mean+(qt(.975,df=n)*resLkh100$xi$sd/sqrt(n)),col="red")
      