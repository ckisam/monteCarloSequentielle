##################################################
##### PROJET HMM ET SMC
##### Inference bayesienne par PMCMC
##### --- TESTS
##### (28/12/2017)
##################################################

rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
setwd(wd)
getwd()
source("inferenceBayesienne.R")

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
simulTheta <-
  genThetaPosterior(theta0,
                    d,
                    simulPhotoCount,
                    N,
                    nb,
                    genNewProposalSimpleIid)
simulTheta.format <- formatResThetaPosterior(simulTheta)
plot(
  simulTheta.format$beta2,
  type = "l",
  col = "navy",
  main = "Valeur de beta2 au cours des itérations",
  xlab = "Itération",
  ylab = "beta2"
)
