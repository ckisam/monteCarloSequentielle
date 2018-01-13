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
  startUp <- list(
    beta1 = 5,
    beta2 = 5,
    delta = 5,
    rho = .2
  )
  startDown <- list(
    beta1 = -5,
    beta2 = -5,
    delta = -5,
    rho = .2
  )
  
  thetaRef <- list(
    beta1 = 0,
    beta2 = 1,
    delta = .1,
    rho = .5
  )
  
  par(mfrow = c(2, 2))
  for (nbPart in c(10, 50, 100, 200)) {
    particleList <- lapply(1:10, function(i) {
      res <-
        likelihoodBootstrapParticleFilter(d,
                                          startUp,
                                          simulY,
                                          nbPart,
                                          #algoResample = residualResampling,
                                          exportXi = TRUE)
    })
    allData <- c(simulXi, do.call(c, lapply(1:10, function(i) {
      return(particleList[[i]]$xi$mean)
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

startUp <- list(
  beta1 = 5,
  beta2 = 5,
  delta = 5,
  rho = .2
)
startDown <- list(
  beta1 = -5,
  beta2 = -5,
  delta = -5,
  rho = .2
)

thetaRef <- list(
  beta1 = 0,
  beta2 = 1,
  delta = .1,
  rho = .5
)

burnin <- 1000
pmcmcSize <- 3000
N <- 100
covariance <- diag(c(.1, .1, .1, .1))
posteriorStartUpThetaRef <-
  genThetaPosterior(startUp,
                    d,
                    genPhotonCount(d, thetaRef, n)$y,
                    N,
                    pmcmcSize,
                    covariance = covariance)
posteriorStartDownThetaRef <-
  genThetaPosterior(startDown,
                    d,
                    genPhotonCount(d, thetaRef, n)$y,
                    N,
                    pmcmcSize,
                    covariance = covariance)
posteriorStartUpThetaRefCompWise <-
  genThetaPosterior(
    startUp,
    d,
    genPhotonCount(d, thetaRef, n)$y,
    N,
    pmcmcSize,
    componentWise = TRUE,
    covariance = covariance
  )
posteriorStartDownThetaRefCompWise <-
  genThetaPosterior(
    startDown,
    d,
    genPhotonCount(d, thetaRef, n)$y,
    N,
    pmcmcSize,
    componentWise = TRUE,
    covariance = covariance
  )
posteriorStartUpThetaRefResample <-
  genThetaPosterior(
    startUp,
    d,
    genPhotonCount(d, thetaRef, n)$y,
    N,
    pmcmcSize,
    algoResample = residualResampling,
    covariance = covariance
  )
posteriorStartDownThetaRefResample <-
  genThetaPosterior(
    startDown,
    d,
    genPhotonCount(d, thetaRef, n)$y,
    N,
    pmcmcSize,
    algoResample = residualResampling,
    covariance = covariance
  )
posteriorStartUpThetaRefCompWiseResample <-
  genThetaPosterior(
    startUp,
    d,
    genPhotonCount(d, thetaRef, n)$y,
    N,
    pmcmcSize,
    algoResample = residualResampling,
    componentWise = TRUE,
    covariance = covariance
  )
posteriorStartDownThetaRefCompWiseResample <-
  genThetaPosterior(
    startDown,
    d,
    genPhotonCount(d, thetaRef, n)$y,
    N,
    pmcmcSize,
    algoResample = residualResampling,
    componentWise = TRUE,
    covariance = covariance
  )
plotSimulResult(posteriorStartUpThetaRef, 1, pmcmcSize) # OK - 2m10
plotSimulResult(posteriorStartDownThetaRef, 1, pmcmcSize) # OK - 2m10
plotSimulResult(posteriorStartUpThetaRefCompWise, 1, pmcmcSize) # OK - 10 min
plotSimulResult(posteriorStartDownThetaRefCompWise, 1, pmcmcSize) # OK - 9 min
plotSimulResult(posteriorStartUpThetaRefResample, 1, pmcmcSize) # Pb
plotSimulResult(posteriorStartDownThetaRefResample, 1, pmcmcSize) # OK - 1m40
plotSimulResult(posteriorStartUpThetaRefCompWiseResample, 1, pmcmcSize) # Pb
plotSimulResult(posteriorStartDownThetaRefCompWiseResample, 1, pmcmcSize) # OK - 12 min

printMeanAndVariance(posteriorStartDownThetaRef, 2000, 3000)
# [1] "Beta1 = -0.0829228027374298 - Beta2 = 1.08849766505897
# - Delta = 0.0586912867078934 - rho = 0.502502020583862"
# [2] "Beta1 = 0.00910785328294363 - Beta2 = 0.0189724875220476
# - Delta = 0.012629811066568 - rho = 0.0341755192059229"
printMeanAndVariance(posteriorStartUpThetaRef, 2000, 3000)
# [1] "Beta1 = -0.0133024577986279 - Beta2 = 0.979937043440443
# - Delta = 0.0707748225827671 - rho = 0.273032101731567"
# [2] "Beta1 = 0.00617256602560364 - Beta2 = 0.0140999154410152
# - Delta = 0.0129860162379289 - rho = 0.0460113249783819"
printMeanAndVariance(posteriorStartDownThetaRefCompWise, 2000, 3000)
# [1] "Beta1 = 0.0161758936613691 - Beta2 = 0.950449829870693
# - Delta = -0.00662497180056832 - rho = 0.470665739825585"
# [2] "Beta1 = 0.00939046606511996 - Beta2 = 0.0132710094463102
# - Delta = 0.0188409659562685 - rho = 0.0359275625267412"
printMeanAndVariance(posteriorStartUpThetaRefCompWise, 2000, 3000)
# [1] "Beta1 = -0.0703026604734542 - Beta2 = 1.03301050083062
# - Delta = -0.00619558098948477 - rho = 0.424406085628546"
# [2] "Beta1 = 0.0122634131617299 - Beta2 = 0.0198791181570819
# - Delta = 0.0154969563142585 - rho = 0.0356930053594359"
printMeanAndVariance(posteriorStartUpThetaRefResample, 2000, 3000)
printMeanAndVariance(posteriorStartDownThetaRefResample, 2000, 3000)
# [1] "Beta1 = -5.7932618573769 - Beta2 = -3.65355548226205
# - Delta = -0.432720031410151 - rho = 0.986923146205778"
# [2] "Beta1 = 0 - Beta2 = 0 - Delta = 0 - rho = 0"
printMeanAndVariance(posteriorStartUpThetaRefCompWiseResample, 2000, 3000)
printMeanAndVariance(posteriorStartDownThetaRefCompWiseResample, 2000, 3000)
# [1] "Beta1 = -0.15348644005593 - Beta2 = 1.3121308928661
# - Delta = -0.0151397782887926 - rho = 0.661796838760115"
# [2] "Beta1 = 0.0151891523734407 - Beta2 = 0.0238601777355272
# - Delta = 0.00571777566039334 - rho = 0.0276327499033897"

################################################################################
########## ECHANTILLONNAGE MCMC AVEC GIBBS STANDARD

theta <- list(
  beta1 = 0,
  beta2 = 1,
  delta = .1,
  rho = .5
)
n <- 100
d <- 12
d <- rep(d, n)

photonCount <- genPhotonCount(d, theta, n)
y <- photonCount$y
xi <- photonCount$xi

startUp <- list(
  beta1 = 5,
  beta2 = 5,
  delta = 5,
  rho = .2
)
startDown <- list(
  beta1 = -5,
  beta2 = -5,
  delta = -5,
  rho = .2
)

burnin <- 2000
pmcmcSize <- 3000
# Up
gibbsSimulUp <-
  standardGibbsSamplerA(d, y, startUp, pmcmcSize) # OK 1m50
plotSimulResult(gibbsSimulUp, 1, 300)
printMeanAndVariance(gibbsSimulUp, burnin, pmcmcSize)
# [1] "Beta1 = -0.0289099822690155 - Beta2 = 1.03887640693943
# - Delta = 0.0932744991082598 - rho = 0.412180098310289"
# [2] "Beta1 = 0.00144189497925433 - Beta2 = 0.00224081024030777
# - Delta = 0.000207063878879559 - rho = 0.000299861605264446"
# Down
gibbsSimulDown <-
  standardGibbsSamplerA(d, y, startDown, pmcmcSize) # OK 1m12
plotSimulResult(gibbsSimulDown, 1, 300)
printMeanAndVariance(gibbsSimulDown, burnin, pmcmcSize)
# [1] "Beta1 = -0.0289761239876972 - Beta2 = 1.03863291752434
# - Delta = 0.0931687209121405 - rho = 0.412300970860029"
# [2] "Beta1 = 0.0015144481833517 - Beta2 = 0.00214089013021137
# - Delta = 0.000184221933683781 - rho = 0.000241460521292388"
