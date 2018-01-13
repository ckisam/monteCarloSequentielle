rm(list = ls())
wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
# wd <- "/Users/Samuel/Desktop/ENSAE - 3A/HMM/Projet"
setwd(wd)
getwd()
source("CodeSource.R")

library(metRology)
library(mvtnorm)

standardGibbsSamplerA <- function(d, y, theta, nb) {
  start <- getStart()
  pct <- nb %/% 10
  logStart("ALgorithme Gibbs Standard (Schema A)")
  logWithTime(start,
              paste("On cherche à engendrer un echantillon MCMC de taille ", nb, sep = ""))
  thetaCur <- theta
  thetaCur$delta <- abs(thetaCur$delta)
  res <- list()
  res[[1]] <- thetaCur
  xiCur <-
    c(rnorm(1, sd = (thetaCur$delta / sqrt(1 - theta$rho ^ 2))))
  for (i in 2:length(y)) {
    xiCur <-
      append(xiCur, rnorm(1, mean = theta$rho * xiCur[length(xiCur)], sd = thetaCur$delta))
  }
  for (i in 1:nb) {
    xiCur <- step1(d, xi, theta, y)
    thetaCur <- step2A(d, thetaCur, xi, y)
    thetaCur <- step3S(d, thetaCur, xi, y)
    res[[i + 1]] <- thetaCur
    if (i %% pct == 0) {
      logWithTime(start, paste((i / pct) * 10, "% généré", sep = ""))
    }
  }
  logTotalTime(start)
  return(res)
}

step3S <- function(d, theta, xi, y) {
  rhoMode <-
    sum(xi[2:length(xi)] * xi[1:(length(xi) - 1)]) / sum((xi[2:(length(xi) -
                                                                  1)]) ^ 2)
  deltaMode <-
    (1 - rhoMode ^ 2) * xi[1] ^ 2 + sum((xi[2:length(xi)] - rhoMode * xi[1:(length(xi) -
                                                                              1)]) ^ 2)
  khi2 <- rchisq(1, df = (length(xi) - 2))
  newDelta <- deltaMode / khi2
  newRho <-
    rnorm(1, mean = rhoMode, sd = (newDelta / sum(xi[2:(length(xi) - 1)] ^
                                                    2)))
  if (abs(newRho) <= .99) {
    theta$delta <- sqrt(newDelta)
    theta$rho <- newRho
  }
  return(theta)
}

step2A <- function(d, theta, xi, y) {
  require(mvtnorm)
  require(expm)
  beta <- c(theta$beta1, theta$beta2)
  mode <-
    optim(
      beta,
      fn = logLikelihoodBeta,
      method = "BFGS",
      hessian = TRUE,
      d = d,
      y = y,
      xi = xi
    )
  scale <- solve(sqrtm(mode$hessian))
  betaMode <- mode$par
  newBeta <-
    as.vector(betaMode + scale %*% as.vector(rmvt(
      1,
      sigma = diag(2),
      df = 5,
      delta = rep(0, 2)
    )))
  ratio <-
    logLikelihoodBeta(beta, d, y, xi) - logLikelihoodBeta(newBeta, d, y, xi)
  ratio <-
    ratio + dmvt(beta,
                 delta = betaMode,
                 sigma = scale,
                 df = 5) - dmvt(newBeta,
                                delta = betaMode,
                                sigma = scale,
                                df = 5)
  proba <- runif(1)
  if (log(proba) <= ratio) {
    beta <- newBeta
  }
  theta$beta1 <- beta[1]
  theta$beta2 <- beta[2]
  return(theta)
}

logLikelihoodBeta <- function(beta, d, y, xi) {
  logLkh <- do.call(c, lapply(1:length(xi), function(i) {
    xb <- beta[1] + beta[2] * (i / length(xi))
    res <- y[i] * (xb) - d[i] * exp(xb + xi[i])
    return(res)
  }))
  logLkh <- sum(logLkh)
  return(-logLkh)
}

step1 <- function(d, xi, theta, y) {
  muSigma <-
    computeMuAndSigma(xi, list(
      y = y,
      delta = theta$delta,
      rho = theta$rho
    ))
  mode <- lapply(1:length(xi), function(i) {
    m <-
      optim(
        par = xi[i],
        fn = logLikelihoodXi,
        method = "BFGS",
        hessian = TRUE,
        mu = muSigma$mu[i],
        sigma = muSigma$sigma[i],
        d = d[i],
        t = i,
        totalTime = length(xi),
        beta1 = theta$beta1,
        beta2 = theta$beta2
      )
    return(m)
  })
  hessian <- do.call(c, lapply(mode, function(m) {
    return(m$hessian)
  }))
  scale <- 1 / sqrt(hessian)
  xiMode <- do.call(c, lapply(mode, function(m) {
    return(m$par)
  }))
  xiNew <- xiMode + rt(n = length(xi), df = 5) * scale
  proba <- runif(length(xi))
  res <- do.call(c, lapply(1:length(xi), function(i) {
    ratio <-
      logLikelihoodXi(xi[i],
                      muSigma$mu[i],
                      muSigma$sigma[i],
                      d[i],
                      i,
                      length(xi),
                      theta$beta1,
                      theta$beta2)
    ratio <- ratio - logLikelihoodXi(xiNew[i],
                                     muSigma$mu[i],
                                     muSigma$sigma[i],
                                     d[i],
                                     i,
                                     length(xi),
                                     theta$beta1,
                                     theta$beta2)
    ratio <-
      ratio - dt.scaled(xiNew[i],
                        df = 5,
                        mean = xiMode[i],
                        sd = scale[i])
    ratio <-
      ratio + dt.scaled(xi[i],
                        df = 5,
                        mean = xiMode[i],
                        sd = scale[i])
    if (log(proba[i]) <= ratio) {
      return(xiNew[i])
    } else{
      return(xi[i])
    }
  }))
  return(res)
}

computeMuAndSigma <- function(xi, param) {
  y <- param$y
  delta <- param$delta
  rho <- param$rho
  
  res <- list()
  res$mu <- c()
  res$sigma <- c()
  T <- length(xi)
  for (t in 1:T) {
    if (t == 1) {
      mu <- y[1] * delta ^ 2 + xi[2] * rho
      sigma <- delta ^ 2
    } else if (t == T) {
      mu <- y[T] * delta ^ 2 + xi[T - 1] * rho
      sigma <- delta ^ 2
    } else {
      mu <-
        (y[t] * delta ^ 2 + (xi[t - 1] + xi[t + 1]) * rho) / (1 + rho ^ 2)
      sigma <- (delta ^ 2) / (1 + rho ^ 2)
    }
    res$mu <- append(res$mu, mu)
    res$sigma <- append(res$sigma, sigma)
  }
  return(res)
}

logLikelihoodXi <-
  function(xi, mu, sigma, d, t, totalTime, beta1, beta2) {
    logLkh <- -((xi - mu) ^ 2) / (2 * sigma ^ 2)
    logLkh <- logLkh - d * exp(beta1 + beta2 * (t / totalTime) + xi)
    return(-logLkh)
  }

################################################################################
########## TEST

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
gibbsSimulUp <- standardGibbsSamplerA(d, y, startUp, pmcmcSize) # OK 1m50
plotSimulResult(gibbsSimulUp, 1, 300)
printMeanAndVariance(gibbsSimulUp, burnin, pmcmcSize)
# [1] "Beta1 = -0.0289099822690155 - Beta2 = 1.03887640693943
# - Delta = 0.0932744991082598 - rho = 0.412180098310289"        
# [2] "Beta1 = 0.00144189497925433 - Beta2 = 0.00224081024030777
# - Delta = 0.000207063878879559 - rho = 0.000299861605264446"
# Down
gibbsSimulDown <- standardGibbsSamplerA(d, y, startDown, pmcmcSize) # OK 1m12
plotSimulResult(gibbsSimulDown, 1, 300)
printMeanAndVariance(gibbsSimulDown, burnin, pmcmcSize)
# [1] "Beta1 = -0.0289761239876972 - Beta2 = 1.03863291752434
# - Delta = 0.0931687209121405 - rho = 0.412300970860029"       
# [2] "Beta1 = 0.0015144481833517 - Beta2 = 0.00214089013021137
# - Delta = 0.000184221933683781 - rho = 0.000241460521292388"

