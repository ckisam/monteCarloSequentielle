







rm(list = ls())
# wd <- "C:/Users/Samuel/Documents/ENSAE - HMM/monteCarloSequentielle"
wd <- "/Users/Samuel/Desktop/ENSAE - 3A/HMM/Projet"
setwd(wd)
getwd()
source("CodeSource.R")

library(metRology)

step1 <- function(xi, theta, y) {
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

# logLikelihoodXi <- function(xi, param) {
logLikelihoodXi <-
  function(xi, mu, sigma, d, t, totalTime, beta1, beta2) {
    # mu <- param$mu
    # sigma <- param$sigma
    # d <- param$d
    # t <- param$t
    # T <- param$T
    # beta1 <- param$beta1
    # beta2 <- param$beta2
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

resStep1 <- step1(xi, theta, y)
