################################"
### scheame A : standard gibbs sampler
###

wd <- "C:/Users/Dell/Desktop/cours ENSAE/chaines de Markov/projet"
setwd(wd)
getwd()

source("modele.R")

#############################

p_beta_post <- function(beta , xi, Y, dt, t) {
  print(c(beta[1] , beta[2], xi[1], xi[2], Y[1], dt, t))
  l <- 1
  for (i in 1:t){
    lambda <- dt*exp(beta[1] + (beta[2] * i / t) + xi[i] )
    l <- l*(dpois(x = Y[i], lambda = lambda))
  }
  return(l)
}

p_rho_delta_post <- function(rho, delta, xi, t ){
  #print(c(rho, delta, t, xi[1], xi[2]))
  l <- dnorm(x = xi[1] ,mean = 0, sd = sqrt(delta^2/(1-rho^2)))
  if( t > 1 ){
    for(i in 2:t){
      l <- l*dnorm( x = xi[i], mean = (rho*xi[i-1]), sd = abs(delta))
    }
  }
  
  return(l)
}

# conditional density for metropolis-hasting
q <- function(x,y){
  v <- dnorm( x = x, mean = y, sd = 1)
  return(v)
}

f1 <- function(i, t){
  return(i/t) 
}


genBeta <- function(t_nu, t_xi, t){
  dt = data.frame(y=(t_xi[1:(t-1)]-t_nu[1:(t-1)]),x1=rep(1,each=(t-1)), x2=sapply(1:(t-1), FUN=f1, t=(t-1)))
  print(dt)
  lmMod <- lm(y ~ x1 + x2, data=dt)
  coefs <- coefficients(lmMod)
  print(c(coefs[[1]],coefs[[2]]))
  return(c(coefs[[1]],coefs[[2]]))
}

genRhoDelta <- function( rhoDelta, xi, t ){
  #print('genRhoDelta')
  #Yt <- rnorm(n = 2, mean = rhoDelta, sd = 1)
  Yt <- c(runif(n = 1, min = -0.99, max = 0.99),runif(n = 1, min = -0.99, max = 0.99))
  #Yt <- MASS:::mvrnorm(n = 1, mu = rhoDelta, Sigma=cbind(c(1,0),c(0,1)))
  #print(Yt)
  #print(rhoDelta)
  #print(p_rho_delta_post(rhoDelta[1], rhoDelta[2], xi, t ))
  #print(p_rho_delta_post(Yt[1], Yt[2], xi, t ))
  p1 <- p_rho_delta_post(Yt[1], Yt[2], xi, t )
  p2 <- p_rho_delta_post(rhoDelta[1], rhoDelta[2], xi, t )
  #print(p1)
  #print(p2)
  if( is.nan(p1/p2) ){
    r <- 1
  }else{
    r <- min( (p1/p2), 1 )
  }
  v <- runif(1, min = 0, max = 1)
  #print(r)
  #print(v)
  
  y <- rhoDelta
  if( v < r ){
    y <- Yt
  } 
  return(y)
} 

genRho <- function( rhoDelta, xi, t ){
  #print('genRhoDelta')
  #Yt <- rnorm(n = 2, mean = rhoDelta, sd = 1)
  Yt <- runif(n = 1, min = -0.99, max = 0.99)
  #Yt <- MASS:::mvrnorm(n = 1, mu = rhoDelta, Sigma=cbind(c(1,0),c(0,1)))
  #print(Yt)
  #print(rhoDelta)
  #print(p_rho_delta_post(rhoDelta[1], rhoDelta[2], xi, t ))
  #print(p_rho_delta_post(Yt[1], Yt[2], xi, t ))
  p1 <- p_rho_delta_post(Yt, rhoDelta[2], xi, t )
  p2 <- p_rho_delta_post(rhoDelta[1], rhoDelta[2], xi, t )
  #print(p1)
  #print(p2)
  if( is.nan(p1/p2) ){
    r <- 1
  }else{
    r <- min( (p1/p2), 1 )
  }
  v <- runif(1, min = 0, max = 1)
  #print(r)
  #print(v)
  
  y <- rhoDelta[1]
  if( v < r ){
    y <- Yt
  } 
  return(y)
} 


genDelta <- function( rhoDelta, xi, t ){
  #print('genRhoDelta')
  #Yt <- rnorm(n = 2, mean = rhoDelta, sd = 1)
  Yt <- runif(n = 1, min = -0.99, max = 0.99)
  #Yt <- MASS:::mvrnorm(n = 1, mu = rhoDelta, Sigma=cbind(c(1,0),c(0,1)))
  #print(Yt)
  #print(rhoDelta)
  #print(p_rho_delta_post(rhoDelta[1], rhoDelta[2], xi, t ))
  #print(p_rho_delta_post(Yt[1], Yt[2], xi, t ))
  p1 <- p_rho_delta_post(rhoDelta[1], Yt, xi, t )
  p2 <- p_rho_delta_post(rhoDelta[1], rhoDelta[2], xi, t )
  #print(p1)
  #print(p2)
  if( is.nan(p1/p2) ){
    r <- 1
  }else{
    r <- min( (p1/p2), 1 )
  }
  v <- runif(1, min = 0, max = 1)
  #print(r)
  #print(v)
  
  y <- rhoDelta[2]
  if( v < r ){
    y <- Yt
  } 
  return(y)
} 


#############################"
# initialisation des parametres
#

t <- 1000
beta <- c(0.5,0.5)
rho <- 0.1
delta <- 1
dt <- 1

Y <- genPhotonCount(d, beta, delta, rho, n)
plot(Y)

t_xi <- vector(length=length(Y))
t_nu <- vector(length=length(Y))


xi <- rnorm(n = 1, mean = 0, sd = sqrt(delta^2 / (1 - rho^2)))
t_xi[1] <- xi


nu <- rnorm(n = 1, mean =beta[0]+beta[1]*(1/t), sd = delta^2 / (1 - rho^2))
t_nu[1] <- xi



rho <- genRho(c(rho,delta), t_xi, 1 )
delta <- genDelta(c(rho,delta), t_xi, 1 )
rho_delta <- c(rho,delta)


for (i in 2:100) {
  rho <- rho_delta[1]
  delta <- rho_delta[2]
  
  #print("xi params")
  #print(c(rho,xi,abs(delta)))
  xi <- rnorm( n = 1, mean = rho*xi, sd = abs(delta))
  t_xi[i] <- xi
  
  beta <- genBeta(t_nu, t_xi, i)
  
  nu <- xi + beta[0] + beta[1]*(1/t)
  t_nu[i] <- nu
  
  
  print("beta")
  print(beta)
  
  rho <- genRhoDelta( rho_delta, t_xi, i )
  rho_delta <- c(rho, rho_delta[2])
  delta <- genRhoDelta( rho_delta, t_xi, i )
  rho_delta <- c(rho, delta)
  
  #print("rho_delta")
  #print(rho_delta)
}

print(t_xi)
plot(t_xi)

