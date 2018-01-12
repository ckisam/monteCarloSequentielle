################################"
### scheame A : standard gibbs sampler
###

wd <- "C:/Users/Dell/Desktop/cours ENSAE/chaines de Markov/projet"
setwd(wd)
getwd()

source("modele.R")

#############################

f <- function( xi, mu, sigma, beta, t){
  dt <- 0.1
  ret <- dnorm( x = xi, mean = mu, sd = sigma )*exp(-dt*exp(beta[1] + (beta[2] * i / t)+xi))
  return(ret)
}

schA_genXi <- function(beta, delta, rho, mu, sigma, y, t){
  U <- 2
  X <- 0
  g <- 1
  f1 <- 1
  while( U > f1/g ){
    
    X <- rnorm( n = 1, mean = mu, sd = sigma)
    f1 <- f(X, mu, sigma, beta, t)
    U <- runif(1, min=0, max=1)
    g <- dnorm( x = X, mean = mu, sd = sigma )

  }
  return(X)
  
}


schA_genRhoDelta <- function(beta, t_xi, Y, dt){
  
  rho_new <- 10
  
  while( rho_new > 0.99 || rho_new < -0.99){
    s1 <- 0
    s2 <- 0
    for(i in 2:length(Y)){
      s1 <- s1+(t_xi[i-1]*t_xi[i])
    }
    for(i in 2:(length(Y)-1)){
      s2 <- s2+(t_xi[i]^2)
    }
    rho <- s1/s2
    ##
    delta <- (1-rho^2)*t_xi[1]^2
    for(i in 2:length(Y)){
      delta <- delta + (t_xi[i] - rho*t_xi[i-1])^2
    }
    delta <- sqrt(delta)
    ##
    
    khi2 <- rchisq(n = 1, df=length(Y)-2)
    
    delta_new <- sqrt(delta^2 / khi2)
    
    sum <- 0
    for(i in 2:length(Y)){
      sum <- sum + t_xi[i]^2
    }

    rho_new <- rnorm(n=1, mean = rho, delta^2/sum)
    print(rho_new)
    
  }

  return(c(rho_new,delta_new))
}



schA_genXiTimeSerie <- function(beta, delta, rho, t_xi, Y){
  t_xi_old <- t_xi
  for(i in 1:length(Y)){
    if(i == 1){
      mu <-  Y[i]*delta^2 + t_xi[2]*rho
    }else{
      if(i == 1000){
        mu <- Y[i]*delta^2 + t_xi[i-1]*rho 
      }else{
        mu <- (Y[i]*delta^2 + (t_xi[i-1] + t_xi[i+1])*rho)/(1+rho^2)
      }
    }
    if(i==1 | i==1000){ 
      sigma <- delta
    }else{
      sigma <- sqrt(delta^2 / (1+rho^2))
    }
    
    t_xi[i] <- schA_genXi(beta, delta, rho, mu, sigma, Y[i], i)
  }
  return(t_xi)
}


f_beta <- function( beta, t_xi, rho, delta, Y, dt ){
  s <- 0
  for(i in 1:length(Y)){
    s <- s+ Y[i]*(beta[1] + beta[2]*(i/length(Y))) - dt*exp(beta[1] + beta[2]*(i/length(Y)) + t_xi[i])
  }
  
  return(exp(s))
}

g_beta <-  function( beta, t_xi, rho, delta, Y, dt ){
  s <- 0
  for(i in 1:length(Y)){
    s <- s+ Y[i]*(beta[1] + beta[2]*(i/length(Y))) 
  }
  
  return(exp(s))
}

schA_genBeta <- function(beta, delta, rho, t_xi, Y){
  
}  

############################
t <- 1000
beta <- c(0.05,0.05)
rho <- 0.1
delta <- 0.1
dt <- 0.1

Y <- genPhotonCount(d, beta, delta, rho, n)
print(Y)
plot(Y)


#### intialisation de la série xi #######
t_xi <- vector(length=length(Y))
t_xi[1] <- rnorm(n=1, mean=0, sd=sqrt(delta^2/(1-rho^2)))
for(i in 2:length(Y)){
  t_xi[i] <- rnorm(n=1, mean=rho*t_xi[i-1])
}
#print(t_xi)

#########################################

for( i in 1:1){
  
  t_xi = schA_genXiTimeSerie(beta, delta, rho, t_xi, Y)
  
  #beta = genBeta(t_xi, Y)
  
  rho_delta = schA_genRhoDelta(beta, t_xi, Y, dt)
  
  
  
}


plot(t_xi)
