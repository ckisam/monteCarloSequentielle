################################"
### scheame A : standard gibbs sampler
###

wd <- "C:/Users/Dell/Desktop/cours ENSAE/chaines de Markov/projet"
setwd(wd)
getwd()

source("modele.R")

#############################

f <- function( xi, mu, sigma, beta, t){
  dt <- 12
  ret <- dnorm( x = xi, mean = mu, sd = sigma )*exp(-dt*exp(beta[1] + (beta[2] * i / t)+xi))
  return(ret)
}

schB_genXi <- function(xi, beta, delta, rho, mu, sigma, y, t){
  U <- 2
  X <- 0
  g <- 1
  f1 <- 1
  
  X <- rnorm( n = 1, mean = mu, sd = sigma)
  f1 <- f(X, mu, sigma, beta, t)
  U <- runif(1, min=0, max=1)
  g <- dnorm( x = X, mean = mu, sd = sigma )
  
  if( U > f1/g ){
    return(xi)
  }else{
    return(X)
  }
  
  
  
}


schB_genRhoDelta <- function(beta, t_xi, Y, dt){
  
  rho_new <- 10
  
  
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
  
  if( rho_new < 0.99 && rho_new > -0.99){
    return(c(rho_new,delta_new))
  }else{
    return(c(rho,delta))
  }
  
  return(c(rho_new,delta_new))
}



schB_genXiTimeSerie <- function(beta, delta, rho, t_xi, Y){
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
    
    t_xi[i] <- schB_genXi(t_xi[i], beta, delta, rho, mu, sigma, Y[i], i)
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

schB_genBeta <- function(beta, delta, rho, t_eta, Y){
  ZT <- matrix(0,ncol=0,nrow=2)
  for( i in 1:length(Y)){
    if( i == 1 ){
      Z <- c(sqrt(1-rho^2),sqrt(1-rho^2)*(i/length(Y)))
    }else{
      Z <- c(1-rho,(i/length(Y))+rho*(i-1)/length(Y))
    }
    
    ZT <- cbind(ZT,Z)
  }
  
  
  
  eta_tild <- vector(length=length(Y))
  for( i in 1:length(Y)){
    if( i == 1 ){
      eta_tild[i] <- sqrt(1-rho^2)*t_eta[i]
    }else{
      eta_tild[i] <- t_eta[i]-rho*t_eta[i-1]
    }
    
  }
  
  sol <- solve(ZT%*%t(ZT))
  
  beta_tild <- (sol%*%ZT)%*%eta_tild
    
  beta_new <- MASS:::mvrnorm(n = 1, mu = beta_tild, Sigma=sol*delta^2)
    
  return(beta_new)
}  

schB_genEta <- function(beta, delta, rho, t_xi, Y){
  t_eta <- vector(length=length(Y))
  for(i in 1:length(Y)){
    t_eta[i] <- t_xi[i] + beta[1] + beta[2]*(i/length(Y))
  }
  return(t_eta)
}

schB_genXiFromEtaBeta <- function(t_eta, beta, Y){
  t_xi <- vector(length=length(Y))
  for( i in 1:length(Y)){
    t_xi[i] <- t_eta[i] - beta[1] - beta[2]*(i/length(Y))
  }
  return(t_xi)  
}



############################
t <- 1000
beta <- c(0.1,0.1)
rho <- 0.1
delta <- 0.1
dt <- 10
d <- rep(10, t)


Y <- genPhotonCount(d, beta, delta, rho, t)
print(Y)
plot(Y)


time1 <- Sys.time()

beta <- c(0,1)
rho <- 0
delta <- 0.01


#### intialisation de la série xi #######
t_xi <- vector(length=length(Y))
t_xi[1] <- rnorm(n=1, mean=0, sd=sqrt(delta^2/(1-rho^2)))
for(i in 2:length(Y)){
  t_xi[i] <- rnorm(n=1, mean=rho*t_xi[i-1])
}
#print(t_xi)

t_beta1 <- vector(length=length(Y))
t_beta2 <- vector(length=length(Y))

t_rho <- vector(length=length(Y))
t_delta <- vector(length=length(Y))


#########################################

for( i in 1:1000){
  
  time2 <- Sys.time()
  
  t_xi <- schB_genXiTimeSerie(beta, delta, rho, t_xi, Y)
  
  t_eta <- schB_genEta(beta, delta, rho, t_xi, Y)
  
  beta <- schB_genBeta(beta, delta, rho, t_eta, Y)
  t_beta1[i] <- beta[1]
  t_beta2[i] <- beta[2]
  
  t_xi <- schB_genXiFromEtaBeta(t_eta, beta, Y)
  
  rho_delta <- schB_genRhoDelta(beta, t_xi, Y, dt)
  t_rho[i] <- rho_delta[1]
  t_delta[i] <- rho_delta[2]
  
  time3 <- Sys.time()
  
}

time4 <- Sys.time()


print((time4 - time1))

print(var(t_beta1[500:1000]))
print(var(t_beta2[500:1000]))
print(var(t_rho[500:1000]))
print(var(t_delta[500:1000]))

#plot(t_xi, type="l")
plot(t_beta1, type="l")
#plot(t_beta2, type="l")
#plot(t_rho, type="l")
#plot(t_delta, type="l")
