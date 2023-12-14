library(pracma)
library(pspline)
library(mgcv)
library(condSURV)
###############################################################################
#FUNCTIONS---------------------------------------------------------------------
simdata <- function(n,rho,cor,pairplot){
  library(Rlab)   
  library(MASS)
  library(tidyverse)
  library(GGally)
  # length of rho means the degree of AR process
  #numexp denotes the number of explanatory variables
  #n sample size
  
  s     <- 0   
  s2    <- 0                 #Explanatory variable for the nonparametric component
  y     <- 0                 #Response variable 
  
  for (i in 1:n){            #Generating uniform nonparametic variable
    s[i]<-1.2*(i-0.5)/n
  }
  for (i in 1:(n/5)){        #Generating uniform nonparametic variable
    s2[i]<-(i-0.5)/(n/5)
  }
  e<-arima.sim(model=list(ar=rho), n=n)           #generating error terms
  
  sigma<-rbind(c(1,cor), c(cor,1))
  mu<-c(0, 0) 
  df<-as.data.frame(mvrnorm(n=n, mu=mu, Sigma=sigma))
  if (pairplot==TRUE){
    ggpairs(df)
  }
  x<-as.matrix(df)      #Parametric component
  beta<-c(-1,1.5)                                   #parametric coefficients
  #g1<-6*(s*sin(6*s))+10                              #Nonparametric curve no.1
  #Nonparametric season effect generation------------------------------------
  const <- 6
  seasoneff<-const*s2*sin(s2)^2               #const is a constant for a wideness of x axis.
  g2<-t(repmat(seasoneff,1,5))
  g<-g2
  
  y     <-x%*%beta+g+e                        #generating completely observed data 
  
  data<-new.env()
  data$y      <-y                 #Output of the function 1: Completely observed
  data$z      <-s                 #Output of the function 6: Nonparametric covarites
  data$real_f <-g                 #Output of the function 8: Real smooth function to be estimated (nonparametric component)
  data$errors <-e                 #AR(1) error terms
  data$beta   <-beta              #Parametric coefs.
  data$exp.var<-x                 #parametric covariates n X 3 (p=3)
  data$rho    <-rho               #autor.corr param. for AR(p) errors.
  data$AR_deg <-length(rho)
  return(data)  
}