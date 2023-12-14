SPkernelts<-function(x,z,y,bw,no.lag,plt,k){
  library(condSURV)
  n    <-length(y)
  #STEP 1. PLM ESTIMATION IGNORING AUTOCORRELATION--------------------------------
  z_seq <- seq(min(z)-0.1,max(z)+0.1,length.out=n)      #Sequence for computation of Kernel
  W      <- matrix(0,n,n)
  for (j in 1:n){
    W[j,]<- NWW(z,z_seq[j],kernel="gaussian",bw)         #Nadaraya-Watson Wieghts 
  }
  xtil  <- (diag(n)-W)%*%x                                         #Partial Residual of X 
  ytil  <- (diag(n)-W)%*%y                                         #Partial Residual of Y
  betap <- solve(t(xtil)%*%xtil+k*diag(ncol(xtil)),tol=1e-100)%*%t(xtil)%*%ytil       #beta_pilot estimation
  fp    <- W%*%(y-x%*%betap)                                       #f pilot estimation
  yhatp <- x%*%betap+fp                                            #pilot yhat estimation
  res   <- y-yhatp                                                 #REsiduals for testing autocorrelation
  rho1 <- acf(res,pl=FALSE,lag=no.lag)                             #Autocorrelation parameter estimation
  rho  <- rho1$acf[2:(no.lag+1)]
  lr <- length(rho)+1
  rhot <- c(1,rho,rep(0,(n-lr))) 
  R<-toeplitz(rhot)        
  #E    <-(varu/(1-rho^2))*R                                        #Var-Cov Matrix
  invR <- solve(R)                                                 #Inverse of R matrix
  #------------------------------------------------------------------------------
  
  z_seq  <- seq(min(z)-0.01,max(z)+0.01,length.out=n)                             #Sequence for computation of Kernel
  W      <- matrix(0,n,n)
  for (j in 1:n){
    W[j,]<- NWW(z,z_seq[j],kernel="gaussian",bw)                                  #Nadaraya-Watson Wieghts 
  }
  xtil   <- (diag(n)-W)%*%x                                                       #Partial Residual of X 
  ytil   <- (diag(n)-W)%*%y                                                       #Partial Residual of Y
  betahat<- solve(t(xtil)%*%invR%*%xtil+k*diag(ncol(xtil)),tol=1e-100)%*%t(xtil)%*%invR%*%ytil       #beta_HAT estimation
  fhat   <- W%*%(y-x%*%betahat)                                                   #f estimation
  yhat   <- x%*%betahat+fhat                                                      #Fitted Values
  #HAT MATRIX-------------------------------------------------------------------
  H<-W+(diag(n)-W)%*%xtil%*%solve(t(xtil)%*%invR%*%xtil+k*diag(ncol(xtil)),tol=1e-100)%*%t(xtil)%*%invR%*%(diag(n)-W)
  
  if (plt==TRUE){
    plot(fhat,ylim=c(min(data$real_f-1),max(data$real_f+1)),type="l",lty=2,col=2,ylab="f(z)",xlab="x",main="Estimation of f(z)")
    par(new=TRUE)
    plot(data$real_f,type="l",col=1,ylab="f(x)",xlab="z",ylim=c(min(data$real_f-1),max(data$real_f+1)))
    legend("topleft",legend=c("KS estimate of f(z)","Real Function"),col=c(2,1),lty=c(2,1))
    grid()
  }
  KS           <-new.env()
  KS$yhat      <-yhat            #Estimated fitted values
  KS$W         <-W               #Nadaraya-Watson weights (based on Gaussian Kernel)
  KS$betahat   <-betahat         #Estimated Regression coefficients
  KS$fhat      <-fhat            #Estimated nonparametric function
  KS$H_mat     <-H               #Hat matrix of the estimated semiparametric time-serires model
  KS$R         <-R               #Autocorrelation matrix
  KS$invR      <-invR            #Inverse of R matrix
  return(KS)
}