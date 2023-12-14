Selection_bw_k<-function(x,z,y,bw_seq,k_seq,no.lag,plt){
  library(condSURV)
  library(psych)
  n    <-length(y)
  #STEP 1. PLM ESTIMATION IGNORING AUTOCORRELATION--------------------------------
  z_seq <- seq(min(z)-0.1,max(z)+0.1,length.out=n)      #Sequence for computation of Kernel
  W      <- matrix(0,n,n)
  for (j in 1:n){
    W[j,]<- NWW(z,z_seq[j],kernel="gaussian",bw=mean(bw_seq))         #Nadaraya-Watson Wieghts 
  }
  xtil  <- (diag(n)-W)%*%x                                         #Partial Residual of X 
  ytil  <- (diag(n)-W)%*%y                                         #Partial Residual of Y
  betap <- solve(t(xtil)%*%xtil,tol=1e-100)%*%t(xtil)%*%ytil       #beta_pilot estimation
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
  no.try <- length(bw_seq)
  no.k   <- length(k_seq)
  
  GCV  <- matrix(0,no.try,no.k)
  AICc <- matrix(0,no.try,no.k)
  RECP <- matrix(0,no.try,no.k)
  BIC  <- matrix(0,no.try,no.k)
  
  for (jt in 1:no.try){
    for (j2 in 1:no.k){
      z_seq  <- seq(min(z)-0.01,max(z)+0.01,length.out=n)                             #Sequence for computation of Kernel
      W      <- matrix(0,n,n)
      for (j in 1:n){
        W[j,]<- NWW(z,z_seq[j],kernel="gaussian",bw_seq[jt])                                  #Nadaraya-Watson Wieghts 
      }
      xtil   <- (diag(n)-W)%*%x                                                       #Partial Residual of X 
      ytil   <- (diag(n)-W)%*%y                                                       #Partial Residual of Y
      betahat<- solve(t(xtil)%*%invR%*%xtil+k_seq[j2]*diag(ncol(xtil)),tol=1e-100)%*%t(xtil)%*%invR%*%ytil       #beta_HAT estimation
      fhat   <- W%*%(y-x%*%betahat)                                                   #f estimation
      yhat   <- x%*%betahat+fhat                                                      #Fitted Values
      #HAT MATRIX-------------------------------------------------------------------
      H<-W+(diag(n)-W)^2%*%xtil%*%solve(t(xtil)%*%invR%*%xtil+k_seq[j2]*diag(ncol(xtil)),tol=1e-100)%*%t(xtil)%*%invR
      
      #varp <- t((diag(n)-H)%*%ytil)%*%((diag(n)-H)%*%ytil)/(tr(H))
      varp <- var(H%*%ytil-ytil)
      
      GCV[jt,j2]  <- ((1/n)*(t((diag(n)-H)%*%y)%*%((diag(n)-H)%*%y)))/((1/n)*(tr(diag(n)-H))^2)
      AICc[jt,j2] <- log(abs((t((diag(n)-H)%*%y)%*%((diag(n)-H)%*%y))/(n-tr(H))))+1+((2*(tr(H)+1))/(n-tr(H)-2))
      BIC[jt,j2]  <- ((t((diag(n)-H)%*%y)%*%((diag(n)-H)%*%y))/n)+(log(n)/n)*tr(H)
      RECP[jt,j2] <- (1/n-tr(H))*(t((H-diag(n))%*%y)%*%((H-diag(n))%*%y)+varp*tr(H%*%t(H)))
      
    }
  }
  
  for (j3 in 1:no.try){
    for (j4 in 1:no.k){
      if (GCV[j3,j4]==min(GCV)){
        bw_GCV <- bw_seq[j3]
        k_GCV  <- k_seq[j4] 
      }
      if (AICc[j3,j4]==min(AICc)){
        bw_AICc <- bw_seq[j3]
        k_AICc  <- k_seq[j4] 
      }
      if (BIC[j3,j4]==min(BIC)){
        bw_BIC <- bw_seq[j3]
        k_BIC  <- k_seq[j4] 
      }
      if (RECP[j3,j4]==min(RECP)){
        bw_RECP <- bw_seq[j3]
        k_RECP  <- k_seq[j4] 
      }
    }
  }
  
  pairs <- data.frame(bw_GCV,k_GCV,bw_AICc,k_AICc,bw_BIC,k_BIC,bw_RECP,k_RECP)
  if (plt==TRUE){
    par(mfrow=c(2,2))
    par(mar=c(3,2,1,1))
    persp(bw_seq, k_seq, GCV,
          main="GCV",
          zlab = "GCV score",
          theta = 25, phi = 20,
          shade = 0.1,expand = 0.5, col = "lightblue", r = 5, d =1)
    
    persp(bw_seq, k_seq, AICc,
          main="AICc",
          zlab = "AICc score",
          theta = 25, phi = 20,
          shade = 0.1,expand = 0.5, col = "lightblue", r = 5, d =1)
    
    persp(bw_seq, k_seq, BIC,
          main="BIC",
          zlab = "BIC score",
          theta = 25, phi = 20,
          shade = 0.1,expand = 0.5, col = "lightblue", r = 5, d =1)
    
    persp(bw_seq, k_seq, RECP,
          main="RECP",
          zlab = "RECP score",
          theta = 25, phi = 20,
          shade = 0.1,expand = 0.5, col = "lightblue", r = 5, d =1)
  }
  sel          <-new.env()
  sel$pairs    <-pairs            
  sel$GCV      <-GCV               
  sel$AICc     <-AICc         
  sel$BIC      <-BIC            
  sel$RECP     <-RECP               
  
  return(sel)
}