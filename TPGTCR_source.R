
library(hypergeo) 
library(ggplot2)

E_i <- function(yi, theta){ 
  rt <- theta[1]; alpht <- theta[2]; betat <- theta[3]     
  at <- alpht + 1/betat   
  bit <-  alpht + abs(yi)^betat / ( 2*(1+rt*sign(yi))^betat )  
  d1 <- at / bit  
  d2 <- digamma(at) - log(bit)    
  return( list(Ezi = d1, Elogzi = d2, at = at, bit = bit) )     
}   
 
SGT_rab_EM <-function(sample){
  log.likelihood_value <-c() 
  n_0=length(sample) ; delta=1e-7
  r_star <- (sum(sample>=mean(sample)) - sum(sample<mean(sample) ) ) / n_0   
  alph_star <-1.5 ; beta_star <-  2
  theta_star=c(r_star,alph_star,beta_star)
  theta=theta_star
  itr=0  
  
  Tr <- 0 ; Ta <- 1.5; Tb <- 2  
  likelihood_values <- numeric() 
  ####### E #################  
  while(1){
    theta0 = theta
    likelihood_value0 <- sum(log(fx(x=sample, r=theta0[1], a = theta0[2], b = theta0[3] )))
    Eresult <- E_i(sample, theta0) 
    Ei1 <- Eresult$Ezi ;  Ei2 <- Eresult$Elogzi ; at <- Eresult$at ; bit <- Eresult$bit
    
  ###### M ###############
    eq1 <- function(x){
      eta <- 1/theta0[3]
      ret <- sum( abs(sample)^(1/eta)*sign(sample) / bit / (1 + x * sign(sample) )^(1+1/eta) )
      return(abs(ret)^2)
    }
    
    theta[1] <- optim(theta0[1], eq1, lower = -0.9, upper= 0.9, method = 'L-BFGS-B' )$par
    Eresult <- E_i(sample, theta)
    at <- Eresult$at ; bit <- Eresult$bit
    
    eq2 <- function(x){
      ret <- log(x) + 1 -digamma(x) + digamma(at) - mean( log(bit) ) - at / n_0 * sum( 1/ bit)
      return((ret)^2)
    } 
    theta[2]  <- optim(theta0[2], eq2, lower = 1.5 ,upper= 100, method = 'L-BFGS-B' )$par
    Eresult <- E_i(sample, theta) 
    at <- Eresult$at ; bit <- Eresult$bit 
    
    eq3 <- function(x){
      M_i <- abs(sample) / (1+ theta0[1]*sign(sample) ) 
      ret <- x +log(2) + digamma(1/x) - digamma(at) + mean( log(bit) ) - at *x^2 / 2 / n_0 * sum( M_i^x / bit *log(M_i)  )
      return((ret)^2) 
    }  
    theta[3] <- optim(Tb-0.5, eq3, lower = 1.5, upper= 5, method = 'L-BFGS-B' )$par  
   
    likelihood_value1= sum(log(fx(x=sample, r=theta[1], a = theta[2], b = theta[3] )))  
    itr<-itr+1 
    
    if (abs(likelihood_value1-likelihood_value0) <= delta|itr>=500){
      break
    }
    likelihood_values[itr] <- likelihood_value0
  }
  plot(likelihood_values)     
  log.likelihood_value <- sum(log(fx(x=sample, r=theta[1], a = theta[2], b = theta[3] )))
  result=list( theta_MLE= theta,
               log.likelihood_value=log.likelihood_value )
  result
}


sgt.logL <-function(theta,y,X,vi, ratio,small1)    
{      
  n=length(y)            
  p=dim(X)[2]   
  beta<-signif(theta[1:p],digits=7)
  mu=X%*%beta  
  sigma<-signif(theta[p+1],digits=7) 
  r<-signif(theta[p+2],digits=7) 
  a<-signif(theta[p+3],digits=7)
  b<-signif(theta[p+4],digits=7)       
  z1<- (vi-mu)/sigma    
  
  if(ratio!=0){    
    big1 <- setdiff(seq(1,n,by=1), small1)      
    vero<-signif( c( Fx_par(z1[small1],r,a,b), fx(z1[big1],r,a,b)/sigma ),digits=7 )  
    return(sum(log(vero)))  
  }  
  if(ratio==0){ 
    vero<-signif( fx(z1,r,a,b)/sigma,digits=7 )  
    return(sum(log(vero)))   
  } 
}  
 
sgt.find.log <- function(theta){
  d <- sgt.logL(theta,y,X,vi, ratio, small1) 
  return(-d)  
}

fsgt_five<-function(x,mut,sigmat,rt,alpht,betat){
  chang<-betat/(2*sigmat*(2*alpht)^(1/betat)*beta(alpht,1/betat))
  d<-chang*(1+abs(x-mut)^betat/(2*alpht*sigmat^betat*(1+rt*sign(x-mut))^betat))^(-alpht-1/betat)
  return(d) 
}   

Fsgt_five <- function(x,mut,sigmat,r,a,b){
  Fx_par( (x-mut)/sigmat, r, a, b) 
} 
 
Fx_par <-function(x,r,a,b){  
  d<-c() 
  for (ii in 1:length(x)){
    if(x[ii]<0 ) { 
      ux<-(1+abs(x[ii])^b/(2*a*(1+r*sign(x[ii]))^b ) )^(-1) 
      Iu<-pbeta(ux,a,1/b) 
      d[ii]<- (1-r)/2*Iu } else {
        ux<-(1+abs(x[ii])^b/(2*a*(1+r*sign(x[ii]))^b ) )^(-1) 
        Iu<-pbeta(ux,a,1/b) 
        d[ii]<- 1 - (1+r)/2*Iu 
      }
  } 
  return(d)   
}  

rsgt <- function(n, location=0, scale=1, r=0, a=1,b=2, dp=NULL) 
{ 
  if(!is.null(dp)) {      
    if(!missing(r))   
      stop("You cannot set both component parameters and dp")
    location <- dp[1] 
    scale <- dp[2] 
    r <- dp[3]
    a  <- dp[4] 
    b  <- dp[5] 
  }   
  if(scale<0 | a<0 | b<0) {
    stop("Parameter scale must be positive") 
  } 
  if(abs(r)>1) {
    stop("Parameter r must be bounded") 
  }   
  p=(1+r)/2        
  w1<-rbinom(n,1,p)*2-1+r  
  y1<-rgamma(n,shape=1/b,rate=1)  
  z1<-rgamma(n,shape=a,rate=a)
  error<- 2^(1/b)*w1*y1^(1/b)*z1^(-1/b)  
  y<-location+scale*error 
  return(y)         
}


EXI <- function(i, n, r, alph, beta) {
  intfzheng <- function(x) {
    ux <- (1 + abs(x)^beta / (2 * alph * (1 + r * sign(x))^beta))^(-1)
    Iu <- pbeta(ux, alph, 1 / beta)
    Fzhengx <- 1 + (-0.5 * r - 0.5) * Iu
    d <- x * Fzhengx^(i - 1) * (1 - Fzhengx)^(n - i) * fx(x, r, alph, beta)
    return(d)
  }
  
  intffu <- function(x) {
    ux <- (1 + abs(x)^beta / (2 * alph * (1 + r * sign(x))^beta))^(-1)
    Iu <- pbeta(ux, alph, 1 / beta)
    Ffux <- 0.5 * (1 - r) * Iu
    d <- x * Ffux^(i - 1) * (1 - Ffux)^(n - i) * fx(x, r, alph, beta)
    return(d)
  }
  
  i1 <- integrate(intfzheng, 0, Inf)$value
  i2 <- integrate(intffu, -Inf, 0)$value
  d <- i * choose(n, i) * (i1 + i2)
  return(d) 
}

fx <- function(x, r, alph, beta){
  con1 <- beta / (2^(1 + 1 / beta) * alph^(1 / beta) * beta(alph, 1 / beta))
  d <- (1 + abs(x)^beta / (2 * alph * (1 + r * sign(x))^beta))^(-alph - 1 / beta)
  return(con1 * d)                
}  

Fx <- function(x, r, alph, beta) {
  ux <- (1 + abs(x)^beta / (2 * alph * (1 + r * sign(x))^beta))^(-1)
  Iu <- pbeta(ux, alph, 1 / beta)
  d <- ifelse(x >= 0, 1, 0) + (-0.5 * r + ifelse(x >= 0, -0.5, 0.5)) * Iu
  return(d)
}

cdf_inv<- function(F,r,a,b){
  F0 <- (1-r)/2  ; d<-c()
  for (ii in 1:length(F)) {
    if(F[ii]<=F0 ) {
      yx<- qbeta( (2*F[ii]/(1-r)), a, 1/b )
      d[ii]<- ( (1/yx-1)*2*a )^(1/b)*( 1+r*sign( F[ii]-(1-r)/2 ) ) * sign( F[ii]-(1-r)/2 )  } else {
        yx<- qbeta( 2*(1 - F[ii])/(1+r),a, 1/b )
        d[ii]<- ( (1/yx-1)*2*a )^(1/b)*( 1+r*sign( F[ii]-(1-r)/2 ) ) * sign( F[ii]-(1-r)/2 )
      } }
  return(d) 
} 

h0deri <- function(z, r, alph, beta) {
  con1 <- -sign(z) * beta * (alph * beta + 1) / (2 * (2 * alph)^(1 / beta + 1) * (1 + sign(z) * r)^beta * beta(alph, 1 / beta))
  ux <- (1 + abs(z)^beta / (2 * alph * (1 + r * sign(z))^beta))^(-1)
  Iu <- pbeta(ux, alph, 1 / beta) 
  Fx_val <- ifelse(z >= 0, 1, 0) + (-0.5 * r + ifelse(z >= 0, -0.5, 0.5)) * Iu
  fxderi_val <- con1 * abs(z)^(beta - 1) * ux^(alph + 1 / beta + 1) 
  result <- fxderi_val / Fx_val - (fx(z, r, alph, beta) / Fx_val)^2
  return(result)
}

h1deri <- function(z, r, alph, beta) {
  abs_z <- abs(z)
  abs_z_beta <- abs_z^beta
  sign_z <- sign(z)
  uz <- (1 + abs_z_beta / (2 * alph * (1 + r * sign_z)^beta))^(-1)
  d1 <- abs_z^(beta - 2) * uz * (beta * uz - 1)
  con1 <- beta / (2 * alph * (1 + sign_z * r)^beta) 
  result <- con1 * d1
  return(result)
} 

h0<-function(z, r, alph, beta){ 
  fx(z, r, alph, beta)/Fx(z, r, alph, beta)  
}   

h1<-function(z, r, alph, beta){
  con1<-sign(z)*beta*abs(z)^(beta-1)/( 2*alph*(1+sign(z)*r)^beta ) 
  uz<-(1+abs(z)^beta/(2*alph*(1+r*sign(z))^beta ))^(-1) 
  return(con1*uz) 
}  

a0i <- function(i, n, r, alph, beta){ 
  exi_val <- EXI(i, n, r, alph, beta)
  h0_val <- h0(exi_val, r, alph, beta)
  h0deri_val <- h0deri(exi_val, r, alph, beta)
  result <- h0_val - h0deri_val * exi_val
  return((result))
}  

a1i <- function(i, n, r, alph, beta) {
  exi_val <- EXI(i, n, r, alph, beta)
  h1_val <- h1(exi_val, r, alph, beta) 
  h1deri_val <- h1deri(exi_val, r, alph, beta)
  result <- h1_val - h1deri_val * exi_val
  return(result)
}

b0i <- function(i, n, r, alph, beta) {
  exi_val <- EXI(i, n, r, alph, beta)
  h0deri_val <- h0deri(exi_val, r, alph, beta)
  result <- -h0deri_val
  return(result) 
}

b1i <- function(i, n, r, alph, beta) {
  exi_val <- EXI(i, n, r, alph, beta)
  h1deri_val <- h1deri(exi_val, r, alph, beta)
  result <- -h1deri_val
  return(result)
}

EX2I<-function(i, n, r, alph, beta){
  intfzheng <- function(x) {
    ux <- (1 + abs(x)^beta / (2 * alph * (1 + r * sign(x))^beta))^(-1)
    Iu <- pbeta(ux, alph, 1 / beta)
    Fzhengx <- 1 + (-0.5 * r - 0.5) * Iu
    d <- x^2 * Fzhengx^(i - 1) * (1 - Fzhengx)^(n - i) * fx(x, r, alph, beta)
    return(d)
  } 
  intffu <- function(x) {
    ux <- (1 + abs(x)^beta / (2 * alph * (1 + r * sign(x))^beta))^(-1)
    Iu <- pbeta(ux, alph, 1 / beta)
    Ffux <- 0.5 * (1 - r) * Iu
    d <- x^2 * Ffux^(i - 1) * (1 - Ffux)^(n - i) * fx(x, r, alph, beta)
    return(d)
  } 
  i1<-integrate( intfzheng,0,Inf )$value
  i2<-integrate( intffu,-Inf,0 )$value
  d<-i*choose(n,i)*(i1+i2) 
  return(d)  
}  

HessianaQ <- function(X,zi,theta,small1){
  n <- nrow(X)    
  p <- ncol(X)     
  beta0 <- signif(theta[1],digits=7)  
  beta1 <- signif(theta[2:(p+1)],digits=7)      
  sigma <- signif(theta[p+2],digits=7)  
  r    <-signif(theta[p+3],digits=7) 
  alph <-signif(theta[p+4],digits=7)
  beta <-signif(theta[p+5],digits=7)    
  
  small1  <- sort(small1) 
  big1 <- setdiff(seq(1,n,by=1), small1)   
  ranks <- floor(rank(zi) )     
  
  tti <-  sapply(1:n, function(i){ EXI(i, n, r, alph, beta) } )  
  eti2 <- sapply(1:n, function(i){ EX2I(i, n, r, alph, beta) } )  
  weight <- 1/ sigma^2   
  I11 <- weight* sum( sapply((ranks[small1]) , b0i, n=n, r=r, alph=alph, beta=beta) ) - (Ta+1/Tb) * sum( sapply((ranks[big1]), b1i, n=n, r=r, alph=alph, beta=beta) ) 
  I12 <- c() ; L2 <- L3 <- c() ; sub_I23 <- c() 
  for (j in 1:p) { 
    new <- weight* ( sum( sapply((ranks[small1]) , b0i, n=n, r=r, alph=alph, beta=beta)*X[small1,j] ) - (Ta+1/Tb) * sum( sapply((ranks[big1]), b1i, n=n, r=r, alph=alph, beta=beta)*X[big1,j]  ) ) 
    I12 <- c(I12, new)      
    new_L2 <- ( sum( ( sapply((ranks[small1]) , a0i, n=n, r=r, alph=alph, beta=beta)-sapply((ranks[small1]) , b0i, n=n, r=r, alph=alph, beta=beta)*tti[(ranks[small1])] )*X[small1,j] ) - (Ta+1/Tb) * sum(
      ( sapply((ranks[big1]) , a1i, n=n, r=r, alph=alph, beta=beta)- sapply((ranks[big1]) , b1i, n=n, r=r, alph=alph, beta=beta)*tti[(ranks[big1])] ) *X[big1,j]  ) ) 
    L2 <- c(L2, new_L2) 
    newsub_I23 <- ( sum( sapply((ranks[small1]) , b0i, n=n, r=r, alph=alph, beta=beta)*(X[small1,j])*tti[(ranks[small1])] ) - (Ta+1/Tb) * sum( sapply((ranks[big1]), b1i, n=n, r=r, alph=alph, beta=beta)*(X[big1,j])*tti[(ranks[big1])] ) ) 
    sub_I23 <- c(sub_I23, newsub_I23)  
  }  
  I23 <- -weight*L2 + weight*sub_I23   
  
  I13 <- weight* ( sum( sapply((ranks[small1]) , b0i, n=n, r=r, alph=alph, beta=beta)*tti[(ranks[small1])] ) - (Ta+1/Tb) * sum( sapply((ranks[big1]), b1i, n=n, r=r, alph=alph, beta=beta)*tti[(ranks[big1])]  ) ) 
  I22 <- weight* ( t(X[small1,])%*% diag( sapply((ranks[small1]) , b0i, n=n, r=r, alph=alph, beta=beta) ) %*%X[small1,] - (Ta+1/Tb) * t(X[big1,])%*% diag( sapply((ranks[big1]) , b1i, n=n, r=r, alph=alph, beta=beta) ) %*%X[big1,]  )  
  I33 <- -weight* ( n-length(small1) + sum( 2*sapply((ranks[small1]) , a0i, n=n, r=r, alph=alph, beta=beta)*tti[(ranks[small1])] - 3*sapply((ranks[small1]) , b0i, n=n, r=r, alph=alph, beta=beta)*eti2[(ranks[small1])]) - 
                      (Ta+1/Tb) * sum(  2*sapply((ranks[big1]) , a1i, n=n, r=r, alph=alph, beta=beta)*tti[(ranks[big1])] - 3*sapply((ranks[big1]) , b1i, n=n, r=r, alph=alph, beta=beta)*eti2[(ranks[big1])] )  )
  
  MatrizQ<-matrix(0,nrow=(p+2),ncol=(p+2) )    
  MatrizQ[1,1] <- I11 
  MatrizQ[1,2:(p+1)] <- I12 
  MatrizQ[1,(p+2)] <- I13  
  MatrizQ[2:(p+1),2:(p+1)] <- I22  
  MatrizQ[2:(p+1),(p+2)] <- I23    
  MatrizQ[(p+2),(p+2)] <- I33   
  MatrizQ[lower.tri(MatrizQ)] <- t(MatrizQ)[lower.tri(MatrizQ)]  
  return(MatrizQ)  
} 

MSE <- function(est ){ 
  d<- ( sum(( (est[2:(p+1)] - Ttheta[2:(p+1)]) )^2)  )   
  return( d )  }   

Rbias <- function(est,theta){ 
  ret <- mean( abs(est-theta)/theta )    
  return( ret )   
}   

Bias_abs <- function(est,theta){ 
  ret <- mean( abs(est-theta) )    
  return( ret )   
}   

Bias <- function(est,theta){ 
  ret <- mean( (est-theta) )     
  return( ret )   
}     

RMSE<-function(est,theta){
  d<-sqrt( mean((est-theta)^2)  )  
  return( d )  
}      

max_vec <- function(x,y){
  iftrue <- { x<y }  
  return( y^iftrue * x^(1-iftrue)  )   
}   

dSGT <- function(x,mu,sigma,r,a,b){
  if( is.null(x) ){ return(c())  }else{
    fsgt_five(x,  mu,sigma,r,a,b )    
  } 
}  
pSGT <-function(q,mu,sigma,r,a,b){ 
  Fsgt_five(q, mu,sigma,r,a,b)    
}     
qSGT <-function(p,mu,sigma,r,a,b){ 
  ret <- mu + sigma * cdf_inv(p,r,a,b )    
  return(ret)  
}   

compute_rmi <- function(zi, rho, r, alph, beta) {
  rmi <- rep(NA, length(y)) 
  Syi <- 1- Fx_par( zi, r, alph, beta  )  
  for (i in 1:length(y)) {
    if (rho[i] == 1) {
      rmi[i] <- log(Syi[i]) 
    } else {
      rmi[i] <- 1 + log(Syi[i])   
    }
  }
  return(rmi)
}

compute_rdi <- function(rmi, rho) {
  rdi <- sign(rmi) * sqrt(-2 * (rmi + (1 - rho) * log(1 - rho - rmi))) 
  return(scale(rdi))  
}    


 