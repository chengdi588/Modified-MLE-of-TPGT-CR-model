 

library(CensRegMod)    
library(censReg)   
library(SMNCensReg) # Garay 2017    
data(wage.rates)      

##  data ----------------------------------------------   

data0 <- wage.rates   
y <- data0$wage ; vi <- ci <- y          

x1 <- data0$age   
x2 <- data0$educ  
x3 <- data0$kidslt6 # kids < 6 years     
x4 <- data0$kidsge6 # kids 6-19          
X <- cbind(1, x1, x2, x3, x4)    
colnames(data0)       

cc <- rho <- (data0$hours==0) + 0   
p <- ncol(X)-1    
n <- nrow(X) ; small1 <- which(cc==1) ; ratio <- sum(cc) / n  ; ratio   
Normal <- em.cens(cc, -X, -y,dist="Normal",diagnostic=FALSE ) 
TCRmodel <- em.cens(cc, -X, -y,dist="T", nu=2, diagnostic=FALSE )     
   
mu <- X%*%Normal$betas    
erro <- -1*(vi - mu) / sqrt(TCRmodel$sigma2)    
RAB <- SGT_rab_EM(erro)$theta_MLE   ; RAB   
 
theta_ini <- c(Normal$betas, sqrt(Normal$sigma2), RAB)    
Max_log_like <- optim(par=theta_ini, fn=sgt.find.log, method = "L-BFGS-B",control = list(maxit=500),
                      lower= c( theta_ini[1:(p+1)]- rep(5,(p+1)),0.001, -0.9, 1.5, 1.7 ) ,
                      upper= c( theta_ini[1:(p+1)]+ rep(5,(p+1)), 20, 0.9, 10, 10 )  )  

theta_best1 <- Max_log_like $par  
RAB <- theta_best1[ ( length(theta_best1)-2): length(theta_best1)]  
RAB 

# -------------- solve MMLE -------------        
{   
  Tr <- r <- RAB[1] ; Ta <- alph <- RAB[2]  ; Tb <- beta <- RAB[3]    
  sigma0 <- sqrt( TCRmodel$sigma2)  
  big1 <- setdiff(seq(1,n,by=1), small1)  
  mu <- X%*%(TCRmodel$betas)    
  
  erro <- (vi - mu)  
  zi <- as.numeric( erro/(sigma0) )    
  ranks <- floor(rank(zi))    
  small_ranks <- ranks[small1]  
  big_ranks <- ranks[big1]  
  
  a_star <- sum(sapply(small_ranks, a0i, n=n, r=r, alph=alph, beta=beta)) -   
    (Ta + 1/Tb) * sum(sapply(big_ranks, a1i, n=n, r=r, alph=alph, beta=beta))  
  
  b_star <- sum(sapply(small_ranks, b0i, n=n, r=r, alph=alph, beta=beta)) -   
    (Ta + 1/Tb) * sum(sapply(big_ranks, b1i, n=n, r=r, alph=alph, beta=beta))  
 
  y_bar <- ( sum( sapply( small_ranks, b0i, n=n, r=r, alph=alph, beta=beta)*vi[small1] ) - (Ta+1/Tb) * sum( 
    sapply(big_ranks, b1i, n=n, r=r, alph=alph, beta=beta) *vi[big1] ) ) / b_star  

  x_bar_j <- numeric(p) 
  for (j in 1:p){ 
    x_bar_j[j] <- (sum(sapply(small_ranks, b0i, n=n, r=r, alph=alph, beta=beta) * X[small1, j]) - 
                     (Ta + 1 / Tb) * sum(sapply(big_ranks, b1i,n=n, r=r, alph=alph, beta=beta) * X[big1, j])) / b_star
  } 
 
  x_bar <- matrix(x_bar_j, nrow = p, ncol = 1)
  
  a0_vec <- sapply(1:n , a0i,n=n, r=r, alph=alph, beta=beta)[ranks] ;  a1_vec <- sapply(1:n , a1i,n=n, r=r, alph=alph, beta=beta)[ranks]
  b0_vec <- sapply(1:n , b0i,n=n, r=r, alph=alph, beta=beta)[ranks] ;  b1_vec <- sapply(1:n , b1i,n=n, r=r, alph=alph, beta=beta)[ranks]    
  M <- X[,-1] - kronecker( matrix(1, n, 1) ,t(x_bar)) 
  B1 <- t(M)%*% ( rho*a0_vec  ) - (Ta+1/Tb) * t(M)%*% ( (1-rho)*a1_vec  )       
  B2 <- t(X[,-1]) %*% ( rho*b0_vec*( ci-y_bar)  ) - (Ta+1/Tb) * t(X[,-1])%*%(diag(n)-diag(rho)) %*% ( b1_vec*(y-y_bar) ) 
  D <- t(M) %*% diag(rho*b0_vec) %*% M -  (Ta+1/Tb) * t(M)%*% diag(b1_vec*(1-rho) ) %*% M  
  K <- solve(D)%*%B2 ;  L <- -solve(D)%*%B1      
  
  ## σ    
  ncens <- floor(n*ratio)    
  J0 <- n-ncens   
  J1 <- t(rho*a0_vec) %*% (ci-y_bar + M%*%K ) -(Ta+1/Tb) *as.numeric( t((1-rho)*a1_vec) %*% (vi-y_bar+ M%*%K ) )
  J2 <- t(rho*b0_vec) %*% (ci-y_bar - M%*%K )^2 -(Ta+1/Tb) *as.numeric( t((1-rho)*b1_vec) %*% (vi-y_bar- M%*%K )^2 )
  
  delta.sigma<- J1^2 + 4*J0*J2       
  sigma_hat <- ( -J1 + sqrt(delta.sigma)) / ( 2*sqrt(J0*(J0-(ncol(X))) ) )  
  sigma_hat <- as.numeric(sigma_hat) 
  
  beta1_hat <- K - L%*%sigma_hat    
  beta0_hat <- y_bar - t(x_bar)%*%beta1_hat - a_star/b_star *sigma_hat         
  theta_hat <- c(beta0_hat, beta1_hat, sigma_hat) 
  
  theta_1 <- c(theta_hat, RAB)  
  SGTCR_like <- sgt.logL(theta_1, (y), X, vi, ratio, small1)  
  print(SGTCR_like)          
}   


## ------  SE -------------------------------------------------------------- 

zi_func <- function(vi,XX,theta_1){   
  beta0<-theta_1[1] ; beta1<-theta_1[2:(p+1)] ; sigma <- theta_1[p+2]      
  r<-signif(theta_1[p+3],digits=7)
  alph<-signif(theta_1[p+4],digits=7)
  beta<-signif(theta_1[p+5],digits=7)   
  mu= XX%*%matrix( c(beta0, beta1), ncol=1 )   
  
  zi <- (vi-mu) / sigma      
  return(zi) 
}   
zi_value <- zi_func(vi,X,theta_1)   

 
H <- HessianaQ(X[,-1], zi_value, theta_1, small1)       
se <-  sqrt( diag( solve(H) ) ) ; se    

max_like <- sgt.logL(theta_1, y, X, vi, ratio, small1)    
log0 <- round(max_like, 4)   
k_num = ncol(X)+1    
AIC<- -2*log0+2*k_num    
BIC<- -2*log0+log(n)*k_num   
EDC<- -2*log0+0.2*sqrt(n)*k_num  
N <- (n + 2)/24     
ABIC <- -2* log0 + k_num*log(N)  


## ------ 1. Compare: AIC BIC EDC -------------------------------------------------------------- 
# SN-CR (2018) 
SMSN_model <- SAEM_EST(y,X,cc,cens="left" ,LS=NULL,precisao=0.0001,MaxIter=500, M=20, pc=0.35, nu=NULL,show.convergence="TRUE",
                       dist="SN",  nu.fixed=FALSE)
 

col1 <- round(matrix(c(theta_1, log0, AIC, BIC , EDC), ncol=1 ), 4 )      
col2 <- matrix(se, ncol= 1)    

col3 <- round(matrix(c(TCRmodel$betas, sqrt(TCRmodel$sigma2), TCRmodel$nu, TCRmodel$logver, TCRmodel$AIC, TCRmodel$BIC, TCRmodel$EDC), ncol=1 ), 4 )    
col4 <- matrix(TCRmodel$SE, ncol= 1)    

col5 <- round(matrix(c(Normal$betas, sqrt(Normal$sigma2), Normal$logver, Normal$AIC, Normal$BIC, Normal$EDC), ncol=1 ), 4 )    
col6 <- matrix(Normal$SE, ncol= 1)    

col7 <- round(matrix(SMSN_model$teta, ncol=1), 4 )  
col8 <- round(matrix(SMSN_model$EP, ncol=1), 4) 

max_length <- max(length(col1), length(col2), length(col3), length(col4), length(col5), length(col6), length(col7), length(col8))

result <- array(NA, dim = c(max_length, 8)) 
 
result[1:length(col1), 1] <- col1
result[1:length(col2), 2] <- col2
result[1:length(col3), 3] <- col3
result[1:length(col4), 4] <- col4
result[1:length(col5), 5] <- col5
result[1:length(col6), 6] <- col6
result[1:length(col7), 7] <- col7
result[1:length(col8), 8] <- col8 
rownames(result) <- paste("Row", 1:max_length) 
print(result) 


## ------------- 2. Test ---------------------------------------------   

log_normal <-  Normal$logver 
log_t <- TCRmodel$logver 
log_sncr <- round(SMSN_model$likelihood, 4)    
-2*( log_normal - log0 )   
-2*( log_t - log0 )  
-2*( log_sncr - log0 )  

## ------------- 3. Plot ---------------------------------------------   

library(fitdistrplus)  
mu0 <- X[cc==0,]%*%matrix(  theta_1[1:5], ncol=1 )   
fit.sgt <- fitdist( y[cc==0],"SGT", method = "mle", start = list( mu=1.3, 
                                                                  sigma=theta_1[6] ,
                                                                  r= theta_1[7],
                                                                  a= theta_1[8],
                                                                  b= theta_1[9] ) )   

summary(fit.sgt)   
dev.new()   
par(mar= (c(7,9,4,1) + 0.1), mgp = c(5.5, 2, 0) )       

denscomp(fit.sgt, addlegend=FALSE, main = "", 
         xlab = "y", fitlwd = 5, cex.axis = 3., 
         cex.lab = 3., cex.main = 3. )  

ppcomp(fit.sgt, addlegend=FALSE, main = "" , 
       lwd = 3, cex.axis = 3.,  
       cex.lab = 3., cex.main = 3. )     

# Residual   

rmi <- compute_rmi(zi_value, rho, r, alph, beta)   
rdi <- compute_rdi(rmi, rho)    
 
plot(rdi, ylim=c(-3.2,3.2),
     xlab = "Index", ylab = "Modified deviance residual",
     lwd = 3, cex.axis = 3., cex.axis = 3., 
     cex.lab = 3., cex.main = 3. )  
abline(h = -3, lty=3,  lwd=3 )
abline(h = 3, lty=3,  lwd=3)


# T residual  

# TCR_mu <- X%*%(TCRmodel$betas)    
# TCR_resid <- as.numeric( (vi - mu)/sqrt(TCRmodel$sigma2) )  
#                   
# plot(TCR_resid, ylim=c(-3.2, 8.2),
#      xlab = "Index", ylab = "Residuals of T-CR",
#      lwd = 3, cex.axis = 3., cex.axis = 3., 
#      cex.lab = 3., cex.main = 3. )  
# abline(h = -3, lty=3,  lwd=3 )
# abline(h = 3, lty=3,  lwd=3)  


