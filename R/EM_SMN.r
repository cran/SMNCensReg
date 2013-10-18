################################################################################
###     EM algorithm for the Normal/independent Censored Regression Model    ###
################################################################################
### type: NI distribution to be consideraded
### 1.- Normal
### 2.- T (student-t)
### 3.- PearsonVII  ---- in this case, delta=1 is fixed
### 4.- Slash 
### 5.- NormalC (contaminated Normal) ---- in this case, nu is fixed
################################################################################
### cens: type of censoring   
###  1: left  (Y_i <= Ki)
###  2: right   (Y_i >= Ki)
###  3: Interv (Y_i >= Ki ou Y_i <= Ki)
################################################################################

###### EM algorithm

EM.Cens.Int <- function(cc, x,y,LS=NULL,nu=3,delta=0,cens="1",type="T",error=0.0001, iter.max=200)
{
  x <- as.matrix(x)
  y <- as.vector(y)
  p <- ncol(x)
  n <- nrow(x)

  reg <- lm(y ~ -1 + x[,1:p])
  betas <- as.vector(coefficients(reg),mode="numeric")
  sigma2 <- sum((y-x%*%betas)^2)/(n-p)

  ERRO<-1e-6
  TOLERANCIA<-1e-6
  MAX_NU<-150
  MIN_NU <- 1.01
  count <- 0
  criterio <- 1
  iter.max = 150
  SE1 <- matrix(0,1,p)
  res <- NULL
  Lim1  <- Lim2 <- c()

  if (cens=="1")
  {
  Lim1 <- rep(-Inf,n)
  Lim2 <- y
  }
  if (cens=="2")
  {
  Lim1 <- y
  Lim2 <- rep(Inf,n)
  }
  if (cens=="3")
  {
  Lim1 <- y
  Lim2 <- LS
  }
  if (type == "T")
  {
  param  <- matrix(c(betas,sigma2,nu),ncol=1)      
  while(criterio > error || count<=iter.max)
     {
         count <- count + 1
            mu <- x%*%(betas)
     # Passo E
     
     NCensEUY <- NCensurEsperUY(y,mu,sigma2,nu,0,type=type)    
        u0 <- NCensEUY$EUY0
        u1 <- NCensEUY$EUY1
        u2 <- NCensEUY$EUY2
     
      CensEUY <- CensEsperUY1(mu[cc==1],sigma2=sigma2,nu=nu,delta=delta,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type=type, cens=cens)
     u0[cc==1]<- CensEUY$EUY0
     u1[cc==1]<- CensEUY$EUY1
     u2[cc==1]<- CensEUY$EUY2
            S1 <- matrix(0, n)
            S2 <- matrix(0, n)
            S3 <- matrix(0, n)
 			
         
         suma1 <- t(t(x)%*%u1)
         
         u0.matriz <- Diagonal(n,as.numeric(sqrt(u0)))
         xnovo <- as.matrix(u0.matriz%*%x)
         suma2 <- t(xnovo)%*%xnovo
         
       
       betas <- t(solve(suma2)%*%t(suma1))
   		sigma2 <- sum(u2-2*u1*mu+mu^2*u0)/n
 			 auxf0 <- (y-x%*%t(betas))/sqrt(sigma2)
  			auxf <- (Lim1-x%*%t(betas))/sqrt(sigma2)
        auxf1 <- (Lim2-x%*%t(betas))/sqrt(sigma2)
         
         ft1 <- function(nu){sum(log(dt(auxf0[cc==0],df=nu)/sqrt(sigma2)))+ sum(log(pt(auxf0[cc==1],df=nu)))}
         
         ft2 <- function(nu){sum(log(dt(auxf0[cc==0],df=nu)/sqrt(sigma2)))+ sum(log(pt(-auxf0[cc==1],df=nu)))}
         
         ft3 <- function(nu){sum(log(dt(auxf0[cc==0],df=nu)/sqrt(sigma2)))+ sum(log(pt(auxf1[cc==1],df=nu)-pt(auxf[cc==1],df=nu)))}
         
         if (cens=="1"){	  
           nu <- optimize(f=ft1, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum 
           logver=ft1(nu) 
         }
         if (cens=="2"){
           nu <- optimize(f=ft2, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum 
           logver=ft2(nu) 
         }
         if(cens=="3"){
    		  nu <- optimize(f=ft3,interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,tol = TOLERANCIA, maximum = TRUE)$maximum 
    		  logver=ft3(nu)
          }
         tetas <- matrix(c(as.vector(betas),sigma2,nu),ncol=1)      
    criterio <- sqrt(t(tetas-param)%*%(tetas-param))
       betas <- t(betas)
       param <- tetas
       }
    AIC <- (-2)*logver + 2*(p+2)   ## p (Colunas de X) + Sigma2 + nu 
    BIC <- (-2)*logver + (p+2)*log(n)
    EDC <- (-2)*logver + (p+2)*0.2*sqrt(n)    
    
    SE1 <- IM.Cens(cc, x,y,LS, BETA=betas,SIGMA2=sigma2,NU=nu,delta,cens,type)
    obj.out <- list(betas = betas, sigma2 = sigma2, nu = nu, delta=delta,logver = logver, count=count, AIC=AIC, BIC=BIC, EDC=EDC, SE=SE1)
   class(obj.out) <- type
   return(obj.out)
  }
   
   
  if (type == "PearsonVII")
  {
	delta<-1
  param  <- matrix(c(betas,sigma2,nu,delta),ncol=1)  
      
    u0 <- u1 <- u2 <- c()
   while(criterio > error || count<=iter.max)
     {
         count <- count + 1
            mu <- x%*%(betas)
        # Passo E

    NCensEUY <- NCensurEsperUY(y,mu,sigma2,nu,delta,type=type)    
           u0 <- NCensEUY$EUY0
        	 u1 <- NCensEUY$EUY1
           u2 <- NCensEUY$EUY2
     
      CensEUY <- CensEsperUY1(mu[cc==1],sigma2=sigma2,nu=nu,delta=delta,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type=type, cens=cens)
     u0[cc==1]<- CensEUY$EUY0
     u1[cc==1]<- CensEUY$EUY1
     u2[cc==1]<- CensEUY$EUY2
            S1 <- matrix(0, n)
            S2 <- matrix(0, n)
            S3 <- matrix(0, n)
         suma1 <- t(t(x)%*%u1)
         
         u0.matriz <- Diagonal(n,as.numeric(sqrt(u0)))
         xnovo <- as.matrix(u0.matriz%*%x)
         suma2 <- t(xnovo)%*%xnovo
         
        betas <- t(solve(suma2)%*%t(suma1))
 		sigma2 <- sum(u2-2*u1*mu+mu^2*u0)/n
		 auxf0 <- (y-x%*%t(betas))/sqrt(sigma2)
 			auxf <- (Lim1-x%*%t(betas))/sqrt(sigma2)
     auxf1 <- (Lim2-x%*%t(betas))/sqrt(sigma2)
         
         fp1 <- function(nu){sum(log(dPearsonVII(auxf0[cc==0],0,1,nu,delta)/sqrt(sigma2)))+ sum(log(PearsonVII(auxf0[cc==1],0,1,nu,delta)))}
         fp2 <- function(nu){sum(log(dPearsonVII(auxf0[cc==0],0,1,nu,delta)/sqrt(sigma2)))+ sum(log(PearsonVII(-auxf0[cc==1],0,1,nu,delta)))}
         fp3 <- function(nu){sum(log(dPearsonVII(auxf0[cc==0],0,1,nu,delta)/sqrt(sigma2)))+ sum(log(PearsonVII(auxf1[cc==1],0,1,nu,delta)-PearsonVII(auxf[cc==1],0,1,nu,delta)))}
         
         if (cens=="1"){
           nu <- optimize(f=fp1,interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,tol = TOLERANCIA, maximum = TRUE)$maximum 
           logver= fp1(nu)
         }
         
         if (cens=="2"){
           nu <- optimize(f=fp2,interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,tol = TOLERANCIA, maximum = TRUE)$maximum 
           logver= fp2(nu)
         }
         
         if (cens=="3"){
           nu <- optimize(f=fp3,interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,tol = TOLERANCIA, maximum = TRUE)$maximum 
           logver= fp3(nu)
         }
         
   	   tetas <- matrix(c(as.vector(betas),sigma2,nu,delta),ncol=1)      
    criterio <- sqrt(t(tetas-param)%*%(tetas-param))
       betas <- t(betas)
       param <- tetas
       }
    AIC <- (2)*logver + 2*(p+3)   ## p (Colunas de X) + Sigma2 + nu + delta
    BIC <- (2)*logver + (p+3)*log(n)
    EDC <- (2)*logver + (p+3)*0.2*sqrt(n)    
    SE1 <- IM.Cens(cc, x,y,LS, BETA=betas,SIGMA2=sigma2,NU=nu,delta,cens,type)
    obj.out <- list(betas = betas, sigma2 = sigma2, nu = nu, delta=delta,logver = logver, count=count, AIC=AIC, BIC=BIC, EDC=EDC, SE=SE1)
  class(obj.out) <- type
  return(obj.out) 
  }
   
  if (type == "Slash")
  {
  u0 <- u1 <- u2 <- c()
 		param  <- matrix(c(betas,sigma2,nu),ncol=1)                 

   while(criterio > error || count<=iter.max)
     {
         count <- count + 1
            mu <- x%*%(betas)
        # Passo E
  NCensEUY <- NCensurEsperUY(y,mu,sigma2,nu,0,type=type)    
           u0 <- NCensEUY$EUY0
        	 u1 <- NCensEUY$EUY1
           u2 <- NCensEUY$EUY2
     
      CensEUY <- CensEsperUY1(mu[cc==1],sigma2=sigma2,nu=nu,delta=delta,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type=type, cens=cens)
     u0[cc==1]<- CensEUY$EUY0
     u1[cc==1]<- CensEUY$EUY1
     u2[cc==1]<- CensEUY$EUY2

            S1 <- matrix(0, n)
            S2 <- matrix(0, n)
            S3 <- matrix(0, n)
 		
         suma1 <- t(t(x)%*%u1)
         
         u0.matriz <- Diagonal(n,as.numeric(sqrt(u0)))
         xnovo <- as.matrix(u0.matriz%*%x)
         suma2 <- t(xnovo)%*%xnovo
         
         
  		 betas <- t(solve(suma2)%*%t(suma1))
   		sigma2 <- sum(u2-2*u1*mu+mu^2*u0)/n

		 auxf0 <- (y-x%*%t(betas))/sqrt(sigma2)
 			auxf <- (Lim1-x%*%t(betas))/sqrt(sigma2)
     auxf1 <- (Lim2-x%*%t(betas))/sqrt(sigma2)
         
         fs1 <- function(nu){sum(log(dSlash(auxf0[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumSlash(auxf0[cc==1],0,1,nu)))}
         
         fs2 <- function(nu){sum(log(dSlash(auxf0[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumSlash(-auxf0[cc==1],0,1,nu)))}
         
         fs3 <- function(nu){sum(log(dSlash(auxf0[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumSlash(auxf1[cc==1],0,1,nu)-AcumSlash(auxf[cc==1],0,1,nu)))}
         
         if (cens=="1")
         {
           nu <- optimize(fs1, c(1.01,30), tol = TOLERANCIA, maximum = TRUE)$maximum
           logver=fs1(nu)
         }
         if (cens=="2")
         {
           nu <- optimize(fs2, c(1.01,30), tol = TOLERANCIA, maximum = TRUE)$maximum
           logver=fs2(nu)
         }
         if (cens=="3")
         {
           nu <- optimize(fs3, c(1.01,30), tol = TOLERANCIA, maximum = TRUE)$maximum
           logver=fs3(nu)
         }
         
       tetas <- matrix(c(as.vector(betas),sigma2,nu),ncol=1)      
    criterio <- sqrt(t(tetas-param)%*%(tetas-param))
       betas <- t(betas)
        param <- tetas
       }
    AIC <- (-2)*logver + 2*(p+2)   ## p (Colunas de X) + Sigma2 + nu
    BIC <- (-2)*logver + (p+2)*log(n)
    EDC <- (-2)*logver + (p+2)*0.2*sqrt(n)    
    SE1 <- IM.Cens(cc, x,y,LS, BETA=betas,SIGMA2=sigma2,NU=nu,delta=NULL,cens,type)
   obj.out <- list(betas = betas, sigma2 = sigma2, nu = nu, logver = logver, count=count, res=res, AIC=AIC, BIC=BIC, EDC=EDC, SE=SE1)
   class(obj.out) <- type
   return(obj.out) 
  }
  if (type == "NormalC")
  {
   # Inicializa beta e sigma2 com os estimadores de mínimos quadrados
     R1 <- c()
  count <- 0
 param  <- matrix(c(betas,sigma2,nu),ncol=1)      
     u0 <- u1 <- u2 <- c()           
 
   while(criterio > error || count<=iter.max)
     {
         count <- count + 1
            mu <- x%*%(betas)
   # Passo E
   
  NCensEUY <- NCensurEsperUY(y,mu,sigma2,nu=nu,0,type=type)    
           u0 <- NCensEUY$EUY0
        	 u1 <- NCensEUY$EUY1
           u2 <- NCensEUY$EUY2
     
      CensEUY <- CensEsperUY1(mu[cc==1],sigma2=sigma2,nu=nu,delta=delta,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type=type, cens=cens)
     u0[cc==1]<- CensEUY$EUY0
     u1[cc==1]<- CensEUY$EUY1
     u2[cc==1]<- CensEUY$EUY2
            S1 <- matrix(0, n)
            S2 <- matrix(0, n)
            S3 <- matrix(0, n)
 		
         suma1 <- t(t(x)%*%u1)
         
         u0.matriz <- Diagonal(n,as.numeric(sqrt(u0)))
         xnovo <- as.matrix(u0.matriz%*%x)
         suma2 <- t(xnovo)%*%xnovo
 
       betas <- t(solve(suma2)%*%t(suma1))
   		sigma2 <- sum(u2-2*u1*mu+mu^2*u0)/n

		 auxf0 <- (y-x%*%t(betas))/sqrt(sigma2)
 			auxf <- (Lim1-x%*%t(betas))/sqrt(sigma2)
     auxf1 <- (Lim2-x%*%t(betas))/sqrt(sigma2)
         
         fnc1 <- f <- function(nu2){
           a1 <- exp(nu2[1])/(1+exp(nu2[1]))
           a2 <- exp(nu2[2])/(1+exp(nu2[2]))      
           ver <-(sum(log(dNormalC(auxf0[cc==0],0,1,c(a1,a2))/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf0[cc==1],0,1,c(a1,a2)))))
           return(-ver)
         }
         
         f11 <- function(nu){
           ver1 <-(sum(log(dNormalC(auxf0[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf0[cc==1],0,1,nu))))
           return(ver1)
         }
              
         
         fnc2 <- function(nu2){
           a1 <- exp(nu2[1])/(1+exp(nu2[1]))
           a2 <- exp(nu2[2])/(1+exp(nu2[2]))      
           ver <- (sum(log(dNormalC(auxf0[cc==0],0,1,c(a1,a2))/sqrt(sigma2)))+ sum(log(AcumNormalC(-auxf0[cc==1],0,1,c(a1,a2)))))
           return(-ver)
         }
         f12 <- function(nu){
           ver1 <-(sum(log(dNormalC(auxf0[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf0[cc==1],0,1,nu))))
           return(ver1)
         }
         
         
         fnc3 <- function(nu2){
           a1 <- exp(nu2[1])/(1+exp(nu2[1]))
           a2 <- exp(nu2[2])/(1+exp(nu2[2]))      
           ver <-(sum(log(dNormalC(auxf0[cc==0],0,1,c(a1,a2))/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf1[cc==1],0,1,c(a1,a2))-AcumNormalC(auxf[cc==1],0,1,c(a1,a2)))))
           return(-ver)
         }
         f13 <- function(nu){
           ver1 <-(sum(log(dNormalC(auxf0[cc==0],0,1,nu)/sqrt(sigma2)))+ sum(log(AcumNormalC(auxf1[cc==1],0,1,nu)-AcumNormalC(auxf[cc==1],0,1,nu))))
           return(ver1)
         }
         
         
         
         if(cens=="1"){
           Art <- optim(c(0.2,0.2),fnc1, method=c("BFGS"),control=list(maxit=20000))$par
           nu3 <- min(round(exp(Art[1])/(1+exp(Art[1]))+0.05,1),0.9)
           nu4 <- min(round(exp(Art[2])/(1+exp(Art[2]))+0.05,1),0.9)
           nu <- c(nu3,nu4) 
           logver=f11(nu)
         }
         
         if(cens=="2"){
           Art <- optim(c(0.2,0.2),fnc2, method=c("BFGS"),control=list(maxit=20000))$par
           nu3 <- min(round(exp(Art[1])/(1+exp(Art[1]))+0.05,1),0.9)
           nu4 <- min(round(exp(Art[2])/(1+exp(Art[2]))+0.05,1),0.9)
           nu <- c(nu3,nu4)
           logver=f12(nu)
         }
         
         if(cens=="3"){
           Art <- optim(c(0.2,0.2),fnc3, method=c("BFGS"),control=list(maxit=20000))$par
           nu3 <- min(round(exp(Art[1])/(1+exp(Art[1]))+0.05,1),0.9)
           nu4 <- min(round(exp(Art[2])/(1+exp(Art[2]))+0.05,1),0.9)
           nu <- c(nu3,nu4)
           logver=f13(nu)
         }         
         
        tetas <- matrix(c(as.vector(betas),sigma2,nu),ncol=1)      
       criterio <- sqrt(t(tetas-param)%*%(tetas-param))
       betas <- t(betas)
       param <- tetas       
       }
    AIC <- (-2)*logver + 2*(p+3)   ## p (Colunas de X) + Sigma2 + nu + Gamma
    BIC <- (-2)*logver + (p+3)*log(n)
    EDC <- (-2)*logver + (p+3)*0.2*sqrt(n)    
    SE1 <- IM.Cens(cc,x,y,LS, BETA=betas,SIGMA2=sigma2,NU=nu,delta=NULL,cens,type)
  obj.out <- list(betas = betas, sigma2 = sigma2, nu = nu, logver = logver, count=count, res=res, AIC=AIC, BIC=BIC, EDC=EDC, SE=SE1)
   
   class(obj.out) <- type
   return(obj.out) 
  }
  if (type == "Normal")
  {

       u0 <- u1 <- u2 <- c()   
 		param <- matrix(c(betas,sigma2),ncol=1)
     criterio <- 1      
   while(criterio > error || count<=iter.max)
     {
         count <- count + 1
            mu <- x%*%(betas)
        # Passo E

     NCensEUY <- NCensurEsperUY(y,mu,sigma2,nu=NULL,0,type=type)    
           u0 <- NCensEUY$EUY0
        	 u1 <- NCensEUY$EUY1
           u2 <- NCensEUY$EUY2
     
      CensEUY <- CensEsperUY1(mu[cc==1],sigma2=sigma2,nu=0,delta=0,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type=type, cens=cens)
     u0[cc==1]<- CensEUY$EUY0
     u1[cc==1]<- CensEUY$EUY1
     u2[cc==1]<- CensEUY$EUY2
            S1 <- matrix(0, n)
            S2 <- matrix(0, n)
            S3 <- matrix(0, n)
 	
         suma1 <- t(t(x)%*%u1)
         
         u0.matriz <- Diagonal(n,as.numeric(sqrt(u0)))
         xnovo <- as.matrix(u0.matriz%*%x)
         suma2 <- t(xnovo)%*%xnovo
         
       betas <- t(solve(suma2)%*%t(suma1))
   		sigma2 <- sum(u2-2*u1*mu+mu^2*u0)/n

    	 auxf0 <- (y-x%*%t(betas))/sqrt(sigma2)

 	 auxf <- (Lim1-x%*%t(betas))/sqrt(sigma2)

       auxf1 <- (Lim2-x%*%t(betas))/sqrt(sigma2)	

	if(cens=="3"){		
      logver <- sum(log(dnorm(auxf0[cc==0])/sqrt(sigma2)))+ sum(log(pnorm(auxf1[cc==1])-pnorm(auxf[cc==1])))  
	}
	if(cens=="1"){
	logver <- sum(log(dnorm(auxf0[cc==0])/sqrt(sigma2)))+ sum(log(pnorm(auxf0[cc==1]))) 
	}
	if(cens=="2"){
	logver <- sum(log(dnorm(auxf0[cc==0])/sqrt(sigma2)))+ sum(log(1-pnorm(auxf0[cc==1]))) 
	}
			
       tetas <- matrix(c(as.vector(betas),sigma2),ncol=1)      
    criterio <- sqrt(t(tetas-param)%*%(tetas-param))
       betas <- t(betas)
       param <- tetas       
       }
    AIC <- (-2)*logver + 2*(p+1)   ## p (Colunas de X) + Sigma2 
    BIC <- (-2)*logver + (p+1)*log(n)
    EDC <- (-2)*logver + (p+1)*0.2*sqrt(n)
    
    SE1 <- IM.Cens(cc, x,y,LS,BETA=betas,SIGMA2=sigma2,NU=nu,delta=NULL,cens,type)
    obj.out <- list(betas = betas, sigma2 = sigma2, logver = logver,count=count, res=res, AIC=AIC, BIC=BIC, EDC=EDC, SE=SE1)
   class(obj.out) <- type
   return(obj.out) 
  }
 }

#### Envelope baseado no Residuo Martingal
################################################################################

EnvelopeRMT<-function(cc,x,y,LS,nu=3,delta=NULL,beta=NULL,cens="1",type = "T")
{
  res <- residNI(y,x,cc,LS,nu=nu,delta=delta,cens,type)
  if(type=="Normal")
  {
  chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab=expression("r_MTi"),xlab="Standard normal quantile",main="N-CR",line="quartile")
  }
  if(type=="T")
  {
  chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab=expression("r_MTi"),xlab="Standard t quantile",main="T-CR",line="quartile")
  }
  if(type=="Slash")
  {
  chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab=expression("r_MTi"),xlab="Standard slash quantile",main="SL-CR",line="quartile")
  }
  if(type=="NormalC")
  {
  chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab=expression("r_MTi"),xlab="Standard contaminated normal quantile",main="CN-CR",line="quartile")
  }
  if(type=="PearsonVII")
  {
  chart.QQPlot(res$resmt,distribution='norm',envelope=0.95,ylab=expression("r_MTi"),xlab="Standard Paearson VII quantile",main="Pearson VII-CR",line="quartile")
  }
}

##### Calculando os Erros Padrões dos Betas
################################################################################

IM.Cens <- function(cc,x,y,LS,BETA,SIGMA2,NU,delta=0,cens="1",type="T")
{
     	x <- as.matrix(x)
      y <- as.vector(y)
  betas <- BETA
 sigma2 <- SIGMA2
     nu <- NU
    	#x <- cbind(1,x)
   	  p <- ncol(x)
    	n <- nrow(x)
    mu <- x%*%(betas)
 VarUY1 <- VarUY0 <- CovUY01 <- c()
    VAR <- VarU<- k <- c()
   prim <- sec <- c()
  EUY10 <- EUY11 <- EUY21 <- c()
     IM <- Aux <- IC1 <- A<- A1 <- matrix(0,nrow=p,ncol=p) 
    Var <- matrix(0,nrow=n,ncol=1) 
   Aux1 <- matrix(0,nrow=p,ncol=p)
     se <- matrix(0,1,p)   
   
        if(type=="Normal")
        {
         k <- 1  
        }
        if(type=="T")
        {
        #k <- (nu+1)/(nu+3)
         k <- MonteCarlo(25000,0,1,nu,delta=NULL,type=type)
        }
        if(type=="PearsonVII")
        {
         k <- MonteCarlo(25000,0,1,nu,delta,type=type)
        }
        if(type=="Slash")
        {
         k <- MonteCarlo(25000,0,1,nu,delta=0,type=type)
        }
        if(type=="NormalC")
        {
         k <- MonteCarlo(25000,0,1,nu,delta=0,type=type)
        }

        if (cens=="1")
        {
        Lim1 <- rep(-Inf,n)
        Lim2 <- y
        }
       if (cens=="2")
        {
       Lim1 <- y
       Lim2 <- rep(Inf,n)
        }
       if (cens=="3")
        {
       Lim1 <- y
       Lim2 <- LS
        }
        Aux <- matrix(0,nrow=p,ncol=p)
       Aux1 <- t(x)%*%x
        IC1 <- Aux1
 Var[cc==1] <- CensVar(mu[cc==1],sigma2,nu,delta,Lim1=Lim1[cc==1],Lim2=Lim2[cc==1],type=type,cens=cens)
  			for(i in 1:n)
        {
        Aux <- (Var[i]*x[i,])%*%t(x[i,])
        IM <- IM + Aux 
        }
   IM <- (IM/sigma2)/sigma2
	  IC1 <- (k*IC1)*(1/sigma2)
        A <- IC1  - IM
       A1 <- solve(A)
       se <- sqrt(diag(A1))
return(se)
}

        