###### Distribuições Acumuladas das Normal Independ (NI)
################################################################################
   
PearsonVII <- function(y,mu,sigma2,nu,delta)
     {
     Acum <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
     for (i in 1:length(y))
     { 
      z[i] <- (y[i]-mu)/sqrt(sigma2a)
     Acum[i] <- pt(z[i],df=nu)
     }
    return(Acum)
     }


P <- function(y,mu,sigma2,nu,delta)
     {
  A <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
        n <- length(y)
        i <- 0
     while (i<n)
     { 
      i <- i +1      
   z[i] <- (y[i]-mu)/sqrt(sigma2a)
   A[i] <- pt(z[i],df=nu)
     }
  return(A)
}

AcumSlash <- function(y,mu,sigma2,nu)
  {
  Acum <- z <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
   { 
  z[i] <- (y[i]-mu)/sqrt(sigma2)
    f1 <- function(u) nu*u^(nu-1)*pnorm(z[i]*sqrt(u))
Acum[i]<- integrate(f1,0,1)$value	 	
   }
return(Acum)
     }

AcumNormalC <- function(y,mu,sigma2,nu)
   {
  Acum <- vector(mode = "numeric", length = length(y))
 for (i in 1:length(y))
   { 
  eta  <- nu[1]    
  gama <- nu[2]
Acum[i] <- eta*pnorm(y[i],mu,sqrt(sigma2/gama)) + (1-eta)*pnorm(y[i],mu,sqrt(sigma2))
   }
return(Acum)
   }


###### Densidades das Normais Independentes
################################################################################

dPearsonVII<- function(y,mu,sigma2,nu,delta)
    {
     f <-  z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2*(delta/nu)
     for (i in 1:length(y))
     { 
     z[i] <- (y[i]-mu)/sqrt(sigma2a)
     f[i] <- dt(z[i],df=nu)/sqrt(sigma2a)
     }
     return(f)
     }

dSlash <- function(y,mu,sigma2,nu)
     {
     resp <- z <- vector(mode = "numeric", length = length(y))
     for (i in 1:length(y)) 
     {
     z[i] <- (y[i]-mu)/sqrt(sigma2)
     f1 <- function(u) nu*u^(nu-0.5)*dnorm(z[i]*sqrt(u))/sqrt(sigma2)
     resp[i] <- integrate(f1,0,1)$value	 	
     }
     return(resp)
     }

dNormalC <- function(y,mu,sigma2,nu)
    {
     Acum <- vector(mode = "numeric", length = length(y))
     for (i in 1:length(y))
     { 
    eta  <- nu[1]    
    gama <- nu[2]
   Acum[i]  <- eta*dnorm(y[i],mu,sqrt(sigma2/gama)) + (1-eta)*dnorm(y[i],mu,sqrt(sigma2))
     }
   return(Acum)
    }

###### E(UY) Resultados obtidos na Tese de Osorio para as NI
################################################################################

NCensurEsperUY <-  function(y,mu=2,sigma2=1,nu=0,delta=0,type="Normal")
{
	   EUY0 <- EUY1 <- EUY2 <- c()
     d <- (y-mu)^2/sigma2
     n <- length(y)  	   
    if(type=="T")
    {
    EUY0 <- (nu+1)/(nu+d) 
    EUY1 <- y*(nu+1)/(nu+d) 
    EUY2 <- (y^2)*(nu+1)/(nu+d) 
    }
  	if(type=="Normal")
    {
    EUY0 <- rep(1,n) 
    EUY1 <- y
    EUY2 <- (y^2)
    }
  	if(type=="PearsonVII")
    {
       d <- (y-mu)^2/sigma2
    EUY0 <- (nu+1)/(delta+d) 
    EUY1 <- y*(nu+1)/(delta+d) 
    EUY2 <- (y^2)*(nu+1)/(delta+d) 
    }
    if(type=="Slash")
    {
    Num  <- GamaInc(nu+1.5,0.5*d)*(0.5*d)^(-nu-1.5) 
    Den  <- GamaInc(nu+0.5,0.5*d)*(0.5*d)^(-nu-0.5) 
    EUY0 <- Num/Den 
    EUY1 <- y*Num/Den
    EUY2 <- (y^2)*Num/Den
    }
    if(type=="NormalC")
    {
    Num  <- 1-nu[1]+nu[1]*(nu[2])^(1.5)*exp(0.5*d*(1-nu[2]))
    Den  <- 1-nu[1]+nu[1]*(nu[2])^(0.5)*exp(0.5*d*(1-nu[2]))
    EUY0 <- Num/Den 
    EUY1 <- y*Num/Den
    EUY2 <- (y^2)*Num/Den
    }
return(list(EUY0=EUY0,EUY1=EUY1,EUY2=EUY2))         
}

###### Função Auxiliar 
################################################################################

GamaInc <- function(a,x)
     {
   res <- vector(mode = "numeric", length = length(x))
     f <-function(t) exp(-t)*t^(a-1)
     for  (i in 1:length(x))
     {
     res[i] <- integrate(f,0,x[i])$value
     }
     return(res)
     }

E_phi <- function(r,a, nu=3,delta=0,type=type, cens=cens)
{
        n <- length(a)
        b <- rep(Inf,n)
        b1<- rep(-Inf,n)
                
        if(setequal(a,b)== TRUE | setequal(a,b1)== TRUE)
        {
        resp <- rep(0,n)
        }
        else
        {
        if(type=="Normal")
        {
#         if (a==Inf | a == -Inf)
#         {
#         resp <- 0
#         }
#         else
#         {
         resp <- dnorm(a)
#         }
        }
       if(type=="T")
        {
#         if (a==Inf | a == -Inf)
#         {
#         resp <- 0
#         }
#         else
#        {
         Aux0 <- gamma(0.5*(nu+2*r))
         Aux1 <- gamma(nu/2)*sqrt(2*pi) 
         Aux2 <- Aux0/Aux1
         Aux3 <- (0.5*nu)^(nu/2)
         Aux4 <- (0.5*(a^2+nu))^(-0.5*(nu+2*r))
         resp <- Aux2*Aux3*Aux4 
#         }
        }
        if(type=="PearsonVII")
        {
#         if (a==Inf | a== -Inf)
#         {
#         resp <- 0
#         }
#         else
#         {
         Aux0 <- gamma(0.5*(nu+2*r))
         Aux1 <- gamma(nu/2)*sqrt(2*pi) 
         Aux2 <- Aux0/Aux1
         Aux3 <- (0.5*delta)^(nu/2)
         Aux4 <- (0.5*(a^2+delta))^(-0.5*(nu+2*r))
         resp <- Aux2*Aux3*Aux4 
#         }
        }
        if(type=="Slash")
        {
#         if (a ==Inf | a == -Inf)
#         {
#         resp <- 0
#         }
#         else
#         {
         Aux0 <- nu/sqrt(2*pi) 
         Aux1 <- (0.5*a^2)^(-(nu+r))
         Aux2 <- GamaInc(nu+r,0.5*a^2)
         resp <- Aux0*Aux1*Aux2 
#         }
        } 
      if(type=="NormalC")
        {
 #        if (a ==Inf | a == -Inf)
 #        {
 #        resp <- 0
 #        }
 #        else
 #        {
         Aux0 <- nu[1]*nu[2]^(r)*dnorm(a*sqrt(nu[2]))
         Aux1 <- (1-nu[1])*dnorm(a)
         resp <- Aux0 + Aux1
 #         }
        }
   }      
return(resp)
}

E_Phi <- function(r,a, nu=3,delta=0,type=type)
{
        n <- length(a)
        if(type=="Normal")
        {
        resp <- pnorm(a)
        }        
        if(type=="T")
        {
        Aux0 <- gamma(0.5*(nu+(2*r)))
        Aux1 <- gamma(nu/2) 
        Aux2 <- Aux0/Aux1
        Aux3 <- (0.5*nu)^(-r)
        Aux4 <- PearsonVII(a,0,1,nu+(2*r),nu)
        resp <- Aux2*Aux3*Aux4
	      }
        if(type=="PearsonVII")
        {
        Aux0 <- gamma(0.5*(nu+(2*r)))
        Aux1 <- gamma(nu/2) 
        Aux2 <- Aux0/Aux1
        Aux3 <- (0.5*delta)^(-r)
        Aux4 <- PearsonVII(a,0,1,nu+(2*r),delta)
        resp <- Aux2*Aux3*Aux4  
        }
        if(type=="Slash")
        {
        Aux0 <- nu/(nu+r)
        Aux1 <- AcumSlash(a,0,1,nu+r)
        resp <- Aux0*Aux1
        }        
        if(type=="NormalC")
        {
        Aux0 <- nu[2]^(r)*AcumNormalC(a,0,1,nu)
        Aux1 <- (1-nu[2]^(r))*(1-nu[1])*pnorm(a)
        resp <- Aux0 + Aux1
        } 
return(resp)
}

CensEsperUY1 <- function(mu=0,sigma2=1,nu=3,delta=0,Lim1=1,Lim2=4,type=type, cens=cens)
{
        Lim11 <- (Lim1-mu)/sqrt(sigma2)        
        Lim21 <- (Lim2-mu)/sqrt(sigma2)
            n <- length(Lim11) 
        if(type=="Normal")
        {
            EU <-  1
        FNIb   <-  pnorm(Lim21)
        FNIa   <-  pnorm(Lim11)
        }
        if(type=="T")
        {
            EU <-  1
        FNIb   <-  pt(Lim21,nu)
        FNIa   <-  pt(Lim11,nu)
        }
        if(type=="PearsonVII")
        {
        FNIb   <-  PearsonVII(Lim21,0,1,nu,delta)
        FNIa   <-  PearsonVII(Lim11,0,1,nu,delta)
        }
        if(type=="Slash")
        {
        FNIb   <- AcumSlash(Lim21,0,1,nu)
        FNIa   <- AcumSlash(Lim11,0,1,nu)
        }
        if(type=="NormalC")
        {
           EU <-  (nu[1]*nu[2]) + (1-nu[1])
        FNIb  <- AcumNormalC(Lim21,0,1,nu)
        FNIa  <- AcumNormalC(Lim11,0,1,nu)
        }
        
         if (cens=="1")
         {
         Aux11 <- rep(0,n)
         }else
         {
         Aux11 <- Lim11
         }
         if (cens=="2")
         {
         Aux22 <- rep(0,n)
         }else
         {
         Aux22 <- Lim21
         }
#         if (Lim11==-Inf)
#         {
#         Aux11 <- 0
#         }else
#         {
#         Aux11 <- Lim11
#         }
#         if (Lim21==Inf)
#         {
#         Aux22 <- 0
#         }else
#         {
#         Aux22 <- Lim21
#         }

   K <- 1/(FNIb-FNIa)
EUX0 <- K*(E_Phi(1,Lim21, nu,delta,type)- E_Phi(1,Lim11, nu,delta,type))    
EUX1 <- K*(E_phi(0.5,Lim11,nu,delta,type,cens)- E_phi(0.5,Lim21,nu,delta,type,cens))
EUX2 <- K*(E_Phi(0,Lim21, nu,delta,type)- E_Phi(0,Lim11, nu,delta,type) + Aux11*E_phi(0.5,Lim11,nu,delta,type,cens) - Aux22*E_phi(0.5,Lim21,nu,delta,type,cens))           # Neste c aso r =2   
EUX20 <- K*(E_Phi(2,Lim21, nu,delta,type)- E_Phi(2,Lim11, nu,delta,type))
EUX21 <- K*(E_phi(1.5,Lim11,nu,delta,type,cens)- E_phi(1.5,Lim21,nu,delta,type,cens))
#EUX22 <- K*(E_Phi(1,Lim21, nu,delta,type)- E_Phi(1,Lim11, nu,delta,type) + Aux11*E_phi(1.5,Lim11,nu,delta,type,cens) - Aux22*E_phi(1.5,Lim21,nu,delta,type,cens))           # Neste c aso r =2   

EUY0  <- EUX0          
EUY1  <- mu*EUX0 + sqrt(sigma2)*EUX1
EUY2  <- EUX0*mu^2 + 2*mu*sqrt(sigma2)*EUX1 + sigma2*EUX2
EUY20 <- EUX20
EUY21 <- mu*EUX20 + sqrt(sigma2)*EUX21
#EUY22 <- EUX20*mu^2 + 2*mu*sqrt(sigma2)*EUX21 + sigma2*EUX22
#return(list(EUY0=EUY0,EUY1=EUY1,EUY2=EUY2, EUY20=EUY20, EUY21=EUY21, EUY22=EUY22))
return(list(EUY0=EUY0,EUY1=EUY1,EUY2=EUY2, EUY20=EUY20, EUY21=EUY21))
}

## para o cálculo de Im sob censura (caso t de Student)
CensVar <- function(mu=0,sigma2=1,nu=3,delta=0,Lim1=1,Lim2=4,type="T", cens){
	  Lim11 <- (Lim1-mu)/sqrt(sigma2)        
    Lim21 <- (Lim2-mu)/sqrt(sigma2)
        n <- length(mu)
	if(type=="Normal"){
        FNIb   <-  pnorm(Lim21)
        FNIa   <-  pnorm(Lim11)
      }
      if(type=="T"){
        FNIb   <-  pt(Lim21,nu)
        FNIa   <-  pt(Lim11,nu)
      }
      if(type=="PearsonVII"){
        FNIb   <-  PearsonVII(Lim21,0,1,nu,delta)
        FNIa   <-  PearsonVII(Lim11,0,1,nu,delta)
      }
      if(type=="Slash"){
        FNIb   <- AcumSlash(Lim21,0,1,nu)
        FNIa   <- AcumSlash(Lim11,0,1,nu)
      }
      if(type=="NormalC"){
        FNIb  <- AcumNormalC(Lim21,0,1,nu)
        FNIa  <- AcumNormalC(Lim11,0,1,nu)
      }
         if (cens=="1")
         {
         Aux11 <- rep(0,n)
         }else
         {
         Aux11 <- Lim11
         }
         if (cens=="2")
         {
         Aux22 <- rep(0,n)
         }else
         {
         Aux22 <- Lim21
         }

k <-(1/(FNIb-FNIa))
Dif.X0 <- (E_Phi(2,Lim21,nu,delta,type)- E_Phi(2,Lim11,nu,delta,type) )
Dif.r1 <- (E_Phi(1,Lim21,nu,delta,type)- E_Phi(1,Lim11,nu,delta,type) )

VarUX0 <- k*(Dif.X0-(k*(Dif.r1^2)))   
VarUX1 <- k*(Dif.r1 + Aux11*E_phi(0.5,Lim11, nu,delta,type,cens) - Aux22*E_phi(0.5,Lim21, nu,delta,type,cens) - k*(Dif.r1^2))
Esp <- CensEsperUY1(mu=0,sigma2=1,nu,delta,Lim1=Lim11,Lim2=Lim21,type,cens)
CovU.UX <- Esp$EUY21 - (Esp$EUY0*Esp$EUY1)
   
Termo1 <- mu*VarUX0*mu + sigma2*VarUX1 + 2*mu*sqrt(sigma2)*CovU.UX

Termo2 <- mu*mu*VarUX0

Esp2 <- CensEsperUY1(mu,sigma2,nu,delta,Lim1=Lim1,Lim2=Lim2,type,cens)
Termo3 <- mu*( Esp2$EUY21 - (Esp2$EUY1*Esp2$EUY0)  )

resp <- Termo1 + Termo2 - 2*Termo3
 
return(resp)
}

### Gerar dados
##########################################

rmix <- function(n, p1, rF1, arg1=NULL, arg2=NULL) {
    #Função para gerar misturas de DUAS populações
    #n: qtd de amostras a ser gerada
    #p1: proporção de mistura da população 1 (1-p1 é prop da 2)
    #rF1: deve ser do tipo string - indica a função a ser chamada para gerar dela (por exemplo rnorm, rbeta, rgamma, ...)
    #arg1, arg2: argumentos a serem passados para a função escolhida em rF1, da pop 1 e 2, respectivamente
    x1 <- vector(mode = "numeric", length = n)
    for (i in 1:n){
      u <- runif(1)
      if (u < p1) x1[i] <- do.call("rF1", c(list(1), arg1))
      if (u > p1) x1[i] <- do.call("rF1", c(list(1), arg2))
    }
    return(x1)
    }


geraPearsonVII <- function(n,mu,sigma2,nu,delta)
{
    y <- vector(mode = "numeric", length =n)
    for (i in 1:n)
    {
    y[i] <- mu[i] + (rgamma(1, 0.5*nu, 0.5*delta))^(-1/2)*rnorm(1,0,sqrt(sigma2))
    }
    return(y)
}

geraS <- function(n,mu, sigma2,nu)
{
     y <- vector(mode = "numeric", length =n)
    u2 <- rbeta(n,nu,1) 
    for (i in 1:n)
    {
    y[i] <- mu[i] + (u2[i])^(-1/2)*rnorm(1,0, sqrt(sigma2))
    }
    return(y)
}

geraNC <- function(n,mu,sigma2,nu)
{
     y <- vector(mode = "numeric", length =n)
    p1<- nu[1]
    gama<-nu[2]
    u<-runif(1)
    for(i in 1:n)
    {
    if(u<p1) y[i]<- rnorm(1,mu[i],sqrt(sigma2/gama))
    if(u>p1) y[i]<- rnorm(1,mu[i],sqrt(sigma2))
    }
    return(y)
}

### Gerar dados censurados NI
##########################################

geraDadosR<-function(y,cc,x,betas,sigma2,nu,delta=NULL,cens="1",type="Normal")
{
	x <- cbind(1,x)
  x <- as.matrix(x)
  p <- ncol(x)
	n <- nrow(x)
	mu=x%*%betas
	ind<-matrix(0,n,1)

  if(type=="Normal")
  {
    resp=rnorm(n,mu,sqrt(sigma2))
    if(cens=="1")
    {
    ind[cc==1]<- resp[cc==1]<y[cc==1]+0
    resp[ind==1]<-y[ind==1]
    }
    if(cens=="2")
    {
    ind[cc==1]<- resp[cc==1]>y[cc==1]+0
    resp[ind==1]<-y[ind==1]
    }
  }
  if(type=="PearsonVII")
  {
     resp<- geraPearsonVII(n,mu,sigma2,nu,delta)
     if(cens=="1")
     {
     ind[cc==1]<- resp[cc==1]<y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
     if(cens=="2")
     {
     ind[cc==1]<- resp[cc==1]>y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
  }
    if(type=="T")
  {
     resp<- mu+sqrt(sigma2)*rt(n,df=nu)
     if(cens=="1")
     {
     ind[cc==1]<- resp[cc==1]<y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
     if(cens=="2")
     {
     ind[cc==1]<- resp[cc==1]>y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
  }
  if(type=="Slash")
  {
     resp<- geraS(n,mu,sigma2,nu)
     if(cens=="1")
     {
     ind[cc==1]<- resp[cc==1]<y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
     if(cens=="2")
     {
     ind[cc==1]<- resp[cc==1]>y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
  }
  if(type=="NormalC")
  {
     resp<- geraNC(n,mu,sigma2,nu)
     if(cens=="1")
     {
     ind[cc==1]<- resp[cc==1]<y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
     if(cens=="2")
     {
     ind[cc==1]<- resp[cc==1]>y[cc==1]+0
     resp[ind==1]<-y[ind==1]
     }
  }
return(list(y=resp,cc=ind))
}
#N1  <- geraDadosR(y,cc,x,Normal$betas,Normal$sigma2,0,0,cens="1",type="Normal")
#T1  <- geraDadosR(y,cc,x,T$betas,T$sigma2,T$nu,0,cens="1",type="T")
#S1  <- geraDadosR(y,cc,x,Slash$betas,Slash$sigma2,Slash$nu,0,cens="1",type="Slash")
#NC1 <- geraDadosR(y,cc,x,NormalC$betas,NormalC$sigma2,c(0.5,0.5),0,cens="1",type="NormalC")


### Programa para Gerar dados das NI
######################################

 gera.ni<- function(n,mu,sigma2,nu,delta,type=type)
    {
      if(type=="Normal")
      {
      y<- mu +rnorm(n,0,sqrt(sigma2))
      }
      if(type=="T")
      {
      y<- mu + sqrt(sigma2)*rt(n,df=nu)
      }
      if(type=="PearsonVII")
      {
      y<- geraPearsonVII(n,mu,sigma2,nu,delta)
      }
      if(type=="Slash")
      {
      y<- geraS(n,mu,sigma2,nu)
      }
      if(type=="NormalC")
      {
      y<-geraNC(n,mu,sigma2,nu)
      }
     return(y=y)
    }
    
    
 ### Amostras MonteCarlo  para o IM
######################################

MonteCarlo <- function(m,mu,sigma2,nu, delta, type=type)
{
  y <- u <- w <- d <- vector(mode = "numeric", length = m)
emp <- Num <- Den <- c()
  for (i in 1:m)
  {
    Num <- 0          
    Den <- 0
    if(type=="T")
    {
 y[i] <- mu + sqrt(sigma2)*rt(1,df=nu)
 u[i] <- y[i]^2
 w[i] <- u[i]*((nu+1)/(nu+u[i]))^2
 d[i] <- w[i]
    }
    if(type=="PearsonVII")
    {
 y[i] <- geraPearsonVII(1,mu,sigma2,nu,delta)
 u[i] <- y[i]^2
 w[i] <- (1/4)*u[i]*((nu+1)/(delta+u[i]))^2
 d[i] <- 4*w[i]
    }
    if(type=="Slash")
    {
 y[i] <- geraS(1,mu,sigma2,nu)
 u[i] <- y[i]^2 
  Num <- GamaInc(nu+1.5,0.5*u[i])*(0.5*u[i])^(-nu-1.5) 
  Den <- GamaInc(nu+0.5,0.5*u[i])*(0.5*u[i])^(-nu-0.5) 
 w[i] <- (1/4)*u[i]*(Num/Den)^2
 d[i] <- 4*w[i]
    }
    if(type=="NormalC")
    {
 y[i] <-geraNC(1,mu,sigma2,nu)
 u[i] <- y[i]^2 
 Num  <- 1-nu[1]+nu[1]*(nu[2])^(1.5)*exp(0.5*u[i]*(1-nu[2]))
 Den  <- 1-nu[1]+nu[1]*(nu[2])^(0.5)*exp(0.5*u[i]*(1-nu[2]))
 w[i] <- (1/4)*u[i]*(Num/Den)^2
 d[i] <- 4*w[i]
    }
  }
    emp <- (1/m)*sum(d)
 return(emp=emp)
}
   
#    mu <- 0
#sigma2 <- 1
    nu <- 6
 delta <- NULL
   type=  "T"
 (nu+1)/(4*(nu+3))

#MonteCarlo(m=10000,mu,sigma2,nu, delta, type=type)
#MonteCarlo(m=50000,mu,sigma2,nu, delta, type=type)  
#MonteCarlo(m=100000,mu,sigma2,nu, delta, type=type)

################################################################################
# ANALISE DE RESIDUOS MARTINGALES 
################################################################################
### Este programa calcula o Martingale type residual (m) e o Mantingale type
### residual transformed (mt)
################################################################################

residNI <- function(y,x,cc,LS,nu=4,delta=0,cens="1",type="T")
{
resm <- resmt <- S <- ypad <- vector(mode = "numeric", length = length(y))
   n <- length(y)
  em <- CensReg.SMN(cc, x, y ,LS, nu=nu, delta=NULL,beta= NULL, cens=cens,type=type,show.envelope="FALSE", SE="FALSE",criteria = "FALSE", error=0.0001,iter.max=300)
betas1<- em$betas
sigma21<-em$sigma2
    nu1<- em$nu
delta1<-em$delta
   x1<- cbind(1,x)
   mu<- x1%*%betas1 
  ypad <- (y-mu)/sqrt(sigma21) 
 if (type=="T")
 {
 S <- 1 - pt(ypad,nu1)
 }
 if (type=="Normal")
 {
 S <- 1- pnorm(ypad)
 }
if (type=="PearsonVII")
 {
 S <- 1- PearsonVII(ypad,0,1,nu1,delta1)
 }
 if (type=="Slash")
 {
 S <- 1- AcumSlash(ypad,0,1,nu1)
 }
 if (type=="NormalC")
 {
  S <- 1-AcumNormalC(ypad,0,1,nu1)
 }
 for (i in 1:n)
  {
  resm[i] <-  1-cc[i]+log(S[i])
 resmt[i] <-  sqrt(-2*((1-cc[i])*log(1-cc[i]-resm[i])+resm[i]))*sign(resm[i])
  }
  return(list(resm=resm,resmt=resmt))
}