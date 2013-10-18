\name{wage.rates}
\docType{data}
\alias{wage.rates}
\title{ Wage rates of 753 women}
\description{ Wage rates of 753 married white women with left censoring 
}
\usage{data(wage.rates)}
\format{ 
A data frame with 753 observed wage rates of married white women in 1975. Some wage rates are set equal to zero, this means that those wives did not work in 1975, therefore, these observations are considered left censored at zero. 
}

\encoding{latin1}

\source{ 
Mroz, T.A. 1987. "The sesitivity of an empirical model of married women´s hours of work to economic and statistical assumptions". Full data set is available at 
}

\examples{
\dontrun{
##Load the data
data(wage.rates)

##Set the response y and covariate x
y <- wage.rates$wage
x <- cbind(wage.rates$age,wage.rates$educ,wage.rates$kidslt6,wage.rates$kidsge6)
cc<- c(rep(0,428),rep(1,325))

##Fits a left censored Student-t model to the data

T <- CensReg.SMN(cc,x,y,LS=NULL,nu=3,delta=NULL,cens="1",type="T",
     show.envelope="FALSE",error=0.0001,iter.max=300)

##Fits a left censored Normal model to the data

Normal <- CensReg.SMN(cc,x,y,LS=NULL,nu=NULL,delta=NULL,cens="1",type="Normal",
          show.envelope="FALSE",error=0.0001,iter.max=300)}

}
\keyword{datasets}



