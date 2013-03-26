\name{CensReg.SMN}
\alias{CensReg.SMN}
\title{ Fit univariate right or left censored linear regression model under the Scale Mixtures of Normal distributions }
\description{
Return EM algorithm output for right or left censored regression model (under Normal, Student-t, Pearson VII, Slash or Contaminated Normal distribution), built the corresponding envelope graph and calculate criterias for model selection such as AIC, BIC and EDC
}
\usage{CensReg.SMN(cc,x,y,LS,nu,delta,beta,cens,type,criteria,show.envelope,sd,error,iter.max)}
\arguments{
\item{cc}{ Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored }
\item{x}{ Matrix or vector of covariables }
\item{y}{ Vector of response variable in case of right/left censoring or vector of lower limits of censoring intervals }
\item{LS}{ Vector of upper limits of censoring intervals. Must not be provided in case of left or right censoring }
\item{nu}{ Parameter of the scale variable of the SMN family. Must not be provided in case of Normal distribution. Must be a bidimensional vector in case of Contaminated Normal distribution }
\item{delta}{ Second parameter of Pearson VII. Must not be provided in case of Normal, Student-t or Slash distribution }
\item{beta}{ Initial value for the regression coefficients. If "NULL", MQ estimate is adopted }
\item{cens}{ "1" for left censoring, "2" for right censoring, "3" for interval censoring }
\item{type}{ Distribution to be used: "Normal" for Normal model, "T" for Student-t model, "PearsonVII" for Pearson VII model, "Slash" for slash model and "NormalC" for contaminated Normal model }
\item{criteria}{ TRUE or FALSE. Indicates if log-likelihood, AIC, BIC and EDC should be computed }
\item{show.envelope}{ TRUE or FALSE. Indicates if envelope graph should be built for the fitted model (based on transformed Martingale residuals) }
\item{sd}{ TRUE or FALSE. Indicates if standard deviation of beta should be computed (based on the observed information matrix) }
\item{error}{ the convergence maximum error }
\item{iter.max}{ maximum iterations for the EM algorithm }
}
\details{ For the contaminated Normal distribution, each component of the bidimensional vector "nu" must lie on (0,1). }
\value{ 
\item{beta}{ EM estimatives for the regression coefficients }
\item{sigma2}{ EM estimative for the variance }
\item{logver}{ Value of the log-likelihood under the fitted model }
\item{count}{ Number of interations until convergence }
\item{AIC}{ AIC criteria for model selection. Only returned if criteria=TRUE }
\item{BIC}{ BIC criteria for model selection. Only returned if criteria=TRUE }
\item{EDC}{ EDC criteria for model selection. Only returned if criteria=TRUE }
\item{sd}{ Standard deviation of the estimators of beta (based on the observed information matrix). Only returned if sd=TRUE. }
}
\references{
 Aldo M. Garay, Victor H. Lachos, Heleno Bolfarine, Celso R. Cabral. "Censored Regression Models with Scale Mixture of Normal Distributions".
 }
\author{ Aldo M. Garay \email{aldowmg@ime.usp.br}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br} and Monique Bettio Massuia \email{moniquemassuia@gmail.com} }
\examples{ ##see examples in \code{\link{wage.rates}} }
\keyword{ censored linear regression }
