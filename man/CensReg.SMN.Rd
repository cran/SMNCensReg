\name{CensReg.SMN}
\alias{CensReg.SMN}
\title{ Fit univariate right, left or interval censored linear regression model
under the scale mixtures of Normal (SMN) distributions }
\description{
Return EM algorithm output for right, left or interval censored regression model under SMN distributions
(Normal, Student-t, Pearson VII, Slash or Contaminated Normal), built the corresponding
envelope graph and compute some criteria for model selection, such as, AIC, BIC and EDC.
}
\usage{
CensReg.SMN(cc,x,y,LS=NULL,nu=3,delta=NULL,cens="1",type="T",
show.envelope="FALSE", error=0.0001,iter.max=300)
}
\arguments{
\item{cc}{ Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored }
\item{x}{ Matrix or vector of covariates }
\item{y}{ Vector of responses in case of right/left censoring. Vector of lower limits if censoring is intervalar }
\item{LS}{ Vector of upper limits if interval censoring. Must not be provided in case of left/right censoring.}
\item{nu}{ Initial value of the parameter of the scale variable of the SMN family. Must not
be provided in case of Normal distribution. Must be a bidimensional vector in
case of contaminated Normal distribution (NormalC) }
\item{delta}{ Initial value of the second parameter of Pearson VII. Must not be provided in case of Normal, Student-t or slash distribution }
\item{cens}{ "1" for left censoring, "2" for right censoring, "3" for interval censoring }
\item{type}{ Distribution to be used in fitting: "Normal" for Normal model, "T" for Student-t model, "PearsonVII" for Pearson VII model, "Slash" for slash model and "NormalC" for contaminated Normal model }
\item{show.envelope}{ TRUE or FALSE. Indicates if envelope graph should be built for the fitted model (based on transformed Martingale residuals). Default is FALSE }
\item{error}{ The convergence maximum error }
\item{iter.max}{ The maximum number of iterations of the EM algorithm. Default=300 }
}
\details{ For the contaminated Normal distribution, each component of the bidimensional vector "nu" must
 lie on (0,1) For the Pearson VII distribution, delta=1. The parameters beta and sigma2 are initialized
with the minimum square estimators of the regression x vs y.}
\value{
\item{beta}{ EM estimates for the regression coefficients. }
\item{sigma2}{ EM estimates for the scale parameters. }
\item{logver}{ Returned the value of the log-likelihood under the fitted model. }
\item{count}{ Number of interations until convergence. }
\item{AIC}{ AIC criteria for model selection.  }
\item{BIC}{ BIC criteria for model selection.  }
\item{EDC}{ EDC criteria for model selection.  }
\item{SE}{ Standard error of the estimates of the regression coefficients. }
}
\references{
 Aldo M. Garay, Victor H. Lachos, Heleno Bolfarine, Celso R. Cabral. "Censored Regression Models with Scale Mixture of Normal Distributions".
 Tech. Rep. 13--13. IMECC--UNICAMP.
 }
\author{ Aldo M. Garay \email{aldowmg@ime.usp.br}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br} and Monique Bettio Massuia \email{moniquemassuia@gmail.com} }
\examples{ ##see examples in \code{\link{wage.rates}} }
\keyword{ censored linear regression }
