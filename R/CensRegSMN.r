################################################################################
#####       Analysis of Normal/independent Censored Regression Model     #######
# cc: vector of censoring indicators: 0 if non-censored, 1 if censored
# x: matrix or vector of covariables
# y: vector of responde variable
# nu: initial value for the degree of freedon in case of Student-t model
# delta : initial value for the second parameter os Pearson VII
# beta : initial value for betas. If "NULL" MQ estimate is adopted
# cens : type of censoring, 1: Left, 2:Rigth
# type: Normal, T, Slash, NormalC or PearsonVII
# criteria: TRUE ou FALSE - should compute AIC, BIC and EDC criterias for model selection
# envelope: TRUE ou FALSE - should compute envelops
# sd: TRUE or FALSE - should compute standard deviation of betas´s estimators
# error
# iter.max = 100

CensReg.SMN <- function(cc, x=NULL,y=NULL ,LS=NULL, nu=NULL, delta=NULL,beta=NULL,cens="1",type="T",criteria="TRUE",show.envelope="FALSE", sd="TRUE", error=0.0001,iter.max=300)
{
	## Verify error at parameters specification
	if(ncol(as.matrix(y)) > 1) stop("Only univariate linear regression supported!")
  if(ncol(as.matrix(cc)) > 1) stop("Only univariate linear regression supported!")
	if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
	if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
	if( (type != "T") && (type != "Normal") && (type != "PearsonVII") && (type != "Slash") && (type != "NormalC")) stop("Distribution family not supported. Check documentation!")
	if( (type=="T") && (length(delta) > 0) ) stop("delta parameter must not be provided in case of Student-t distribution.")

	if(length(beta)!=(ncol(as.matrix(x))+1) & length(beta)!=0 ){stop("Wrong dimension for initial values of beta")}

	if( type == "T" | type == "Slash" )
  {
  		if(length(nu) > 1) stop("nu parameter must be a scalar")
  		if(length(nu) == 0) stop("nu parameter must be provided.")
  		if(nu <= 0) stop("nu parameter must be positive.")
  }
  if( type == "PearsonVII" )
  {
  		if(length(nu) > 1) stop("nu parameter must be a scalar")
  		if(length(nu) == 0) stop("nu parameter must be provided in case of Pearson VII distribution.")
  		if(nu <= 0) stop("nu parameter must be positive.")
  		if(length(delta) > 1) stop("delta parameter must be a scalar")
  		if(length(delta) == 0) stop("delta parameter must be provided in case of Pearson VII distribution.")
  		if(delta <= 0) stop("delta parameter must be positive.")
  }
  if(type == "Normal")
  {
  if( (length(nu)!= 0) | length(delta)!= 0 ) stop("nu and delta parameters must not be provided in case of Normal distribution.")
  }
  if(type == "NormalC")
  {
  if(length(nu) !=2) stop("nu must be a bidimensional vector in case of Contaminated Normal distribution")
  if(nu[1] <=0 || nu[1] >= 1) stop("nu[1] must lies in (0,1)")
  if(nu[2] <=0 || nu[2] >= 1) stop("nu[2] must lies in (0,1)")
  }
  if( (cens != "1") && (cens != "2") && (cens != "3"))  stop("Censored type not supported. 1 for left censoring, 2 for right censoring and 3 for intervalar censoring.")
  if(show.envelope=="TRUE")
  {
  EnvelopeRMT(cc,x,y,LS,nu,delta,beta,cens=cens,type=type)
  }
  out <- EM.Cens.Int(cc, x,y,LS,nu,delta,beta,cens,type,criteria, sd, error, iter.max)
  out
}


