## DESCRIPTION
createdata.f creates the pseudo dataset that to be used in the M-step.
Fbar.f, EIu.f,Eterm.f are used internally for createdata.f.
EM.f implements the method we proposed for variable selection of interval-censored data.

## ARGUMENTS
EM.f(indata, lam0, beta0, lasso.lam, ncov, npieces, cutpoints, penalty.function, penalty.factor = NULL, nopenalty.index = NULL, thresh = 10^-6, maxit = 200)

indata:               	a data frame containing variables with names ¡°id¡±, ¡°timeL¡±, ¡°timeR¡±, ¡°x1¡±, ¡­,  to denote the individuals¡¯ id, left endpoint, right endpoint and covaraites.
lam0:                 	initial values of the baseline hazards.
beta0:                	initial values of the coefficients.
lasso.lam:            	the tuning parameter.
ncov:                 	the number of covariates.
npieces:              	the number of pieces used in piecewise constant baseline hazard model.
cutpoints:            	a data frame containing the information of cut points.
penalty.function:     	the penalty function to be used, "lasso", "alasso", or "scad".
penalty.factor:       	separate penalty factors can be applied to each coefficient, this is a number that multiplies lambda to allow differential shrinkage.
nopenalty.index:      	the index of non-penalized covariates. 
thresh:		      	convergence threshold for EM algorithm; default is 1e-6.
maxit:			maximum number of passes over the data for all lambda values; default is 200.


## VALUE 
tol:    the value of tolerance at the last iteration.
iter:   the number of iterations until convergence.
beta:   the values of the coefficients.
lam:    the values of the baseline hazards.


## EXAMPLE
## the dataset have 100 observations and 10 covariates, with first 2 and last 2 covariates significant; the covariates are generated from multivariate normal distribution with mean 0 and covariance matrix with entry in the ith row and jth column equal to 0.5^{i-j}; the event time is generated from the Weibull distribution with shape equal to 1.25 and scale such that P(T < 1 | X = 0) = 0.95; the number of assessments are generated from a Poisson distribution with mean 10 and the assessment times are generated from the uniform distribution [0, 1]

## read data
indata <- read.table("example.data.txt", header = TRUE)

# decide the cutpoints 
fit0 <- survfit(Surv(timeL, timeR, status, type='interval')~1, data = indata)
surv.prob <- summary(fit0)$surv;
npieces <- 4
ncov <- 10
probs <- seq(from = 1-max(surv.prob), to = 1-min(surv.prob), length.out = npieces + 1)
probs <- probs[-c(1, npieces+1)]
cpoints <- quantile(fit0, probs = probs, conf.int = FALSE)
cutpoints <- data.frame(start=c(0,cpoints),
                        stop=c(cpoints,9999),
                        piece=1:(length(cpoints)+1))
  
## proposed method with LASSO penalty
## the penalty parameter is 0.01
fit.lasso <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov), 
     		lasso.lam = 0.01, ncov = ncov, npieces = npieces, 
     		cutpoints = cutpoints, penalty.function = "lasso")
fit.lasso
## proposed method with adaptive LASSO (ALASSO) penalty
fit.alasso <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov), 
                  lasso.lam = 0.01, ncov = ncov, npieces = npieces, 
     		  cutpoints = cutpoints, penalty.function = "alasso")
fit.alasso
## proposed method with SCAD penalty
fit.scad <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov), 
                  lasso.lam = 0.01, ncov = ncov, npieces = npieces,
                  cutpoints = cutpoints, penalty.function = "scad")
fit.scad
