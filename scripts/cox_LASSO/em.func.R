# ==========================================================
# Authors: Ying Wu and Richard Cook
# EM Algorithm for Interval Censored Data
# 
#
# Date: Nov 04, 2014
# ==========================================================

require(MASS)
require(survival)
require(glmnet)


EM.f <- function(indata, lam0, beta0, lasso.lam, ncov, npieces, cutpoints, penalty.function, 
                 penalty.factor = NULL, nopenalty.index = NULL, thresh = 1e-6, maxit = 200,
                 xlabel = "x", Llabel = "timeL", Rlabel = "timeR") {

  cur.lam  <- lam0
  cur.beta <- beta0
  
  if(is.null(penalty.factor)){
    penalty.factor <- penalty.factor.f(beta0=cur.beta, lasso.lam=lasso.lam, 
                                       ncov=ncov, npieces=npieces, 
                                       penalty.function=penalty.function,
                                       nopenalty.index = nopenalty.index)
  }
  
  iter <- 0
  tol  <- 9999
  while ( tol > thresh ) {
    if(penalty.function=="alasso" | penalty.function == "scad"){
      penalty.factor <- penalty.factor.f(beta0=cur.beta, lasso.lam=lasso.lam, 
                                         ncov=ncov, npieces=npieces, 
                                         penalty.function=penalty.function,
                                         nopenalty.index = nopenalty.index)
    }
    iter <- iter + 1
    
    pre.lam  <- cur.lam
    pre.beta <- cur.beta
    
    pseudodata <- createdata.f(indata=indata, lam=pre.lam, beta=pre.beta, ncov=ncov, cutpoints=cutpoints,
                               timeL = Llabel, timeR = Rlabel, X = xlabel)
    pseudodata <- pseudodata[!is.na(pseudodata$EIk),]
    pseudodata <- pseudodata[!is.na(pseudodata$logEwk),]
    
    if ( npieces == 1 ) {
      y <- pseudodata$EIk
      x <- as.matrix( pseudodata[,c(paste("x",1:ncov,sep=""))] )
      adjterm <- pseudodata$logEwk
      
      fit <- glmnet(x=x, y=y, family="poisson", offset=adjterm,
                    penalty.factor=penalty.factor,
                    lambda=lasso.lam, alpha=1)
      cf <- coef(fit)
      cur.lam  <- exp( cf[1] )
      cur.beta <- cf[-1]
      
      y <- NULL; x <- NULL; adjterm <- NULL; fit <- NULL; cf <- NULL
    }
    else {
      pieces <- matrix(0, nrow=nrow(pseudodata), ncol=(npieces-1))
      for (k in 2:npieces) {
        pieces[,k-1] <- ifelse(pseudodata$piece == k, 1, 0)
      }   
      y <- pseudodata$EIk
      x <- as.matrix( cbind(pieces, as.matrix(pseudodata[,c(paste("x",1:ncov,sep=""))])) )
      adjterm <- pseudodata$logEwk
      
      fit <- glmnet(x=x, y=y, family="poisson", offset=adjterm,
                    penalty.factor=penalty.factor,
                    lambda=lasso.lam, alpha=1)
      
      cf <- as.vector(coef(fit))
      cur.beta <- cf[-c(1:npieces)]
      
      cf <- cf[c(1:npieces)]
      cur.lam <- exp( c(0, cf[-1]) + cf[1] )
      
      pieces <- NULL; y <- NULL; x <- NULL; adjterm <- NULL; fit <- NULL; cf <- NULL
    }    
    cur.lam  <- ifelse(is.na(cur.lam), 0, cur.lam)
    cur.beta <- ifelse(is.na(cur.beta), 0, cur.beta)
    
    dif.lam  <- abs( (cur.lam - pre.lam) / pre.lam )
    dif.beta <- abs( (cur.beta - pre.beta) / pre.beta )
    dif.beta <- ifelse(pre.beta == 0, abs(cur.beta - pre.beta), dif.beta)
    
    tol <- max( c(dif.lam, dif.beta) )  
    if ( iter > maxit ) { break }
  }
  
  out <- NULL
  out$tol  <- tol
  out$iter <- iter
  out$beta <- as.vector(cur.beta)
  out$lam  <- as.vector(cur.lam)
  return(out)
  
}


penalty.factor.f <- function(beta0 = NULL, lasso.lam = NULL, ncov, npieces, penalty.function, nopenalty.index = NULL){
  if(penalty.function == "lasso"){
    pf <- rep(1, ncov)
    pf[nopenalty.index] <- 0
    penalty.factor <- c(rep(0, npieces-1), pf)
  }
  if(penalty.function == "alasso"){
    beta0.na <- beta0
    beta0.na[beta0 == 0] <- 10^(-5)
    pf <- 1/abs(beta0.na)
    pf[nopenalty.index] <- 0
    penalty.factor <- c(rep(0, npieces-1), pf)
  }
  scad.pen.deriv <- function(beta0, lasso.lam){
    abs.beta0 <- abs(beta0)
    indicator <- ifelse(abs.beta0 > lasso.lam, 0, 1)
    deriv <- lasso.lam*(indicator + (1-indicator)*ifelse(3.7*lasso.lam - abs.beta0 >= 0, 3.7*lasso.lam - abs.beta0, 0)/2.7/lasso.lam)
    return(deriv)
  }
  if(penalty.function == "scad"){
    deriv <- scad.pen.deriv(beta0=beta0, lasso.lam=lasso.lam)
    deriv[nopenalty.index] <- 0
    penalty.factor <- c(rep(0, npieces-1),deriv)
  }
  return(penalty.factor)
}


