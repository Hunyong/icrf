# ==========================================================
# Authors: Ying Wu and Richard Cook
# EM Algorithm for Interval Censored Data
# Subroutines - create pseudo datasets
#
# Date: Nov 04, 2014
# ==========================================================


# ----------------------------------------------------------
# Create Pseudo Data Set
# ----------------------------------------------------------

createdata.f <- function(indata, lam, beta, ncov, cutpoints, timeL = "timeL", timeR = "timeR", X = "x") {
  outdata <- lapply(1:nrow(indata), function(ith, indata, lam, beta, ncov, cutpoints) {
    out <- Eterm.f(Li=indata[ith, timeL],
                   Ri=indata[ith, timeR],
                   zi=as.vector(unlist(indata[ith, paste(X,1:ncov,sep="")])),
                   lam=lam, beta=beta,
                   ncov=ncov,
                   cutpoints=cutpoints)
    nlen <- length(out$piece)
    
    out$id <- rep(ith, nlen)
    out <- data.frame(out[,c("id","piece","EIk","Ewk","logEwk")])
    
    outcov <- apply(indata[ith,paste(X,1:ncov,sep="")], 2, rep, times = nlen)
    outcov <- matrix(outcov, nrow=nlen, ncol=ncov)
    row.names(outcov) <- 1:nlen
    outcov <- data.frame(outcov)
    
    outdata <- cbind(out, outcov)
    return(outdata)
  }, indata=indata, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
  outdata <- do.call("rbind", outdata)
  outdata <- data.frame(outdata)
  dimnames(outdata)[[2]] <- c("id","piece","EIk","Ewk","logEwk",paste("x",1:ncov,sep=""))
  outdata <- outdata[order(outdata$id),]
  return(outdata)
}     

# ----------------------------------------------------------
# Evaluate Fbar(t[i] | Z[i]; theta)
# ----------------------------------------------------------

Fbar.f <- function(ti, nzi, zi, lam, beta, ncov, cutpoints) {
  cstart <- cutpoints$start
  cstop  <- cutpoints$stop
  npiece <- length(cstart)
  nsubj  <- length(ti)

  beta <- matrix(beta, nrow=ncov, ncol=1)
  zi   <- matrix(zi, nrow=nzi, ncol=ncov)
  zi.times.beta <- as.vector( zi %*% beta )

  Ht <- rep(0,nsubj)
  for (k in 1:npiece) {
    wk <- rep(0,nsubj)
    wk <- ifelse( (ti <= cstop[k]) & (ti >= cstart[k]), ti - cstart[k], wk)
    wk <- ifelse( (ti > cstart[k]) & (ti > cstop[k]), rep(cstop[k] - cstart[k],nsubj), wk)

    Ht <- Ht + ( (wk*lam[k])*exp( zi.times.beta ) )
  }

  Fbar <- exp( (-1)*Ht )
  return(Fbar)
}


# ----------------------------------------------------------
# E-Step: 
# Evaluate E(I[k](u[i]) | D[i]) and E(S[k](u[i]) | D[i]))
# ----------------------------------------------------------

EIu.f <- function(ti, zi, lam, beta, ncov, cutpoints, Fterm) {
  Fbar <- Fbar.f(ti=ti, nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
  term <- Fbar - Fterm
  term <- ifelse(is.na(term),  0, term)
  term <- ifelse(is.nan(term), 0, term)
  return( term )
}
  

Eterm.f <- function(Li, Ri, zi, lam, beta, ncov, cutpoints) {
  npiece <- max(cutpoints$piece)
      
  outi <- cutpoints
  outi$EIk <- rep(0, npiece)
  outi$Ewk <- rep(0, npiece)
  
  if(Ri == 9999) {
    for (k in 1:npiece) {  
      EIk <- 0 
      Ewk <- min(c(Li,outi$stop[k])) - min( c(Li, outi$start[k]) )
      outi$EIk[k] <- EIk
      outi$Ewk[k] <- Ewk
    } 
  }
  else{
    FbarLi <- Fbar.f(ti=Li, nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
    FbarRi <- Fbar.f(ti=Ri, nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
    
    for (k in 1:npiece) {   
      EIk <- 0
      Ewk <- 0
      if (Ri < outi$start[k]) {
        EIk <- 0
        Ewk <- 0
      }
      else if (outi$stop[k] < Li) {
        EIk <- 0
        Ewk <- outi$stop[k] - outi$start[k]
      }
      else {
        FbarLik <- Fbar.f(ti=max(c(Li, outi$start[k])), nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
        FbarRik <- Fbar.f(ti=min(c(outi$stop[k], Ri)),  nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
        
        EIk <- (FbarLik - FbarRik)/(FbarLi - FbarRi)
        
        term1 <- max( c(Li - outi$start[k], 0) )
        
        lower.val <- max( c(Li, outi$start[k]) )
        upper.val <- min( c(outi$stop[k], Ri) )
        upper.val <- ifelse(upper.val == 9999, Inf, upper.val)
        int.term <- integrate(EIu.f, lower=lower.val, upper=upper.val,
                              zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints, Fterm=FbarRi)
        term2 <- int.term$value/(FbarLi - FbarRi)
        Ewk <- term1 + term2
      }
      
      outi$EIk[k] <- EIk
      outi$Ewk[k] <- Ewk
    } 
  }

  outi$logEwk <- log(outi$Ewk)
  outi$logEwk <- ifelse(outi$Ewk == 0, NA, outi$logEwk)
      
  return( outi[,c("piece","EIk","Ewk","logEwk")] )
}
  


