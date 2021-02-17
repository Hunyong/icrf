##### metric functions #####

    # F.hat evaluation at gridline
    # S.hat.grid.i.training <- function(subj = 1, fold = 1, rsf = a, tau = 4, tick = 0.001) {
    #   grid.tmp <- seq(0, tau, by = tick)
    #   S <- 1 - cumsum(rsf$interval_prob_seq[[fold]][subj,])
    #   lin.interpolate.vec(grid.tmp, xvec = rsf$time_interest, yvec = S, return.y.only = TRUE) # F03-ICR-1base.R
    # }

    S.hat.grid.i <- function(S.hat.i, time_interest, grid) {
      lin.interpolate.vec(grid, xvec = time_interest, yvec = S.hat.i, return.y.only = TRUE) # F03-ICR-1base.R
    }

    # S.hat.grid.i.training(fold = 2)
    # tmp <- Muti_ERT_Predict(testset, Forest = a$Forest_seq[[1]], a$SurvMat_seq[[1]], a$time_interest)
    # S.hat.grid.i(S.hat.i = tmp$Surv_predict[1,], time_interest = tmp$time_interest, grid = seq(0, 4, by = 0.001))

    # F evaluation at gridline
    S.grid.i <- function(mu = 1, theta = 1, fun = "exp", grid) {
      if (fun == "exp") {
        return(1 - pexp(grid, rate = 1/mu))
      } else if (fun == "log-normal") {
        return(1 - pnorm(log(grid), mean = log(mu) - theta^2/2, sd = theta))
      } else if (fun == "gamma") {
        return(1 - pgamma(grid, shape = mu, scale = theta))
      } else if (fun == "exp.discrete") {
        s = pexp(grid, rate = 1/mu)
        s.d = pexp(floor(2 * grid)/2, rate = 1/mu)
        return(1 - (s + s.d)/2)
      }
    }

    cond.expt <- function(S.hat.grid, LR, grid) {
      L <- as.numeric(LR["L"])
      R <- as.numeric(LR["R"])
      if (is.infinite(R)) { R <- grid[length(grid)]} # truncated expectation: integration taken only up to tau.

      s.hat.cond <- S.hat.grid[grid >= L & grid <= R]

      n.cond <- length(s.hat.cond)
      S.LR <- range(s.hat.cond)
      if (n.cond <= 1 | S.LR[1] == S.LR[2]) {  # interval is null (essentially a point), or surv curv is flat in that interval.
        return((L + R)/2)
      }
      s.hat.cond <- (s.hat.cond - S.LR[1]) /(S.LR[2] - S.LR[1])
      return (L + sum(s.hat.cond) * (R - L) / n.cond)
    }
    # cond.expt(cox.hat[1,], LR = c(L = 1, R = 2), grid = grid)
    # c.index <- function(ET, ET.interval) {
    #   ord <- order(ET.interval)
    #   ET.interval <- ET.interval[ord]  #smallest to largest in terms of ET.interval
    #   ET <- ET[ord]
    #
    #   n.ET <- length(ET)
    #   result <- 0
    #   for (i in 1:(n.ET - 1)) {
    #     for (j in (i + 1):n.ET) {
    #       result <- result +
    #         if (ET.interval[j] - ET.interval[i]) { ## When T1 and T2 are not tied,
    #           (sign(ET[j] - ET[i]) + 1)/2          # 1 if ET1 > ET2, 0.5 if ET1 == ET2, 0 o/w
    #         } else {                               ## When T1 and T2 are tied,
    #           0.5 +  (ET[j] == ET[i])/2            # 1 if tie, 0.5 if no tie
    #         }
    #     }
    #   }
    #   2 * result / n.ET / (n.ET - 1)
    # }
    # # c.index(ET = cox.ET, ET.interval = cox.ET.cond)



    IBS <- function(S.hat.grid, LR, grid, weight = NULL) {
      # integrated brier score
      L <- as.numeric(LR["L"])
      R <- as.numeric(LR["R"])
      if (is.infinite(R)) { R <- grid[length(grid)]} # truncated expectation: integration taken only up to tau.

      indicator1 <- ifelse(grid <= L, 1, ifelse(grid >= R, 0, NA))
      indicator2 <- indicator1

      n <- length(grid)
      if (is.null(weight)) {weight = rep(1, n)}
      sum.weight1 = sum(weight * !is.na(indicator1), na.rm = TRUE)
      sum.weight2 = sum(weight, na.rm = TRUE)

      n.censor <- sum(is.na(indicator1))
      s.hat.cond <- S.hat.grid[is.na(indicator1)]

      if (length(s.hat.cond) == 0) {
        indicator2[is.na(indicator1)] <- 0.5
      } else {
        S.LR <- range(s.hat.cond)
        den <- S.LR[1] - S.LR[2]
        indicator2[is.na(indicator1)] <- (s.hat.cond - S.LR[1]) / den
# tmp.a <- sum(is.na(indicator1)); tmp.b <- length(s.hat.cond);
# if (tmp.a %% tmp.b) {print(s.hat.cond); print(indicator1); print(c(tmp.a, tmp.b))}
      }

      score1 <- sum((indicator1 - S.hat.grid)^2 * weight, na.rm = TRUE) / sum.weight1
      # note, if every grid points are censored, score1 = NaN. Such objects will be ignored later.
      score2 <- sum((indicator2 - S.hat.grid)^2 * weight, na.rm = TRUE) / sum.weight2

      return (c(score1 = score1, score2 = score2))
    }

    # deviations for a subject (integrated, sup absolute error and integrated squared error)
    dev.i <- function(S.hat.grid, S.grid, tick, LR, grid) {
      diff <- S.hat.grid - S.grid
      int.err <- sum(abs(diff))*tick
      sup.err <- max(abs(diff))
      # int.sq.err <- sum(diff^2)*tick
      brier <- IBS(S.hat.grid = S.hat.grid, LR = LR, grid = grid)
      return(c(int.err = int.err, sup.err = sup.err, IBS1 = brier["score1"], IBS2 = brier["score2"]))
    }

    # lin.interpolate(49.25, xvec =  48:52, yvec = (1:5)/3)
    # S()
    # dev.i(subj = 1, fold = 1, rsf = a, tau = 4, tick = .001, mu = 1, theta = 1, fun = "exp")

    IBS.mat <- function(S.hat.grid.mat, LR.mat, grid, weight = NULL) {
      n = dim(LR.mat)[1]
      apply(
        sapply(1:n, function(s) IBS(S.hat.grid.mat[s, ], LR.mat[s, ], grid = grid, weight = weight)),
        1, mean)
    }


##### base functions #####
    search2 <- function(s, xvec, include = TRUE) {
      # Finding the index of "the largest x smaller or equal to s"
      # Fn is the cdf table (x (should be ordered) and y)
      # index = which (s >= c(-Inf, xvec) & s <= c(xvec, Inf)) - 1
      index = if (include) {which (s >= c(-Inf, xvec)) - 1} else {which (s > c(-Inf, xvec)) - 1}
      max(index)
    }
    # search2(49.25, xvec =  48:52)
    lin.interpolate <- function(s, xvec, yvec = NULL, return.y.only = FALSE) {
      x.len = length(xvec)
      low.index <- search2(s, xvec)
      lower <- xvec[low.index]
      if (low.index == x.len) {
        upper <- lower
        proportion <- 0
      } else {
        upper <- xvec[low.index + 1]
        proportion <- (s - lower) / (upper - lower)
      }
      result <- c(low.index = low.index, proportion = proportion)
      if (!is.null(yvec)) {
        lower.y <- yvec[low.index]
        upper.y <- yvec[min(x.len, low.index + 1)]
        y.interpolate <- lower.y + proportion * (upper.y - lower.y)
        if (return.y.only) return(y.interpolate)
        result["y.interpolate"] <- y.interpolate
      }
      return(result)
    }
    lin.interpolate.vec <- Vectorize(lin.interpolate, vectorize.args = "s")
    # lin.interpolate(49.25, xvec =  48:52)
    S.grid.i <- function(mu = 1, theta = 1, fun = "exp", grid) {
      if (fun == "exp") {
        return(1 - pexp(grid, rate = 1/mu))
      } else if (fun == "log-normal") {
        return(1 - pnorm(log(grid), mean = log(mu) - theta^2/2, sd = theta))
      } else if (fun == "gamma") {
        return(1 - pgamma(grid, shape = mu, scale = theta))
      } else if (fun == "exp.discrete") {
        s = pexp(grid, rate = 1/mu)
        s.d = pexp(floor(2 * grid)/2, rate = 1/mu)
        return(1 - (s + s.d)/2)
      }
    }


##### data functions #####
    dat.PH <- function(n = 200, n.monitor = 3, tau = 4, rho = 0.9, b0 = 0.1, b1 = 0.1, discrete = FALSE) {
      P = 25     # dimension of X
      # n = 200    # training sample size
      # tau = 4
      # n.monitor = 3           # case-n.monitor interval censoring
      # rho = 0.9
      # b0 = 0.1

      Sigma = outer(1:P, 1:P, function(x, y) rho^{abs(x - y)})

      dataX = mvrnorm(n = n, mu = rep(0, P), Sigma = Sigma)
      colnames(dataX) <- paste0("X", 1:P)

      # mu = exp(b0 * sum_i=11^20 x_i - 0.5)
      mu = 
        exp(b0 * (as.matrix(dataX) %*% 
                  matrix(rep(c(0, 1, 0), times = c(10, 10, 5)), ncol = 1)) - b1)

      y = rexp(n, 1/mu)
      if (discrete) {
        y.discrete = ceiling(y * 2)/2
        y = ifelse(rbinom(n, 1, 0.5), y, y.discrete)  # with 50% chance discrete.
      }
      C = t(sapply(1:n, function(x) sort(rexp(n.monitor, 1/mean(mu))))) # matrix of monitoring times
      if (n.monitor == 1) C = t(C)
      C[C > tau] <- Inf   # if censoring time is greater than tau, essentially no censoring or Inf.
      # C = pmin(C, tau)  # this assumption is when monitoring time is too large, monitor early at tau and end study.

      interval.index = sapply(1:(n), function(s) search2(y[s], c(0, C[s, ], Inf), include = FALSE))
      L = cbind(0, C, Inf)[cbind(1:(n), interval.index)]
      R = cbind(0, C, Inf)[cbind(1:(n), interval.index + 1)]

      # censor = (y <= C)

      dataset = data.frame(cbind(dataX, L, R, mu = as.numeric(mu)))

      return(dataset)
    }


    ### 2. violation of PH 1
    dat.NPH1 <- function(n = 200, n.monitor = 3, tau = 6) {
      P = 10     # dimension of X
      # n = 200    # training sample size
      # tau = 6
      # n.monitor = 3           # case-n.monitor interval censoring

      dataX = as.data.frame(lapply(1:P, function(x) runif(n)))
      colnames(dataX) <- paste0("X", 1:P)

      # mu = sin(X1 *pi) + 2|X2 - .5| + X3^3
      mu = sin(dataX$X1 * pi) + 2 * abs(dataX$X2 - .5) + dataX$X3^3

      y = rexp(n, 1/mu)
      C = t(sapply(1:n, function(x) sort(runif(n.monitor)* tau))) # matrix of monitoring times
      if (n.monitor == 1) C = t(C)
      C[C > tau] <- Inf   # if censoring time is greater than tau, essentially no censoring or Inf.

      interval.index = sapply(1:(n), function(s) search2(y[s], c(0, C[s, ], Inf), include = FALSE))
      L = cbind(0, C, Inf)[cbind(1:(n), interval.index)]
      R = cbind(0, C, Inf)[cbind(1:(n), interval.index + 1)]

      # censor = (y <= C)

      dataset = data.frame(cbind(dataX, L, R, mu = as.numeric(mu)))

      return(dataset)
    }

    #dat.NPH1()

    ### 3. violation of PH 2
    dat.NPH2 <- function(n = 200, n.monitor = 3, tau = 10, rho = 0.75, theta = 2) {
      P = 25     # dimension of X
      # n = 200    # training sample size
      # tau = 10
      # n.monitor = 3           # case-n.monitor interval censoring
      # rho = 0.75

      Sigma = outer(1:P, 1:P, function(x, y) rho^{abs(x - y)})

      dataX = mvrnorm(n = n, mu = rep(0, P), Sigma = Sigma)
      colnames(dataX) <- paste0("X", 1:P)

      # mu = 0.5 + 0.3 * |x11+ ... + x15|
      mu = 0.5 + 0.3 * abs(apply(dataX[, 11:15], 1, sum))
      y = rgamma(n, mu, scale = theta)
      C = t(sapply(1:n, function(x) sort(runif(n.monitor)* tau * 1.5))) # matrix of monitoring times
      if (n.monitor == 1) C = t(C)
      C[C > tau] <- Inf   # if censoring time is greater than tau, essentially no censoring or Inf.

      interval.index = sapply(1:(n), function(s) search2(y[s], c(0, C[s, ], Inf), include = FALSE))
      L = cbind(0, C, Inf)[cbind(1:(n), interval.index)]
      R = cbind(0, C, Inf)[cbind(1:(n), interval.index + 1)]

      # censor = (y <= C)

      dataset = data.frame(cbind(dataX, L, R, mu = as.numeric(mu)))

      return(dataset)
    }

    #dat.NPH2()


    ### 4. conditional independence
    dat.CI <- function(n = 200, n.monitor = 3, tau = 4, rho = 0.2, theta = 1) {
      P = 25     # dimension of X
      # n = 200    # training sample size
      # tau = 4
      # n.monitor = 3           # case-n.monitor interval censoring
      # rho = 0.75

      Sigma = outer(1:P, 1:P, function(x, y) rho^{abs(x - y)})

      dataX = mvrnorm(n = n, mu = rep(0, P), Sigma = Sigma)
      colnames(dataX) <- paste0("X", 1:P)

      # mu = 0.3 * |x1+ ... + x5| + 0.3 * |x21+ ... + x25| +
      mu = 0.3 * abs(apply(dataX[, 1:5], 1, sum)) + 0.3 * abs(apply(dataX[, 21:25], 1, sum))
      y = exp(rnorm(n, log(mu) - theta^2/2, sd = theta))   # => log-normal with mean = mu
      # mean of log N(m, s) = exp(m + s^2/2).      => m = log(mu) - s^2/2
      C = t(sapply(1:n, function(x) exp(sort(rnorm(n.monitor, mean = mu[x] * .8, sd = 1))))) # matrix of monitoring times
      if (n.monitor == 1) C = t(C)
      C[C > tau] <- Inf   # if censoring time is greater than tau, essentially no censoring or Inf.

      interval.index = sapply(1:(n), function(s) search2(y[s], c(0, C[s, ], Inf), include = FALSE))
      L = cbind(0, C, Inf)[cbind(1:(n), interval.index)]
      R = cbind(0, C, Inf)[cbind(1:(n), interval.index + 1)]

      # censor = (y <= C)

      dataset = data.frame(cbind(dataX, L, R, mu = as.numeric(mu)))

      return(dataset)
    }
    dat.CI.B <- function(n = 200, n.monitor = 3, tau = 4, rho = 0.2, theta = 1) {
      P = 25     # dimension of X
      # n = 200    # training sample size
      # tau = 4
      # n.monitor = 3           # case-n.monitor interval censoring
      # rho = 0.75
      
      Sigma = outer(1:P, 1:P, function(x, y) rho^{abs(x - y)})
      
      dataX = mvrnorm(n = n, mu = rep(0, P), Sigma = Sigma)
      colnames(dataX) <- paste0("X", 1:P)
      
      # mu = 0.1 * |x1+ ... + x5| + 0.1 * |x21+ ... + x25| +
      mu = 0.3 * abs(apply(dataX[, 1:5], 1, sum)) + 0.3 * abs(apply(dataX[, 21:25], 1, sum))
      y = exp(rnorm(n, log(mu) - theta^2/2, sd = theta))   # => log-normal with mean = mu
      # mean of log N(m, s) = exp(m + s^2/2).      => m = log(mu) - s^2/2
      C = t(sapply(1:n, function(x) exp(sort(rnorm(n.monitor, mean = log(mu[x] * .8) - theta^2/2, sd = 1))))) # matrix of monitoring times
      if (n.monitor == 1) C = t(C)
      C[C > tau] <- Inf   # if censoring time is greater than tau, essentially no censoring or Inf.
      
      interval.index = sapply(1:(n), function(s) search2(y[s], c(0, C[s, ], Inf), include = FALSE))
      L = cbind(0, C, Inf)[cbind(1:(n), interval.index)]
      R = cbind(0, C, Inf)[cbind(1:(n), interval.index + 1)]
      
      # censor = (y <= C)
      
      dataset = data.frame(cbind(dataX, L, R, mu = as.numeric(mu)))
      
      return(dataset)
    }
    
    #dat.CI()

    ### 5. dependent
    dat.D <- function(n = 200, n.monitor = 3, tau = 2, rho = 0.75) {
      P = 10     # dimension of X
      # n = 200    # training sample size
      # tau = 10
      # n.monitor = 3           # case-n.monitor interval censoring
      # rho = 0.75

      Sigma = outer(1:P, 1:P, function(x, y) rho^{abs(x - y)})

      dataX = mvrnorm(n = n, mu = rep(0, P), Sigma = Sigma)
      colnames(dataX) <- paste0("X", 1:P)

      # mu = 2 expit(x1 + x2 + x3)
      mu = 2 * plogis(apply(dataX[, 1:3], 1, sum))
      y = rexp(n, 1/mu)
      # censoring mechanism: exponentials with mean = y

      C = t(sapply(1:n, function(x) exp(sort(rnorm(n.monitor, log(y[x]) - 0.5, sd = 1))))) # matrix of monitoring times
      if (n.monitor == 1) C = t(C)
      C[C > tau] <- Inf   # if censoring time is greater than tau, essentially no censoring or Inf.

      interval.index = sapply(1:(n), function(s) search2(y[s], c(0, C[s, ], Inf), include = FALSE))
      L = cbind(0, C, Inf)[cbind(1:(n), interval.index)]
      R = cbind(0, C, Inf)[cbind(1:(n), interval.index + 1)]

      # censor = (y <= C)

      dataset = data.frame(cbind(dataX, L, R, mu = as.numeric(mu)))

      return(dataset)
    }

    #dat.D()
    
    
##### Inducing censoring #####
cens.nlms = function(data, n.monitor = 1, tau, remove.original = FALSE, time.var = "time", delta.var = "delta") {
  n = dim(data)[1]  
  mu = with(data, 1000 + 100 * (10 - age/10 + hhnum))
  C = t(sapply(1:n, function(x) sort(rnorm(n.monitor, mean = mu[x], sd = 300)))) # matrix of monitoring times
  if (n.monitor == 1) C = t(C)
print(summary(C))
  C[C > tau] <- Inf   # if censoring time is greater than tau, essentially no censoring or Inf.
  C[C < 0] <- 0       # if censoring time is less than 0, force it to zero.
  interval.index = sapply(1:(n), function(s) search2(data$time[s], c(0, C[s, ], Inf), include = FALSE))
  data$L = cbind(0, C, Inf)[cbind(1:(n), interval.index)]
  data$R = cbind(0, C, Inf)[cbind(1:(n), interval.index + 1)]
  if (remove.original) {
    data[[time.var]] = data[[delta.var]] = NULL
  }
  data
}
# cens.nlms(nlms.complete, tau = 2192)
