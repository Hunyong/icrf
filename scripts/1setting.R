### 1.setting.R
### An auxiliary code for the simulations (1A, 1B, and 1C). 
### The setting() function uses global assignments deploying the necessary simulation parameters to the environment,
###  and setting up the folders, file names, etc.
### rf(), fu(), cox() are the generic functions made for an easy control of the arguments.

setting <- function(scenario, sim, n.monitor, ntree, pilot = 0, ticksize = 0.01, date = NULL, 
                    n = NULL, tau = 5, b1 = 0.1, 
                    simClass = "main") {
    
    path_output <<- paste0("output/", if (is.null(date)) Sys.Date() else date,"/")
    if (!dir.exists("output/")) dir.create("output/")
    if (!dir.exists(path_output)) {
      dir.create(path_output)
      message(paste0(path_output, " was not available and is now created"))
    }
    if (is.null(n)) ntrain <<- 200 else ntrain <<- n
    tau <<- tau
    
    # file names
    if (simClass == "main") {
      fn_output <<- paste0(path_output, "sim_scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-rep_", sim, ".rds")
      fn_eval   <<- paste0(path_output, "eval_scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-rep_", sim, ".rds")
    } else if (simClass == "size") {
      fn_output <<- paste0(path_output, "sizeSim_size_", ntrain, "-scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-rep_", sim, ".rds")
      fn_eval   <<- paste0(path_output, "sizeEval_size_", ntrain, "-scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-rep_", sim, ".rds")
    } else if (simClass == "honesty") {
      fn_output <<- paste0(path_output, "NonHonestySim_scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-rep_", sim, ".rds")
      fn_eval   <<- paste0(path_output, "NonHonestyEval_scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-rep_", sim, ".rds")
    } else if (simClass == "time") {
      fn_output <<- paste0(path_output, "timeSim_size_", ntrain, "-scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-nFold_", nfold, "-rep_", sim, ".rds")
      fn_eval   <<- paste0(path_output, "timeEval_size_", ntrain, "-scenario_", scenario, "-n.m_", n.monitor,
                           "-nT_", ntree, "-nFold_", nfold, "-rep_", sim, ".rds")
    }
    
    if (pilot == 1) {
      fn_output <<- gsub("(s|S)im", "sim_pilot", fn_output)
      fn_eval   <<- gsub("(e|E)val", "eval_pilot", fn_eval)
    }
    
    if (scenario == 1) {
      titles <<- "scenario 1: proportional hazard"
      P <<- 25                                        # data dimensions
      rho <<- 0.9; b0 <<- 0.1; theta <<- NA           # scenario1 specific params
      dat.gen <<- dat.PH
      dat.args <<- list(n.monitor = n.monitor, tau = tau, rho = rho, b0 = b0, b1 = b1, discrete = FALSE)
      S.fun <<- "exp"
    } else if (scenario == 2) {
      titles <<- "scenario 2: non-proportional hazard 1"
      P <<- 10                                        # data dimensions
      rho <<- NA; b0 <<- NA; theta <<- NA             # scenario1 specific params
      dat.gen <<- dat.NPH1
      dat.args <<- list(n.monitor = n.monitor, tau = tau)
      S.fun <<- "exp"
    } else if (scenario == 3) {
      titles <<- "scenario 3: non-proportional hazard 2"
      P <<- 25                                        # data dimensions
      rho <<- 0.75; b0 <<- NA; theta <<- 2            # scenario1 specific params
      dat.gen <<- dat.NPH2
      dat.args <<- list(n.monitor = n.monitor, tau = tau, rho = rho, theta = theta)
      S.fun <<- "gamma"
    } else if (scenario == 4) {
      titles <<- "scenario 4: NPH, conditionally independent censoring"
      P <<- 25                                        # data dimensions
      rho <<- 0.2; b0 <<- NA; theta <<- 1             # scenario1 specific params
      #dat.gen <<- dat.CI
      dat.gen <<- dat.CI.B
      dat.args <<- list(n.monitor = n.monitor, tau = tau, rho = rho, theta = theta)
      S.fun <<- "log-normal"
    } else if (scenario == 5) {
      titles <<- "scenario 5: NPH, dependent censoring"
      P <<- 10                                        # data dimensions
      rho <<- 0.75; b0 <<- NA; theta <<- NA           # scenario1 specific params
      dat.gen <<- dat.D
      dat.args <<- list(n.monitor = n.monitor, tau = tau, rho = rho)
      S.fun <<- "exp"
    } else if (scenario == 6) {
      titles <<- "scenario 6: proportional hazard, non-smooth"
      P <<- 25                                        # data dimensions
      rho <<- 0.9; b0 <<- 0.1; theta <<- NA           # scenario1 specific params
      if (is.null(b1)) b1 <- 0
      dat.gen <<- dat.PH
      dat.args <<- list(n.monitor = n.monitor, tau = tau, rho = rho, b0 = b0, b1 = b1, discrete = TRUE)
      S.fun <<- "exp.discrete"
    } else {
      stop("scenario was not defined.")
    }
    seed.no <<- sim * 10000
    set.seed(seed.no + 0)
    train  <<- do.call(dat.gen, c(n = ntrain, dat.args))
    set.seed(seed.no + 9)
    Test   <<- do.call(dat.gen, c(n = ntest, dat.args))
    Grid   <<- seq(0, tau, by = ticksize)
    Grid   <<- c(Grid, Inf)
    s.test <<- t(sapply (Test[, "mu"], function(s) S.grid.i(mu = s, theta = theta, fun = S.fun, grid = Grid)))
    s.train <<- t(sapply (train[, "mu"], function(s) S.grid.i(mu = s, theta = theta, fun = S.fun, grid = Grid)))

    mtry   <<- ceiling(sqrt(P))           # tree parameters - data dependent
    train.tmp  <<- do.call(dat.gen, c(n = 2, dat.args))
    form1  <<- as.formula(paste("Surv(L, R, type = 'interval2')~",
                              paste(names(train.tmp)[1:P], collapse = " + ")))

    time.bgn <<- Sys.time()
      
    # print(paste0("file.exist? ", file.exists(fn2)))
    print(paste0("mtry = ", mtry))
    message("##### [", titles, "] - simulation replicate [", sim, " out of ", n.sim, "] #####")
    message("titles, ntrain, P, tau, rho, b0, theta, dat.gen, dat.args, S.fun, seed.no, train, Test, Grid, s.test, mtry, train.tmp, form1, fn_output, fn_err, time.bgn, rf.base.args, and others have been globally assigned!!!")
    message(paste0("Starting time is ", time.bgn))
}

rf <- function(...) {
  icrf:::icrf.default(x = train[, 1:P], L = train$L, R = train$R, timeSmooth = Grid, tau = tau,
                      xtest = Test[,1:P], ytest = s.test,
                      keep.forest = T, proximity = T, ntree = ntree,
                      nodesize = nmin, nfold = nfold, returnBest = TRUE, imse.monitor = 1,
                      ...)
}

Fu <- function(RF = T, smoothing = T, split.rule = "PetoLogrank") {
  print(paste0(ifelse(RF, "FuRF", "FuTree"), ifelse(smoothing, " with", " without"), " smoothing" ))
  if (RF) {
    nSamp = ceiling(.632 * ntrain)
    repl = T
    mTry = mtry
    nTree = ntree
    n.min = nmin
  } else {
    nSamp = ntrain
    repl = F
    mTry = P
    nTree = 1
    n.min = nmin.t
  }

  if (smoothing) {
    bw = NULL
  } else {
    bw = 0
  }
  icrf:::icrf.default(x=train[, 1:P], L = train$L, R = train$R, timeSmooth = Grid, tau = tau,
                       xtest = Test[,1:P], ytest = s.test, initialSmoothing = FALSE,
                       keep.forest = T, proximity = T,
                       nodesize = n.min, bandwidth = bw,
                       sampsize = nSamp, replace = repl, mtry = mTry, ntree = nTree, nfold = 1,
                       split.rule = split.rule, quasihonesty = T, ERT = F)
}

cox <- function(formula, smooth = FALSE) {
  require(icenReg)
  require(dplyr)
  time.bgn = Sys.time()
  mod1 <- icenReg::ic_sp(formula, model = "ph", data = train)
  if (formula[[3]] == 1) nullmod = TRUE else nullmod = FALSE
  hat.fn <- function(n, dat, model, bandwidth, nullModel = FALSE) {
    sapply(1:n, function(i){
      sp_curves <- if (!nullModel) icenReg::getSCurves(model, dat[i, ]) else icenReg::getSCurves(model) 
      sp_int <- t(sp_curves$Tbull_ints[-1, ])
      sp_curve <- sp_curves$S_curves[[1]]
      nn <- length(sp_curve)
      sp_curve[is.na(sp_curve)] <- 0  # force NaN to be zero (S(near end) = 0)
      npmle <- list(
        intmap = sp_int,                           # First row is redundant.
        pf = sp_curve[-nn] - sp_curve[-1]   # Likewise, first row is redundant
      )

      s.hat <- 1 - icrf:::isdSm(LR = matrix(c(L = train$L, R = train$R), ncol=2), 
                                grid.smooth = Grid, btt = c(bandwidth, 0, tau),
                                npmle = npmle)
    }) %>% t
  }
  
  train.hat    <- hat.fn(n = ntrain, dat = train, model = mod1, nullModel = nullmod,
                         bandwidth = if (smooth) NA else 0)
  test.hat     <- hat.fn(n = ntest, dat = Test, model = mod1, nullModel = nullmod,
                         bandwidth = if (smooth) NA else 0)
  
# tmp.train.cox <<- train.hat
# tmp.test.cox <<- test.hat
# tmp.list <<- list(Grid = Grid, tau = tau, L = train$L, R = train$R,
#                   s.test = s.test)

  imseerr <- measure(train.hat, Grid, tau = tau, method = "imse",
                    L = train$L, R = train$R)
  interr <- measure(test.hat, Grid, tau = tau, method = "int",
                    surv.true = s.test)
  print(matrix(interr, ncol = 2))
  time.end = Sys.time()
  list(cox = mod1,
       predictedNO.Sm = train.hat,
       imse.oob = matrix(c(imse.type1 = NaN, imse.type2 = NaN), ncol = 2),
       imse.NO = matrix(imseerr, ncol = 2),
       test = list(predicted = test.hat, testerror = matrix(interr, ncol = 2)),
       runtime = data.frame(begin = time.bgn, end = time.end, elapsed = time.end - time.bgn))
}


coxLASSO <- function(formula, smooth = FALSE) {
  
  fit0 <- survfit(Surv(L, R, rep(3, ntrain), type='interval')~1, data = train)
  surv.prob <- summary(fit0)$surv;
  npieces <- round(sqrt(tau) * 10)
  ncov <- P
  probs <- seq(from = 1-max(surv.prob), to = 1-min(surv.prob), length.out = npieces + 1)
  probs <- probs[-c(1, npieces+1)]
  cpoints <- quantile(fit0, probs = probs, conf.int = FALSE)
  cutpoints <- data.frame(start=c(0,cpoints),
                          stop=c(cpoints, tau * 1.5),
                          piece=1:(length(cpoints)+1))
  
  # times <- sort(unique(c(train$L, train$R)))
  # times <- c(times[times < tau], tau)
  # npieces <- length(times)
  # cutpoints <- data.frame(start= c(0, times),
  #                         stop = c(times, 9999),
  #                         piece= 1:(npieces + 1))
  
  mod1 <- EM.f(indata = train, lam0 = rep(1, npieces), beta0 = rep(0, P), 
               lasso.lam = 0.01, ncov = P, npieces = npieces, 
               cutpoints = cutpoints, penalty.function = "alasso",
               xlabel = "X", Llabel = "L", Rlabel = "R")
tmp.mod1 <<- mod1
## consider "smoothed" Cox Lasso as well!!!!
  
  hat.fn <- function(n, dat, model, bandwidth) {
      sapply(1:n, function(i){
        haz.ratio <- exp(sum(model$beta * dat[i, 1:P]))
        base.haz.cum <- cumsum(model$lam)
        surv <- exp(- haz.ratio * base.haz.cum)
        pf <- c(1, surv[-npieces]) - surv
        nn <- npieces
        npmle <- list(
          intmap = t(cutpoints[, c("start", "stop")]),                           # First row is redundant.
          pf = pf   # Likewise, first row is redundant
        )
        s.hat <- 1 - isdSm(LR = NULL, grid.smooth = Grid, btt = c(bandwidth, 0, tau),
                           npmle = npmle)
      }) %>% t
  }
  if (smooth == "both" | smooth == 0) {
    train.hat1    <- hat.fn(n = ntrain, dat = train, model = mod1, bandwidth = 0)
    test.hat1     <- hat.fn(n = ntest, dat = Test, model = mod1, bandwidth = 0)
    imseerr <- measure(train.hat1, Grid, t0 = 0, tau = tau, method = "imse",
                      L = train$L, R = train$R)
    interr <- measure(test.hat1, Grid, t0 = 0, tau = tau, method = "int",
                      surv.true = s.test)
    result1 <- 
      list(cox = mod1,
         predictedNO.Sm = train.hat1,
         imse.oob = matrix(c(imse.type1 = NaN, imse.type2 = NaN), ncol = 2),
         imse.NO = matrix(imseerr, ncol = 2),
         test = list(predicted = test.hat1, testerror = matrix(interr, ncol = 2)))
  }
  if (smooth == "both" | smooth == 1) {
    train.hat2    <- hat.fn(n = ntrain, dat = train, model = mod1, bandwidth = tau * ntrain^(-1/5)/2)
    test.hat2     <- hat.fn(n = ntest, dat = Test, model = mod1, bandwidth = tau * ntrain^(-1/5)/2)
    imseerr <- measure(train.hat2, Grid, t0 = 0, tau = tau, method = "imse",
                      L = train$L, R = train$R)
    interr <- measure(test.hat2, Grid, t0 = 0, tau = tau, method = "int",
                      surv.true = s.test)
    result2 <- 
      list(cox = mod1,
           predictedNO.Sm = train.hat2,
           imse.oob = matrix(c(imse.type1 = NaN, imse.type2 = NaN), ncol = 2),
           imse.NO = matrix(imseerr, ncol = 2),
           test = list(predicted = test.hat2, testerror = matrix(interr, ncol = 2)))
  }
  print(matrix(interr, ncol = 2))
  return(list(noSmooth = result1, smooth = result2))
}

errVec <- c("imse.type1 (oob)", "imse.type2 (oob)", "imse.type1", "imse.type2", 
            "int.error", "sup.error")
emptyErr <- function(nfold) {
  matrix(NA, nrow = nfold - 1, ncol = length(errVec), dimnames = list(2:nfold, errVec))
}
extrErr <- function(obj, nfold = dim(a)[1]) {
  a <- cbind(obj$imse.oob, obj$imse.NO, obj$test$testerr)
  if (dim(a)[1] < nfold)
    a <- rbind(a, emptyErr(nfold))
  colnames(a) <- errVec
  rownames(a) <- 1:nfold
  a
}
summaryEval <- function(obj) {
  a <- lapply(obj, extrErr, nfold)
  time <- sapply(obj, function(s) s$runtime$elapsed)
  structure(
    array(unlist(a), 
          dim = c(dim(a[[1]]), length(a)),
          dimnames = c(dimnames(a[[1]]), list(names(a)))),
    runtime = time)
}




plotRF <- function(rfobj, i = 1, pred = c("test", "train"),
                   truth = if (pred == "test") s.test[i,] else s.train[i,],
                   truth.time = Grid, tau = tau) {
  require(ggplot2)
  require(dplyr)
  pred <- match.arg(pred)

  if (pred == "test") {
    s = rfobj$test$predicted[i, ]
    time =  truth.time
  } else if (pred == "train") {
    s = rfobj$predicted[i, ]
    time = rfobj$time.points
  }

  if(is.null(truth)) {
    truth.time = NULL
    n.truth = 0
  } else {
    n.truth = length(truth)
    if (is.null(truth.time)) {
      truth.time = time
    }
  }

  data.frame(s = c(s, truth),
             time = c(time, truth.time),
             group = rep(c("est", "truth"), times = c(length(time), n.truth))) %>%
    ggplot(aes(time, s, col = group, group = group)) +
    geom_line() + xlim(c(0, tau))
}
