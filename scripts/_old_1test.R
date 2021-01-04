  ### Get arguments here!!!


  # library(devtools)
  # load_all()
  # setwd("scripts")
  # install.packages("/Users/hycho/Documents/1Research/Rpackages/icrf_0.1.0.tar.gz", repos = NULL, INSTALL_opts = c('--no-lock'))
  {
    rm(list = ls())
    args = commandArgs(trailingOnly=TRUE)  # passed from script
    print(args)
    scenario  = as.numeric(args[1]) # scenario = 1, 2, 3, 4, 5, 6
    sim       = as.numeric(args[2]) # replicates = 1..500
    n.monitor = as.numeric(args[3]) # number of censoring times = 1, 3
    pilot     = as.numeric(args[4]) # if 1 pilot, if 0 real simulation.
    # scenario = 5; sim = 79; n.monitor = 1; pilot = 0
    print("test.R")
    print(if (pilot) "pilot test!" else "real simulation!")
  }
  library(icrf)
  library(ggplot2)
  library(icenReg)
  library(MASS)
  library(dplyr)
  source("0functions.R")
  source("1setting.R")
  source("cox_LASSO/em.func.R")
  source("cox_LASSO/em.func.sub.R")

  { 
    ticksize = 0.01; ntest = 300        # test set size for evaluation
    #ticksize = tau/100   # size of grid for evaluation
    n.sim = 100
    if (pilot == 1) {
      ntree = 2L; nmin = 6; nfold = 2           # tree parameters
    } else {
      ntree = 300L; nmin = 6; nfold = 10           # tree parameters
    }
  }
  setting(scenario, sim, n.monitor, ntree, pilot)
  
  if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))
  #print(readRDS(fn2))

  # ntrain = 30 vs 100
  result <- list()
  
  print("Wilcoxon's RF - 111, 112, 120, 131, 132")
  set.seed(seed.no + 1); result$w111 <- rf(method = "wilcoxon", updateNPMLE = T, ERT = T, uniformERT = T) # 0.223 / 0.214
  set.seed(seed.no + 1); result$w112 <- rf(method = "wilcoxon", updateNPMLE = T, ERT = T, uniformERT = F) # 0.267 / 0.256
  set.seed(seed.no + 1); result$w120 <- rf(method = "wilcoxon", updateNPMLE = T, ERT = F, uniformERT = -1) # 0.307 / 0.231
  set.seed(seed.no + 1); result$w131 <- rf(method = "wilcoxon", updateNPMLE = T, ERT = T, sampsize = ntrain * 0.90, replace = F) # 0.226 / 0.209
  set.seed(seed.no + 1); result$w132 <- rf(method = "wilcoxon", updateNPMLE = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209
  
  print("Peto's RF - 111, 112, 120, 131, 132")
  set.seed(seed.no + 1); result$p111 <- rf(method = "Peto", updateNPMLE = T, ERT = T, uniformERT = T) # 0.228 / 0.212
  set.seed(seed.no + 1); result$p112 <- rf(method = "Peto", updateNPMLE = T, ERT = T, uniformERT = F) # 0.253 / 0.235
  set.seed(seed.no + 1); result$p120 <- rf(method = "Peto", updateNPMLE = T, ERT = F, uniformERT = -1) # 0.231 / 0.205 ****
  set.seed(seed.no + 1); result$p131 <- rf(method = "Peto", updateNPMLE = T, ERT = T, sampsize = ntrain * 0.90, replace = F) # 0.237 / 0.206
  set.seed(seed.no + 1); result$p132 <- rf(method = "Peto", updateNPMLE = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.237 / 0.206
  

  print("logrank RF - 111, 112, 120, 131, 132")
  set.seed(seed.no + 1); result$L111 <- rf(method = "log", updateNPMLE = T, ERT = T, uniformERT = T) # 0.211*** / 0.206 ***
  set.seed(seed.no + 1); result$L112 <- rf(method = "log", updateNPMLE = T, ERT = T, uniformERT = F) # 0.294 / 0.257
  set.seed(seed.no + 1); result$L120 <- rf(method = "log", updateNPMLE = T, ERT = F, uniformERT = -1) # 0.236 / 0.209
  set.seed(seed.no + 1); result$L131 <- rf(method = "log", updateNPMLE = T, ERT = T, sampsize = ntrain * 0.90, replace = F) # 0.224 / 0.214
  set.seed(seed.no + 1); result$L132 <- rf(method = "log", updateNPMLE = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.224 / 0.214
  
  print("Fu's trees - TR1, TR2, RF1, RF2")
  set.seed(seed.no + 1); result$FuTR1 <- Fu(RF = F, smoothing = F) # 0.785
  set.seed(seed.no + 1); result$FuTR2 <- Fu(RF = F, smoothing = T) # 0.742
  set.seed(seed.no + 1); result$FuRF1 <- Fu(RF = T, smoothing = F) # 0.205
  set.seed(seed.no + 1); result$FuRF2 <- Fu(RF = T, smoothing = T) # 0.204
#plotRF(result$FuTR1, i = 1)
  
  print("Cox regression")
  result$cox <- cox(form1, smooth = FALSE) # 0.224 / 0.214
  result$cox.sm <- cox(form1, smooth = TRUE) # 0.224 / 0.214
  
  source("cox_LASSO/em.func.R")
  source("cox_LASSO/em.func.sub.R")
  print("Cox-LASSO")
  # result$cox <- cox(form1) # 0.224 / 0.214
  a3 <- coxLASSO(form1)
  a3 <- coxLASSO(form1)
  
  
  print("Evaluation and saving.")
  result.eval <- summaryEval(result)
  saveRDS(result.eval, fn_eval)
  
  if (sim == 1) saveRDS(result, fn_output)
  Sys.time() - time.bgn  # time elapsed.
