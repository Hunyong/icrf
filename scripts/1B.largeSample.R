### 1B.largeSample.R
### Simulations comparing the accuracies under different sample sizes.

## 0. Settings
## 0.1 simulation parameters
  # rm(list = ls())
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  if (length(args) == 0) {
    warning("argument is not provided. Set as the default values.")
    args = c(n = 100, sim = 1, n.monitor = 1, pilot = 0, date = NA)
  }
  names(args) = c("n", "sim", "n.monitor", "pilot", "date")
  print(args)
  
  n         = as.numeric(args[1])
  sim       = as.numeric(args[2]) # replicates = 1..500
  n.monitor = as.numeric(args[3]) # number of censoring times = 1, 3
  pilot     = as.numeric(args[4]) # if 1 pilot, if 0 real simulation.
  date      = as.character(args[5]) # date for the output folder name
  if (is.na(date)| date == "0") date = NULL
  print("scripts/1B.largeSample.R")
  print(if (pilot) "pilot test!" else "real simulation!")

## 0.2 library
  library(icrf); library(icenReg); 
  library(MASS); library(dplyr); library(ggplot2)
  source("scripts/0functions.R")
  source("scripts/1setting.R")

## 0.3 rest of the settings
  ticksize = 0.01; ntest = 300        # test set size for evaluation
  #ticksize = tau/100   # size of grid for evaluation
  if (pilot == 1) {
    ntree = 10L; nmin = 6; nmin.t = 20; nfold = 2           # tree parameters
  } else {
    ntree = 100L; nmin = round(n^0.7/12 + 1, 0) ; nmin.t = round(n^0.7/3 + 2, 0); nfold = 10
  }
  scenario  = 1
  n.sim = 300

  setting(scenario = scenario, sim, n.monitor, ntree, pilot, n = n, simClass = "size", date = date)

  if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

  
## 1. simulations
  result <- list()
  
  print("1. Wilcoxon's RF (GWRS)")
  set.seed(seed.no + 1); result$w132 <- rf(split.rule = "Wilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F); gc()
  
  print("5. - 8. Fu's trees - TR1, TR2, RF1, RF2")
  set.seed(seed.no + 1); result$FuTR1 <- Fu(RF = F, smoothing = F); gc()
  set.seed(seed.no + 1); result$FuTR2 <- Fu(RF = F, smoothing = T); gc()
  set.seed(seed.no + 1); result$FuRF1 <- Fu(RF = T, smoothing = F); gc()
  set.seed(seed.no + 1); result$FuRF2 <- Fu(RF = T, smoothing = T); gc()
  #plotRF(result$FuTR1, i = 1)
  
  print("Cox regression")
  result$cox <- cox(form1, smooth = FALSE); gc()
  result$cox.sm <- cox(form1, smooth = TRUE); gc()
  
## 2. evaluation and saving
  print("Evaluation and saving.")
  result.eval <- summaryEval(result)
  saveRDS(result.eval, fn_eval)
  
  if (sim == 1) saveRDS(result, fn_output)
  Sys.time() - time.bgn  # time elapsed.
