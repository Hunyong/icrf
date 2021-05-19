### 1A.main.R
### Main simulations comparing icrf with other existing methods.

## 0. Settings
## 0.1 simulation parameters
  # rm(list = ls())
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  if (length(args) == 0) {
    warning("argument is not provided. Set as the default values.")
    args = c(scenario = 1, sim = 1, n.monitor = 1, pilot = 0, date = NA)
  }
  names(args) = c("scenario", "sim", "n.monitor", "pilot", "date")
  print(args)
  
  scenario  = as.numeric(args[1]) # scenario = 1, 2, 3, 4, 5, 6
  sim       = as.numeric(args[2]) # replicates = 1..500
  n.monitor = as.numeric(args[3]) # number of censoring times = 1, 3
  pilot     = as.numeric(args[4]) # if 1 pilot test, if 0 real simulation.
  date      = as.character(args[5]) # date for the output folder name
  if (is.na(date)| date == "0") date = NULL
  print("1A.main.R")
  print(if (pilot) "pilot test!" else "real simulation!")

## 0.2 library
  library(icrf); library(icenReg); 
  library(MASS); library(dplyr); library(ggplot2)
  source("scripts/0functions.R")
  source("scripts/1setting.R")

## 0.3 rest of the settings
  ticksize = 0.01; ntest = 300        # test set size for evaluation
  #ticksize = tau/100   # size of grid for evaluation
  n.sim = 300
  if (pilot == 1) {
    ntree = 10L; nmin = 6; nmin.t = 20; nfold = 3           # tree parameters
  } else {
    ntree = 300L; nmin = 6; nmin.t = 20; nfold = 10          # tree parameters
  }
  setting(scenario, sim, n.monitor, ntree, pilot, date = date, simClass = "main")

  if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))


## 1. simulations
  result <- list()
  
  print("1. Wilcoxon's RF (GWRS)")
  set.seed(seed.no + 1); result$wH <- rf(split.rule = "Wilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
  set.seed(seed.no + 1); result$wE <- rf(split.rule = "Wilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)
  
  print("2. logrank RF (GLR)")
  set.seed(seed.no + 1); result$lH <- rf(split.rule = "logrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
  set.seed(seed.no + 1); result$lE <- rf(split.rule = "logrank", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)
  
  print("3. Peto's Wilcoxon RF (SWRS)")
  set.seed(seed.no + 1); result$swH <- rf(split.rule = "PetoWilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
  set.seed(seed.no + 1); result$swE <- rf(split.rule = "PetoWilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)
  
  print("4. Peto's logrank RF (SLR)")
  set.seed(seed.no + 1); result$slH <- rf(split.rule = "PetoLogrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
  set.seed(seed.no + 1); result$slE <- rf(split.rule = "PetoLogrank", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)
  
  print("5. - 8. Fu's trees - TR1, TR2, RF1, RF2")
  set.seed(seed.no + 1); result$FuTR1 <- Fu(RF = F, smoothing = F)
  set.seed(seed.no + 1); result$FuTR2 <- Fu(RF = F, smoothing = T)
  set.seed(seed.no + 1); result$FuRF1 <- Fu(RF = T, smoothing = F)
  set.seed(seed.no + 1); result$FuRF2 <- Fu(RF = T, smoothing = T)
  #plotRF(result$FuTR1, i = 1)
  
  print("Cox regression")
  result$cox <- cox(form1, smooth = FALSE)
  result$cox.sm <- cox(form1, smooth = TRUE)
  
## 2. Evaluation and saving
  print("Evaluation and saving.")
  result.eval <- summaryEval(result)
  
  saveRDS(result.eval, fn_eval)
  if (sim == 1) saveRDS(result, fn_output)