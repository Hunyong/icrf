{
  rm(list = ls())
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  names(args) = c("scenario", "sim", "n.monitor", "pilot")
  # args = c(1, 1, 1, 0)
  print(args)
  scenario  = as.numeric(args[1]) # scenario = 1, 2, 3, 4, 5, 6
  sim       = as.numeric(args[2]) # replicates = 1..500
  n.monitor = as.numeric(args[3]) # number of censoring times = 1, 3
  pilot     = as.numeric(args[4]) # if 1 pilot, if 0 real simulation.
  print("1C.honesty.R") # to see if partial honesty works.
}
library(icrf); library(icenReg); 
library(MASS); library(dplyr); library(ggplot2)
source("0functions.R")
source("1setting.R")
# source("cox_LASSO/em.func.R")
# source("cox_LASSO/em.func.sub.R")

{
  ticksize = 0.01; ntest = 300        # test set size for evaluation
  #ticksize = tau/100   # size of grid for evaluation
  n.sim = 100
  if (pilot == 1) {
    ntree = 10L; nmin = 6; nmin.t = 20; nfold = 2           # tree parameters
  } else {
    ntree = 300L; nmin = 6; nmin.t = 20; nfold = 10          # tree parameters
  }
}
setting(scenario, sim, n.monitor, ntree, pilot, nonhonesty = TRUE)

if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

result <- list()

print("1. Wilcoxon's RF - 135: non-honest trees")  #updateNPMLE = TRUE: quasihonesty, FALSE : exploitative
set.seed(seed.no + 1); result$w135 <- rf(method = "Wilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209

print("Evaluation and saving.")
result.eval <- summaryEval(result)
saveRDS(result.eval, fn_eval)

if (sim == 1) saveRDS(result, fn_output)
Sys.time() - time.bgn  # time elapsed.