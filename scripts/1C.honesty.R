### 1C.honesty.R
### Simulations comparing quasi-honesty vs. exploitative approaches.

{
  # rm(list = ls())
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  if (length(args) == 0) {
    warning("argument is not provided. Set as the default values.")
    args = c(scenario = 1, sim = 1, n.monitor = 1, pilot = 0)
  }
  names(args) = c("scenario", "sim", "n.monitor", "pilot")
  print(args)
  scenario  = as.numeric(args[1]) # scenario = 1, 2, 3, 4, 5, 6
  sim       = as.numeric(args[2]) # replicates = 1..500
  n.monitor = as.numeric(args[3]) # number of censoring times = 1, 3
  pilot     = as.numeric(args[4]) # if 1 pilot, if 0 real simulation.
  print("1C.honesty.R") # to see if partial honesty works.
}
library(icrf); library(icenReg); 
library(MASS); library(dplyr); library(ggplot2)
source("scripts/0functions.R")
source("scripts/1setting.R")

{
  ticksize = 0.01; ntest = 300        # test set size for evaluation
  #ticksize = tau/100   # size of grid for evaluation
  n.sim = 300
  if (pilot == 1) {
    ntree = 10L; nmin = 6; nmin.t = 20; nfold = 2           # tree parameters
  } else {
    ntree = 300L; nmin = 6; nmin.t = 20; nfold = 10          # tree parameters
  }
}
setting(scenario, sim, n.monitor, ntree, pilot, simClass = "honesty")

if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

result <- list()

print("1. Wilcoxon's RF (GWRS), non-honest (exploitative) trees")  #updateNPMLE = TRUE: quasihonesty, FALSE : exploitative
set.seed(seed.no + 1); result$wE <- rf(split.rule = "Wilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209

print("Evaluation and saving.")
result.eval <- summaryEval(result)
saveRDS(result.eval, fn_eval)

if (sim == 1) saveRDS(result, fn_output)
Sys.time() - time.bgn  # time elapsed.
