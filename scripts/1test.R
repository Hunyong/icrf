# library(devtools)
# load_all()
# setwd("scripts")
# install.packages("/Users/hycho/Documents/1Research/Rpackages/icrf_0.1.2.tar.gz", repos = NULL, INSTALL_opts = c('--no-lock'))
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
  # scenario = 5; sim = 79; n.monitor = 1; pilot = 0
  print("test.R")
  print(if (pilot) "pilot test!" else "real simulation!")
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
    ntree = 10L; nmin = 6; nmin.t = 20; nfold = 3           # tree parameters
  } else {
    ntree = 300L; nmin = 6; nmin.t = 20; nfold = 10          # tree parameters
  }
}
setting(scenario, sim, n.monitor, ntree, pilot)

if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

result <- list()

print("1. Wilcoxon's RF - 135")
set.seed(seed.no + 1); result$w132 <- rf(method = "Wilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209
set.seed(seed.no + 1); result$w135 <- rf(method = "Wilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209

print("2. logrank RF - 135")
set.seed(seed.no + 1); result$l132 <- rf(method = "logrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.224 / 0.214
set.seed(seed.no + 1); result$l135 <- rf(method = "logrank", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.224 / 0.214

print("3. Peto's Wilcoxon RF - 135")
set.seed(seed.no + 1); result$pw132 <- rf(method = "PetoWilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.237 / 0.206
set.seed(seed.no + 1); result$pw135 <- rf(method = "PetoWilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.237 / 0.206

print("4. Peto's logrank RF - 135")
set.seed(seed.no + 1); result$pl132 <- rf(method = "PetoLogrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.237 / 0.206
set.seed(seed.no + 1); result$pl135 <- rf(method = "PetoLogrank", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.237 / 0.206

print("5. - 8. Fu's trees - TR1, TR2, RF1, RF2")
set.seed(seed.no + 1); result$FuTR1 <- Fu(RF = F, smoothing = F)
set.seed(seed.no + 1); result$FuTR2 <- Fu(RF = F, smoothing = T)
set.seed(seed.no + 1); result$FuRF1 <- Fu(RF = T, smoothing = F)
set.seed(seed.no + 1); result$FuRF2 <- Fu(RF = T, smoothing = T)
#plotRF(result$FuTR1, i = 1)

print("Cox regression")
result$cox <- cox(form1, smooth = FALSE)
result$cox.sm <- cox(form1, smooth = TRUE)

# source("cox_LASSO/em.func.R")
# source("cox_LASSO/em.func.sub.R")
# print("Cox-LASSO")
# tmp.a <- coxLASSO(form1, smooth = "both")
# result$cox.LS <- tmp.a$noSmooth
# result$cox.SL.sm <- tmp.a$smooth


print("Evaluation and saving.")
result.eval <- summaryEval(result)
saveRDS(result.eval, fn_eval)

if (sim == 1) saveRDS(result, fn_output)
Sys.time() - time.bgn  # time elapsed.
