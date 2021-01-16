### 1B.largeSample.R
### Simulations comparing the accuracies under different sample sizes.

# library(devtools)
# load_all()
# install.packages("icrf")
{
  # rm(list = ls())
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  if (length(args) == 0) {
    warning("argument is not provided. Set as the default values.")
    args = c(n = 100, sim = 1, n.monitor = 1, pilot = 0)
  }
  names(args) = c("n", "sim", "n.monitor", "pilot")
  print(args)
  n         = as.numeric(args[1]) # number of censoring times = 1, 3
  sim       = as.numeric(args[2]) # replicates = 1..500
  n.monitor = as.numeric(args[3])
  pilot     = as.numeric(args[4]) # if 1 pilot, if 0 real simulation.
  print("scripts/1B.largeSample.R")
  print(if (pilot) "pilot test!" else "real simulation!")
}
library(icrf); library(icenReg); 
library(MASS); library(dplyr); library(ggplot2)
source("scripts/0functions.R")
source("scripts/1setting.R")

{ 
  ticksize = 0.01; ntest = 300        # test set size for evaluation
  #ticksize = tau/100   # size of grid for evaluation
  n.sim = 100
  if (pilot == 1) {
    ntree = 10L; nmin = 6; nmin.t = 20; nfold = 2           # tree parameters
  } else {
    ntree = 300L; nmin = 6; nmin.t = 20; nfold = 10          # tree parameters
  }
  scenario  = 1
}
setting(scenario = 1, sim, n.monitor, ntree, pilot, n = n, largeSample = TRUE)

if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

result <- list()

print("1. Wilcoxon's RF (GWRS)")
set.seed(seed.no + 1); result$w132 <- rf(method = "Wilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209
set.seed(seed.no + 1); result$w135 <- rf(method = "Wilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209

print("5. - 8. Fu's trees - TR1, TR2, RF1, RF2")
set.seed(seed.no + 1); result$FuTR1 <- Fu(RF = F, smoothing = F)
set.seed(seed.no + 1); result$FuTR2 <- Fu(RF = F, smoothing = T)
set.seed(seed.no + 1); result$FuRF1 <- Fu(RF = T, smoothing = F)
set.seed(seed.no + 1); result$FuRF2 <- Fu(RF = T, smoothing = T)
#plotRF(result$FuTR1, i = 1)

print("Cox regression")
result$cox <- cox(form1, smooth = FALSE)
result$cox.sm <- cox(form1, smooth = TRUE)


print("Evaluation and saving.")
result.eval <- summaryEval(result)
saveRDS(result.eval, fn_eval)

if (sim == 1) saveRDS(result, fn_output)
Sys.time() - time.bgn  # time elapsed.
