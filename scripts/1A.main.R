### 1A.main.R
### Main simulations comparing icrf with other existing methods.

# library(devtools)
# load_all()
# install.packages("icrf")
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
  pilot     = as.numeric(args[4]) # if 1 pilot test, if 0 real simulation.
  print("test.R")
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
    ntree = 10L; nmin = 6; nmin.t = 20; nfold = 3           # tree parameters
  } else {
    ntree = 300L; nmin = 6; nmin.t = 20; nfold = 10          # tree parameters
  }
}
setting(scenario, sim, n.monitor, ntree, pilot)

if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

result <- list()

print("1. Wilcoxon's RF (GWRS)")
set.seed(seed.no + 1); result$wH <- rf(method = "Wilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
set.seed(seed.no + 1); result$wE <- rf(method = "Wilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)

print("2. logrank RF (GLR)")
set.seed(seed.no + 1); result$lH <- rf(method = "logrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
set.seed(seed.no + 1); result$lE <- rf(method = "logrank", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)

print("3. Peto's Wilcoxon RF (SWRS)")
set.seed(seed.no + 1); result$swH <- rf(method = "PetoWilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
set.seed(seed.no + 1); result$swE <- rf(method = "PetoWilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)

print("4. Peto's logrank RF (SLR)")
set.seed(seed.no + 1); result$slH <- rf(method = "PetoLogrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F)
set.seed(seed.no + 1); result$slE <- rf(method = "PetoLogrank", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F)

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
