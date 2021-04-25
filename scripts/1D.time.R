### 1D.computation.R
### Measuring the computation time

# install.packages("icrf")
{
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  if (length(args) == 0) {
    warning("argument is not provided. Set as the default values.")
    args = c(n = 100, ntree = 100, nfold = 10, sim = 1)
  }
  names(args) = c("n", "ntree", "nfold", "sim")
  print(args)
  n         = as.numeric(args[1])
  ntree     = as.numeric(args[2])
  nfold     = as.numeric(args[3])
  sim       = as.numeric(args[4]) # replicates = 1..500
  print("scripts/1D.computation.R")
}
library(icrf); library(icenReg); 
library(MASS); library(dplyr); library(ggplot2)
source("scripts/0functions.R")
source("scripts/1setting.R")

scenario = 1; n.monitor = 1; ntest = 300; n.sim = 300
nmin = 6; nmin.t = 20         # tree parameters
setting(scenario = scenario, sim = sim, n.monitor =  n.monitor, ntree = ntree, pilot = 0, n = n, simClass = "time")

if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

result <- list()

print("1. Wilcoxon's RF (GWRS)")
set.seed(seed.no + 1); result$wH <- rf(method = "Wilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F); gc()
set.seed(seed.no + 1); result$wE <- rf(method = "Wilcoxon", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F); gc()
set.seed(seed.no + 1); result$lH <- rf(method = "logrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F); gc()
set.seed(seed.no + 1); result$lE <- rf(method = "logrank", quasihonesty = F, ERT = T, sampsize = ntrain * 0.95, replace = F); gc()

print("5. - 8. Fu's trees - TR1, TR2, RF1, RF2")
set.seed(seed.no + 1); result$FuTR1 <- Fu(RF = F, smoothing = F); gc()
set.seed(seed.no + 1); result$FuTR2 <- Fu(RF = F, smoothing = T); gc()
set.seed(seed.no + 1); result$FuRF1 <- Fu(RF = T, smoothing = F); gc()
set.seed(seed.no + 1); result$FuRF2 <- Fu(RF = T, smoothing = T); gc()

print("Evaluation and saving.")
result.eval <- 
  sapply(result, function(x) x$runtime["elapsed"]) %>% # extracting runtimes
  unlist
result.eval <-
  data.frame(n = n, ntree = ntree, nfold = nfold, sim = sim,
             method = names(result.eval) %>% gsub("\\.elapsed", "", .), 
             time = result.eval)
print(result.eval)
saveRDS(result.eval, fn_eval)

cat("done.")
