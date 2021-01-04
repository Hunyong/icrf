library(icrf); library(icenReg); 
library(MASS); library(dplyr); library(ggplot2)
source("0functions.R")
source("1setting.R")

n.monitor = 1; pilot = 0
ntree = 300L; nmin = 6; nmin.t = 20; nfold = 10
ticksize = 0.01; ntest = 300        # test set size for evaluation
n.sim = 100
for (scenario in 1:6) {
  for (sim in 1:100) {
    cat("scenario = ", scenario, ", sim = ", sim,"\n")
    setting(scenario, sim, n.monitor, ntree, pilot)
    
    # updating the output
    if (file.exists(fn_output)) {
      result <- readRDS(fn_output)
      svRDS <- TRUE
    } else {
      result <- list()
      svRDS <- FALSE
    }
    set.seed(seed.no + 1); result$FuTR1 <- Fu(RF = F, smoothing = F)
    set.seed(seed.no + 1); result$FuTR2 <- Fu(RF = F, smoothing = T)
    if (svRDS) saveRDS(result, fn_output)
    
    # updating the evaluation
    result.eval <- readRDS(fn_eval)
    result.eval[,, "FuTR1"] <- summaryEval(result)[,, "FuTR1"]
    result.eval[,, "FuTR2"] <- summaryEval(result)[,, "FuTR2"]
    saveRDS(result.eval, fn_eval)
  }
}

if (file.exists(fn_eval)) stop(paste0(fn_eval, " file already exits"))

result <- list()

print("1. Wilcoxon's RF - 132")
set.seed(seed.no + 1); result$w132 <- rf(method = "Wilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.226 / 0.209

print("2. logrank RF - 132")
set.seed(seed.no + 1); result$l132 <- rf(method = "logrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.224 / 0.214

print("3. Peto's Wilcoxon RF - 132")
set.seed(seed.no + 1); result$pw132 <- rf(method = "PetoWilcoxon", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.237 / 0.206

print("4. Peto's logrank RF - 132")
set.seed(seed.no + 1); result$pl132 <- rf(method = "PetoLogrank", quasihonesty = T, ERT = T, sampsize = ntrain * 0.95, replace = F) # 0.237 / 0.206

print("5. - 8. Fu's trees - TR1, TR2, RF1, RF2")
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
