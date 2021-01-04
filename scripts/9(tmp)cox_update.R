#dashboard <- as.data.frame(expand.grid(sim = 1:100, scenario = 1:6, n.monitor = c(1,3)))
#dashboard$check <- dashboard$checkFn <- ""
#nruns <- dim(dashboard)[1]

# libraries
  library(icrf)
  library(icenReg)
  library(MASS)
  library(dplyr)
  source("0functions.R")
  source("1setting.R")
  { 
    ticksize = 0.01; ntest = 300        # test set size for evaluation
    ntree = 300L; nmin = 6; nfold = 10           # tree parameters
  }

# updates
#for (i in 2:3) {
for (i in 1:nruns) {
  scenario  = dashboard$scenario[i]
  sim       = dashboard$sim[i]
  n.monitor = dashboard$n.monitor[i]
print(i)  
  setting(scenario, sim, n.monitor, ntree, 0, date = "2019-08-11")
  
  if (!file.exists(fn_eval)) {
    dashboard$check[i] <- "no file yet"
    print(paste0(fn_eval, " file does not exit yet"))
    next
  }
  if (dashboard$check[i] == "updated") {
    print("Already updated.")
  } else {
    result.eval <- readRDS(fn_eval)
    if (all(result.eval[1,3:4,"cox"] == result.eval[1,5:6,"cox"])) {
      tmp.cox <- cox(form1)
      result.eval[1, , "cox"] <- extrErr(tmp.cox)
      saveRDS(result.eval, fn_eval)
    }
    dashboard$check[i] <- "updated"
    print(paste0(fn_eval, " updated"))
  }
    
  if (dashboard$checkFn[i] == "updated" && sim == 1) {
    print("Already updated.")
    next
  } 
  if (file.exists(fn_output)) {
    result <- readRDS(fn_output)
    result$cox <- tmp.cox
    saveRDS(result, fn_output)
    dashboard$checkFn[i] <- "updated"
    print(paste0(fn_output, " updated"))
    rm(result)
    gc()
  }
  rm(tmp.cox, result.eval)
  gc()
}
