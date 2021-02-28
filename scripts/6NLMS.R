# The National Longitudinal Mortality Study (NLMS) is a national, longitudinal, 
# mortality study sponsored by the National Heart, Lung, and Blood Institute, 
# the National Cancer Institute, the National Institute on Aging (all part of 
# the National Institutes of Health), the National Center for Health Statistics 
# part of the Center for Disease Control and Prevention  and the U.S. Census 
# Bureau for the purpose of studying the effects of differentials in demographic 
# and socio-economic characteristics on mortality.

# The data are available through
# https://biolincc.nhlbi.nih.gov/studies/nlms/


### 0. library, set up, and data
  # 0.0 library
  #install.packages("/Users/hunyongcho/Documents/1Research/Rpackages/icrf_1.0.2.tar.gz", repos = NULL, INSTALL_opts = c('--no-lock'))
  library(dplyr)
  library(icrf)
  library(survival)
  source("scripts/0functions.R")
  
  # 0.1 set up
  pilot = FALSE
  
  out_path = "output/nlms"
  if (!dir.exists(out_path)) dir.create(out_path)
  fig_path = "figure/nlms"
  if (!dir.exists(fig_path)) dir.create(fig_path)
  
  i = as.numeric(commandArgs(trailingOnly=TRUE))[1]  # passed from script
  if (file.exists(paste0(out_path, "/nlms_eval_", i,".rds"))) stop("Already done.")
  
  if (is.na(i)) {
    warning("i (the cross-validation index) was not provided. A full loop will be implemented.")
    i = 1; rng = 1:300; parallel = FALSE
  } else {
    rng = i:i; parallel = TRUE
    cat("i = ", i, "\n")
  }
   
  if (pilot) {
    samp.size = 0.1
    ntree= 2
    nfold= 2
    rng = i:i # 1.
  } else {
    samp.size = 0.7 # 3630 * 0.7 = 1815
    ntree = 50     # preliminary. Later increase to ntree = 300, nfold = 10
    nfold = 5
  }
  
  # 0.2 reading the data
  nlms.raw = read.csv("NLMS_PublicUse_Release5b/6c.csv")
  nlms.raw$follow %>% table
  
  # returns TRUE only when there is at least one level not observed in the data.
  check.levels = function(x) is.factor(x) & !all(levels(x) %in% unique(x)) 
  drop.levels = function(x) if (is.factor(x)) droplevels(x) else x
  # na_to_zero = function(x, to_replace = 0) ifelse(is.na(x), to_replace, x)
  tab = function(x)
    if (is.numeric(x) & length(unique(x)) > 15) {
      summary(x)
    } else {
      table(x, useNA = "always")
    }
  nlms.summary = lapply(nlms.raw, tab)
  names(nlms.summary)
  
  nlms =
    nlms.raw %>% 
    transmute(##1. Survival status
              delta = inddea,                            # No missing. 0 = Alive, 1 = Death
              time = follow,                             # No missing. 0...2192. 2192 = Administrative censoring date
                      # table(nlms$follow == 2192, nlms$inddea) # indicates one person died on the final censoring day.
              # cause113                                 # No missing. Cause of death. 0,1,2,...,113. This does not constitute the regression data.
              # dayod (day of week of death), hosp (hospital type), hospd (hospital death indicator), indalg (Indicator for Algorithmic Death)
              #   are not relevant to the regression.
              # L = time, R = ifelse(delta, time, Inf),    # right-censored data
  
              ## 2. demographic
              age, 
              sex,
              # race = race %>% factor,     # NA = unknown, 1,2,3,4,5
              race = ifelse(race > 2, 3, race) %>% as.factor, # combining minorities.
              hisp = hisp %>% factor,     # NA = unknown, 1,2,3,4
              hhnum,                                     # number of people in household
              reltrf = reltrf %>% factor,                # 1,2,3,4,5,6. No missing. Relation to reference person within the household.
              # hhid # household id is not considered in this regression
              
              ## Economic, Educational
              adjinc = adjinc %>% ordered,# NA = unknown (0.1%)   # Adjusted income. Due to small fraction of missing, put it as ordinal.
              # povpct = povpct %>% ordered,# Income as Percent of Poverty Level. Missing fraction 2.6%.
              #                               # Almost collinear with adjinc
              ssnyn,                                     # 0,1. No missing. Presence of SSN,
              tenure = tenure %>% factor, # 1,2,3. No missing. Housing Tenure.
              
              # Health / Health insurance
              wt,                                        # No missing. 0...4446. weight
              health  = health %>% factor,# 1,2,3,4,5. NA = missing. Health (Would you say the person's health is) 1 Excellent, ... 5 Poor
                                                              # Due to high chunk of missing (10%), treat it as nominal rather than ordinal.
              # histatus,                                  # 0,1. No missing. Health insurance status. Almost explained by hitype (i.e., can be aliased)
              hitype = hitype %>% factor,                # 0,1,...,5. No missing. Health insurance type. Medicare, ...
              
              ## Regional / Citizenship
              urban = urban %>% factor,   # 1,2,NA. Urban/Rural status
              citizen = ifelse(citizen == 1, 1, 2), # 1,2,3,4,5. NA = missing. Citizenship. Highly unbalanced. Dichotomized into citizen vs. non-citizen.
              # pob  = pob %>% pmin(900) %>% factor,       # 0 = missing. 101, 102, ..., 111, 901, ..., 956  # Region of Birth
              #                                              # pob greater than 900 is set as 900 (# levels >  53 the max capacity of RF.)
              #                                             # highly concentrated in 900, this induces collinarity with other variables such as citizenship.
              # smsast = smsast %>% factor, # 1,2,3,NA. Is a household located in SMSA (+ central or not) or not? # Almost collinear to urban.
              
              # Below are highly missing (%complete = 0.49 rcow, 0.51 majocc & occ & ind & majind, 0.06 indalg, 0.75 esr & ms & educ, 0.73 vt)
              # occ, ind # the occupation and industry code are too fine (505 levels). Use majocc and majind instead.
              # rcow   = rcow %>% factor,   # 1,2,3,4,5. NA = missing. Class of Worker. Private, government, self, ...
              # majocc = majocc %>% factor, # 1,2,...,11 NA = missing. Major Occupation Code
              # majind = majind %>% factor, # 1,2,...,10,14 NA = missing. Major Industry Code
              # stater = stater %>% factor,                # 11,12,...,51. No missing. State Recode
              # esr    = esr %>% factor,               # 1,2,3,4,5, NA = missing. Employment status recode.
              # vt     = vt %>% na_to_zero(-1) %>% factor  # 0,1. NA = missing. Veteran status.
              # ms   = ms %>% factor,       # NA = unknown, 1,2,3,4,5 # Marital status
              # educ = educ %>% ordered,    # NA = missing or children, 1,2,...,14
              # 
              
              )
  nlms %>% dim # 745,162 x 14
  nlms.complete = nlms %>% na.omit %>% dplyr::filter(age > 80)
  nlms.complete.n = nlms.complete %>% NROW %>% print # 574,795 x 15
  tau = 1500   # The actual study end was 2192. However, for the data we analyze, the failure is very sparse after 1500.
  sapply(nlms.complete, check.levels) # hitype, reltrf
  for (i in 1:dim(nlms.complete)[2]) nlms.complete[[i]] <- drop.levels(nlms.complete[[i]])
  sapply(nlms.complete, nlevels)
  sapply(nlms.complete, check.levels)
  
  nlms.Grid <- c(0, exp(seq(0, log(tau), length.out = 200)), Inf)
  nlms$delta %>% mean # 0.033 # High censoring rate.
  nlms.complete$delta %>% mean # 0.30 # moderate censoring rate.
  
  # nlms.eval.array.i <- 
  #   array(NA, dim = c(8, 3), 
  #         dimnames = list(method  = c("ICRF.H", "icrf.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "Cox1", "Cox2"),
  #                         measure = c("imse.type1", "imse.type2", "runtime")))
  
  

  n.success = 0
### 1. cross-validation
  for (i in rng) {
    cat("\nreplicate i = ", i, " ")  
    a <- Sys.time()
    # 0.4 hold-out sample
    set.seed(i)
    
    # samp1 = train set index
    # samp1 <- sample(1:nlms.complete.n, samp.size) %>% sort
    samp1 <- 
      nlms.complete %>% 
      mutate(no = 1:n()) %>%
      group_by(race, hisp, hitype) %>% # stratified sampling by the important variables
      sample_frac(size = samp.size) %>% 
      "$"(no)
    lack <- sapply(nlms.complete[samp1, ], check.levels) 
    # lack %T>% print %>% any
     
    if (any(lack)) {
      print("-  Not all levels are not present in train set. Go to the next i.")    # The Cox model fails to handle this (Problematic when predicting data with full levels).
      print(lack %>% which)
      saveRDS(nlms.eval.array.i, paste0(out_path, "/nlms_eval_", i,".rds"))
      if (parallel) stop("done.") else next
    } else {
      n.success = n.success + 1
      cat(sprintf("- %sth effective replicate (after skipping the partial samples).", n.success))
    }
    cat("\n")

    # samp2 = test set index
    samp2 <- which(!1:nlms.complete.n %in% samp1)
    samp2.n <- samp2 %>% length # 339
    
    # train and test sets
    nlms.train <- 
      nlms.complete[samp1, ] %>% 
      cens.nlms(n.monitor = 1, tau = tau, remove.original = TRUE) 
# %>% 
#   dplyr::filter(!(L == 0 & is.infinite(R)))  # Removing non-informative subjects, if any???
    
    nlms.test  <- 
      nlms.complete[samp2, ] %>% 
      mutate(L = time, R = ifelse(delta, time, Inf))
    # nlms.train.list <- lapply(samp1.list, function(s) nlms.complete[s, ])
    
    
    args.H = 
      args.E =
      list(# ~ ., data = nlms.train, data.type = "right", right.label = c("time", "delta"),
           ~ ., data = nlms.train, data.type = "interval", interval.label = c("L", "R"),
           tau = tau, proximity = F, importance = FALSE, nPerm = 10,
           nfold = nfold, ntree = ntree, nodesize = 6, mtry = 4,         # tree structure
           replace = F, sampsize = dim(nlms.train)[1] * .95,    # resampling 
           method = "Wilcoxon", ERT = TRUE, uniformERT = TRUE,      # splitting rules
           quasihonesty = TRUE, timeSmooth = nlms.Grid)
    args.E$quasihonesty = FALSE
    
    args.FuTR1 = 
      args.FuTR2 =
      list(#~ ., data = nlms.train, data.type = "right", right.label = c("time", "delta"),
           ~ ., data = nlms.train, data.type = "interval", interval.label = c("L", "R"), 
           tau = tau, proximity = F, importance = FALSE, nPerm = 10,
           nfold = 1, ntree = 1, nodesize = 20, mtry = dim(nlms.train)[2] - 2,         # tree structure
           replace = F, sampsize = dim(nlms.train)[1],    # resampling 
           method = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,      # splitting rules
           quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = 0, 
           timeSmooth = nlms.Grid)
    args.FuTR2$bandwidth = NULL # smoothing
    
    args.FuRF1 = 
      args.FuRF2 =
      list(#~ ., data = nlms.train, data.type = "right", right.label = c("time", "delta"),
           ~ ., data = nlms.train, data.type = "interval", interval.label = c("L", "R"),
           tau = tau, proximity = F, importance = FALSE, nPerm = 10,
           nfold = 1, ntree = ntree, nodesize = 6, mtry = 4,         # tree structure
           replace = T, sampsize = dim(nlms.train)[1] * .632,    # resampling 
           method = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,      # splitting rules
           quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = 0,
           timeSmooth = nlms.Grid)
    args.FuRF2$bandwidth = NULL # smoothing
    
    if (i == 1) { # Measure the variable importance only for the first run.
      args.H$importance = args.E$importance = args.FuTR1$importance = 
        args.FuTR2$importance = args.FuRF1$importance = args.FuRF2$importance = TRUE
    }
    
    cat("1.icrf.H\n");set.seed(i); nlms.icrf.H <- do.call(icrf, args.H)
    cat("2.icrf.E\n");set.seed(i); nlms.icrf.E <- do.call(icrf, args.E)
    cat("3.FuTR1\n"); set.seed(i); nlms.FuTR1  <- do.call(icrf, args.FuTR1)
    cat("4.FuTR2\n"); set.seed(i); nlms.FuTR2  <- do.call(icrf, args.FuTR2)
    cat("5.FuRF1\n"); set.seed(i); nlms.FuRF1  <- do.call(icrf, args.FuRF1)
    cat("6.FuRF2\n"); set.seed(i); nlms.FuRF2  <- do.call(icrf, args.FuRF2)
    
    cat("7.Cox\n");   
    nlms.cox <- 
      icenReg::ic_sp(Surv(L, R, type = 'interval2') ~ ., model = "ph", 
                     data = nlms.train)
    
    if (i %in% 1:10) {
      saveRDS(nlms.icrf.H, sprintf("%s/%snlms_ICRF.H.%s.rds", out_path, ifelse(pilot, "pilot_", ""), i))
      saveRDS(nlms.icrf.E, sprintf("%s/%snlms_ICRF.E.%s.rds", out_path, ifelse(pilot, "pilot_", ""), i))
      saveRDS(nlms.FuTR1, sprintf("%s/%snlms_FuTR1.%s.rds", out_path, ifelse(pilot, "pilot_", ""), i))
      saveRDS(nlms.FuTR2, sprintf("%s/%snlms_FuTR2.%s.rds", out_path, ifelse(pilot, "pilot_", ""), i))
      saveRDS(nlms.FuRF1, sprintf("%s/%snlms_FuRF1.%s.rds", out_path, ifelse(pilot, "pilot_", ""), i))
      saveRDS(nlms.FuRF2, sprintf("%s/%snlms_FuRF2.%s.rds", out_path, ifelse(pilot, "pilot_", ""), i))
      #saveRDS(nlms.cox.list, sprintf("%s/%snlms_cox.rds", out_path, ifelse(pilot, "pilot_", ""), i))
      saveRDS(nlms.cox, sprintf("%s/%snlms_cox.rds", out_path, ifelse(pilot, "pilot_", ""), i))
    }
    
    
    
    ## 2. Test set evaluation - Trees
    # prediction for trees
    nlms.pred <- lapply(list(nlms.icrf.H, nlms.icrf.E, nlms.FuTR1, nlms.FuTR2, nlms.FuRF1, nlms.FuRF2),
                        function(s) predict(s, newdata = nlms.test, smooth = T))
    names(nlms.pred) <- c("ICRF.H", "ICRF.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2")
    # test set evaluation for trees
    nlms.eval <- sapply(nlms.pred, function(s) measure(surv.hat = s, tau = tau, 
                                                       timepoints = nlms.icrf.H$time.points.smooth, method = "imse", 
                                                       L = nlms.test$time, 
                                                       R = ifelse(nlms.test$delta, nlms.test$time, Inf))) %>% t
    runtime = c(nlms.icrf.H$runtime$elapsed, nlms.icrf.E$runtime$elapsed, 
                nlms.FuTR1$runtime$elapsed, nlms.FuTR2$runtime$elapsed, 
                nlms.FuRF1$runtime$elapsed, nlms.FuRF2$runtime$elapsed)
    nlms.eval = cbind(nlms.eval, runtime = runtime)
    # print(nlms.eval)
    
    

    ## 3. Cox model test prediction
    nlms.cox.hat.pred <-
      lapply(1:samp2.n, function(i){
        sp_curves <- icenReg::getSCurves(nlms.cox, nlms.test[i, ])
        sp_int    <- t(sp_curves$Tbull_ints[-1, ])
        # When the final interval should be unbounded as long as there is at lest one unbounded interval in data.
        if (is.infinite(max(nlms.train$R))) sp_int[2, dim(sp_int)[2]] = Inf 
        sp_curve  <- sp_curves$S_curves[[1]]
        nn        <- length(sp_curve)
        sp_curve[is.na(sp_curve)] <- 0  # force NaN to be zero (S(near end) = 0)
        npmle <- list(
          intmap = sp_int,                           # First row is redundant.
          pf = sp_curve[-nn] - sp_curve[-1]   # Likewise, first row is redundant
        )
  
        s.hat1 <- 1 - icrf:::isdSm(LR = matrix(c(L = nlms.test$L, R = nlms.test$R), ncol=2),
                                   grid.smooth = nlms.Grid, btt = c(0, 0, tau),
                                   npmle = npmle)
        s.hat2 <- 1 - icrf:::isdSm(LR = matrix(c(L = nlms.test$L, R = nlms.test$R), ncol=2),
                                   grid.smooth = nlms.Grid, btt = c(nlms.icrf.H$method$bandwidth, 0, tau),
                                   npmle = npmle)
        matrix(c(s.hat1, s.hat2), ncol = 2)
      })
    nlms.cox.hat.pred.nonsmooth <- sapply(1:samp2.n, function(s) nlms.cox.hat.pred[[s]][,1]) %>% t
    nlms.cox.hat.pred.smooth <- sapply(1:samp2.n, function(s) nlms.cox.hat.pred[[s]][,2]) %>% t
  
    
    # evaluation
    nlms.cox.imse.pred.nonsmooth <- 
      measure(nlms.cox.hat.pred.nonsmooth, nlms.Grid, tau = tau, 
              method = "imse", L = nlms.test$L, R = nlms.test$R)
    nlms.cox.imse.pred.smooth <- 
      measure(nlms.cox.hat.pred.smooth, nlms.Grid, tau = tau, 
              method = "imse", L = nlms.test$L, R = nlms.test$R)
    # putting cox all in one
    nlms.cox.pred.list <- 
      list(cox = nlms.cox,
           predictedNO = nlms.cox.hat.pred.nonsmooth,
           predictedNO.Sm = nlms.cox.hat.pred.smooth,
           imse.oob = matrix(c(imse.type1 = NaN, imse.type2 = NaN), ncol = 2),
           imse.NO = matrix(nlms.cox.imse.pred.nonsmooth, ncol = 2),
           imse.NO.Sm = matrix(nlms.cox.imse.pred.smooth, ncol = 2),
           test = NULL)
    
    nlms.eval <- 
      rbind(nlms.eval,
            Cox1 = c(nlms.cox.imse.pred.nonsmooth, runtime = NA),
            Cox2 = c(nlms.cox.imse.pred.smooth, runtime = NA))
    print(nlms.eval)
    saveRDS(nlms.eval, paste0(out_path, "/nlms_eval_", i,".rds"))
    gc()
  }
    

  
  if (i == 1) {
  
    # reading back the first replicate
    nlms.icrf.H <- readRDS(paste0(out_path, "/nlms_ICRF.H.rds"))
    nlms.icrf.E <- readRDS(paste0(out_path, "/nlms_ICRF.E.rds"))
    nlms.FuTR1 <- readRDS(paste0(out_path, "/nlms_FuTR1.rds"))
    nlms.FuTR2 <- readRDS(paste0(out_path, "/nlms_FuTR2.rds"))
    nlms.FuRF1 <- readRDS(paste0(out_path, "/nlms_FuRF1.rds"))
    nlms.FuRF2 <- readRDS(paste0(out_path, "/nlms_FuRF2.rds"))
    nlms.cox <- readRDS(paste0(out_path, "/nlms_cox.rds"))
    
    
    ## 4. variable importance   #Health insurance type, age are the most important factors according to IMSE criterion (for node impurity, weight and health is the most important).
    print(nlms.icrf.H$importance %>% apply(c(1,3), mean) %>% as.data.frame %>% arrange(desc(`%IncIMSE1`)))
    print(nlms.icrf.E$importance %>% apply(c(1,3), mean) %>% as.data.frame %>% arrange(desc(`%IncIMSE1`)))
    c(nlms.icrf.H$importance[,nlms.icrf.H$bestFold$bestFold,"%IncIMSE1"], nlms.icrf.E$importance[,nlms.icrf.E$bestFold$bestFold,"%IncIMSE1"]) %>% {./max(.)}
    c(nlms.icrf.H$importance[,nlms.icrf.H$bestFold$bestFold,"%IncIMSE2"], nlms.icrf.E$importance[,nlms.icrf.E$bestFold$bestFold,"%IncIMSE2"]) %>% {./max(.)}
    
    # ## 5. survival prediction
    # data.grid <-
    #   expand.grid(Country = nlms.complete$Country %>% levels, 
    #               GroupActivity = nlms.complete$GroupActivity %>% levels, 
    #               BurialDepth = seq(1, 700, length.out = 100))
    # data.grid.pred.icrf.H <- 
    #   predict(nlms.icrf.H, data.grid, smooth = TRUE)
    # 
    # data.grid.pred.icrf.E <- 
    #   predict(nlms.icrf.E, data.grid, smooth = TRUE)
    # 
    # data.grid.pred.FuRF <- 
    #   predict(nlms.FuRF2, data.grid, smooth = TRUE)
    # 
    # data.grid.pred.FuTR <- 
    #   predict(nlms.FuTR2, data.grid, smooth = TRUE)
    # 
    # all(data.grid.pred.icrf.H[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    # all(data.grid.pred.icrf.E[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    # all(data.grid.pred.FuTR[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    # all(data.grid.pred.FuRF[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    # 
    # data.grid$ICRF.H <-
    #   (as.matrix(data.grid.pred.icrf.H)[, -201] - as.matrix(data.grid.pred.icrf.H)[, -1]) %*% 
    #   nlms.icrf.H$time.points.smooth[-201] + 
    #   as.matrix(data.grid.pred.icrf.H)[, 201] * nlms.icrf.H$time.points.smooth[200] # virtually zero.
    # 
    # data.grid$ICRF.E <-
    #   (as.matrix(data.grid.pred.icrf.E)[, -201] - as.matrix(data.grid.pred.icrf.E)[, -1]) %*% 
    #   nlms.icrf.E$time.points.smooth[-201] + 
    #   as.matrix(data.grid.pred.icrf.E)[, 201] * nlms.icrf.H$time.points.smooth[200] # virtually zero.
    # 
    # data.grid$STIC <-
    #   (as.matrix(data.grid.pred.FuTR)[, -201] - as.matrix(data.grid.pred.FuTR)[, -1]) %*% 
    #   nlms.icrf.H$time.points.smooth[-201] + 
    #   as.matrix(data.grid.pred.FuTR)[, 201] * nlms.icrf.H$time.points.smooth[200] # virtually zero.
    # 
    # 
    # data.grid$SFIC <-
    #   (as.matrix(data.grid.pred.FuRF)[, -201] - as.matrix(data.grid.pred.FuRF)[, -1]) %*% 
    #   nlms.icrf.H$time.points.smooth[-201] + 
    #   as.matrix(data.grid.pred.FuRF)[, 201] * nlms.icrf.H$time.points.smooth[200] # virtually zero.
    # 
    # 
    # data.grid$Cox <-
    #   predict(nlms.cox, data.grid, type = "response")  
    # 
    # data.grid.long <-
    #   data.grid %>% tidyr::gather(key = "method", value = "ET", ICRF.H, ICRF.E, STIC, SFIC, Cox) %>% 
    #   mutate(method = factor(method, levels = c("ICRF.H", "ICRF.E", "STIC", "SFIC", "Cox"),
    #                          labels = c("ICRF (quasi honest)", "ICRF (exploitative)", "STIC (smooth)", "SFIC (Smooth)", "Cox (Smooth)")))
    # 
    # lab.country <- c("Switzerland", "Canada")
    # names(lab.country) <- c("CH", "CND")
    # 
    # lab.activity <- c("Backcountry Skiing", "Out-of-Bounds Skiing", "Mountaineering/Ice Climbing",
    #                   "Mechanized Skiing", "Snowmobiling", "Other Recreational", "Non-Recreational")
    # 
    # data.grid.long %>% 
    #   mutate(GroupActivity = factor(GroupActivity, levels = lab.activity)) %>% 
    #   ggplot(aes(BurialDepth, y = ET, col = GroupActivity)) +
    #   geom_line() +
    #   facet_grid(method ~ Country, labeller = labeller(Country = lab.country)) +
    #   geom_count(data = nlms, mapping = aes(x = BurialDepth, y = 1, col = GroupActivity), 
    #              alpha = 0.5, position = ggstance::position_dodgev(height=0.3)) +
    #   xlab("burial depth (in cm)") +
    #   ylab("expected log survival time (in minutes)") +
    #   theme_bw() +
    #   theme(legend.position = "bottom", legend.box = "vertical") +
    #   guides(size = guide_legend("number of data points"),
    #          color = guide_legend("group activity"))
    # ggsave(paste0(fig_path, "/fig_nlms_ET.png"), width = 20, height = 30, units = "cm")
    # 
    # data.grid.long %>% 
    #   mutate(GroupActivity = factor(GroupActivity, levels = lab.activity)) %>% 
    #   ggplot(aes(BurialDepth, y = ET, col = GroupActivity)) +
    #   geom_line() +
    #   facet_grid(method ~ Country, labeller = labeller(Country = lab.country)) +
    #   geom_count(data = nlms, mapping = aes(x = BurialDepth, y = 1, col = GroupActivity), 
    #              alpha = 0.5, position = ggstance::position_dodgev(height=0.3)) +
    #   xlab("burial depth (in cm)") +
    #   ylab("expected log survival time (in minutes)") +
    #   theme_bw() +
    #   theme(legend.box = "vertical") +
    #   guides(size = guide_legend("number of data points"),
    #          color = guide_legend("group activity"))
    # ggsave(paste0(fig_path, "/fig_nlms_ET_.png"), width = 30, height = 19, units = "cm")
    
  }
  
  
  if (FALSE) {
    library(ggplot2)
    n.eval = 300
    mth = c("ICRF.H", "icrf.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "Cox1", "Cox2")
    nlms.eval.array <-
      array(NA, dim = c(length(mth), 3, n.eval),
            dimnames = list(method  = mth,
                            measure = c("imse.type1", "imse.type2", "runtime"),
                            replicate = 1:n.eval))  # i.
    i.total = 0
    for (i in 1:n.eval) {
      tmp <- readRDS(paste0(out_path, "/nlms_eval_", i, ".rds"))
      if (all(is.na(tmp))) next
      cat(i, " ")
      i.total = i.total + 1
      nlms.eval.array[ , , i] <- readRDS(paste0(out_path, "/nlms_eval_", i, ".rds"))
      # if (i.total >= 100) break
    }
    
    # mean
    nlms.eval.m <- apply(nlms.eval.array, 1:2, mean, na.rm = TRUE) %>% print
    # sd
    nlms.eval.sd <- apply(nlms.eval.array, 1:2, sd, na.rm = TRUE) %>% print

    ## WRS312 plot
    lvs1 <- c("Cox1", "Cox2", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "ICRF.H", "icrf.E")
    lbs1 <- c("Cox", "Cox (smooth)", "STIC", "STIC (smooth)", "SFIC", "SFIC (smooth)", "ICRF (quasi-honest)", "ICRF (exploitative)")
    lvs3 <- c("imse.type1", "imse.type2")
    lbs3 <- c("IMSE1", "IMSE2")
    # lvs5 <- lbs5 <- samp.size[4:1]

    nlms.eval.summary <-
      data.frame(expand.grid(c(dimnames(nlms.eval.m)))) %>%
      mutate(mean = as.vector(nlms.eval.m),
             sd = as.vector(nlms.eval.sd)) %>%
      mutate(method = factor(method, levels = lvs1, labels = lbs1),
             measure = factor(measure, levels = lvs3, labels = lbs3))
    # pd <- position_dodge(0.3)
    # ggplot(nlms.eval.summary, aes(method, mean, col = method, shape = method, group = method)) +
    #   geom_line(position = pd) +
    #   geom_point(position = pd) +
    #   geom_errorbar(aes(ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd), width = 0.1,
    #                 position = pd, alpha = 0.4) +
    #   facet_grid(. ~ measure) +
    #   ylab ("Mean error with 95% confindence intervals") +
    #   xlab ("training sample size")
    # ggsave(paste0(fig_path, "/fig_nlms_size.png"), width = 20, height = 15, units = "cm")
    
    
    nlms.eval.df = 
      nlms.eval.array[, 1,] %>% t %>% as.data.frame %>% # IMSE1 and 2 are equal in this data set.
      tidyr::gather(key = "method", value = "IMSE") %>% 
      mutate(method = factor(method, levels = lvs1, labels = lbs1))
    ggplot(nlms.eval.df, aes(method, IMSE, col = method)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, height = 0, alpha = 0.2) +
      theme(axis.text.x = element_text(angle = 40)) +
      geom_point(data = nlms.eval.summary, aes(method, mean), col = "black", shape = "square")
    ggsave(paste0(fig_path, "/fig_nlms_box.png"), width = 20, height = 15, units = "cm")
  }




# 
# nlms.cox2 <- coxph(Surv(time, delta) ~ age + race + hisp, data = nlms.train)
# nlms.cox2 <- coxph(Surv(time, delta) ~ ., data = nlms.train)
# # Due to singularity (educ == 0 <=> ms == 0), one of the variables (ms) is dropped.
# 
# car::vif(lm(time ~ age + race + hisp + hhnum + reltrf + adjinc + ssnyn + tenure + wt + health + hitype + urban + citizen, 
#             data = nlms.train, singular.ok = TRUE))
# car::vif(lm(time ~ ms + educ, 
#             data = nlms.train, singular.ok = TRUE))
# 
# names(nlms.train[,-(1:2)]) %>% paste(collapse = " + ")
# sapply(nlms.raw, function(x) mean(!is.na(x)))
# 
# table(nlms.raw$urban, nlms.raw$citizen)
# table(nlms.train$pob, nlms.train$citizen)
