# Avalanche data analysis
#install.packages("/Users/hycho/Documents/1Research/Rpackages/icrf_0.1.2.tar.gz", repos = NULL, INSTALL_opts = c('--no-lock'))
library(readxl)
library(dplyr)
library(magrittr)
library(icrf)
require(icenReg)
library(xtable)
library(ggplot2)

out_path = "../output/avalanche"
if (!dir.exists(out_path)) dir.create(out_path)
fig_path = "../figure/avalanche"
if (!dir.exists(fig_path)) dir.create(fig_path)

i = as.numeric(commandArgs(trailingOnly=TRUE))[1]  # passed from script
nfold = as.numeric(commandArgs(trailingOnly=TRUE))[2]
if (is.na(nfold)) nfold = 10
if (file.exists(paste0(out_path, "/avalanche_eval_", i,".rds"))) stop("Already done.")
cat("i = ", i, ", nfold = ", nfold)

# 0.0 reading
  avalanche = read_xlsx("../Canada_Swiss_Avalanche Data_092112.xlsx", sheet = 2)

# 0.1 data cleansing
  
  # 0.1.1 GroupActivity
  # 7 levels: Backcountry Skiing (602) / Out-of-Bounds Skiing (310) / Mountaineering/Ice-Climbing (79)
  #           Mechanized Skiing (69) / Other-Recreational (66) / Snowmobiling (65) / Non-Recreational (56) 
  avalanche$GroupActivity %<>% as.factor
  
  # 0.1.2 CoD
  # 4 levels: Unknown (521) / Asphyxia (116) / Trauma (27) / NA (survivors)
  avalanche$CoD %<>% factor(levels = c("Asphyxia", "Trauma", "Unknown"))
  
  # 0.1.3 Country
  # 2 levels: CH (946) / CND (301)
  avalanche$Country %<>% as.factor
  
  # 0.1.4 BurialDepth
  avalanche$BurialDepth %<>% as.numeric
  sum(is.na(avalanche$BurialDepth)) # 120 out of 1247 are missing
  
  # 0.1.5 survival outcomes into intervals
  avalanche$L = ifelse(avalanche$Survival, avalanche$BurialTime, 0)
  avalanche$R = ifelse(avalanche$Survival, Inf, avalanche$BurialTime)
  
# 0.2 Setting tau = 8760
  # distribution of burial time
  avalanche$BurialTime %>% summary
  avalanche$BurialTime %>% sort %>% plot
  mean(avalanche$BurialTime <= 60 * 24) # 93% of burial times were within 1 day (1440 minutes).
  mean(avalanche$BurialTime <= 60 * 24 * 10) # 98% of burial times were within 10 days (14,400 minutes).
  tau = 14400
  
# 0.3 complete data
  aval.complete <- avalanche[!is.na(avalanche$BurialDepth), ]
  # aval.complete.time <- aval.complete[, c("L", "R")] %>% as.matrix %>% sort %>% unique
  aval.complete.n <- dim(aval.complete)[1]
  aval.Grid <- c(seq(0, log(1 + tau), length.out = 200), Inf)
  
  # check if cox model levels are present in training sets
  cnt <- levels(aval.complete$Country)
  act <- levels(aval.complete$GroupActivity)
  
  
  samp.size <- ceiling(aval.complete.n * 0.7 * 2^-(0:3))
  aval.eval.array.i <- 
    array(NA, dim = c(9, 2, 4), 
          dimnames = list(method  = c("ICRF.best", "ICRF.H", "icrf.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "Cox1", "Cox2"),
                          measure = c("imse.type1", "imse.type2"),
                          size    = samp.size)) # j. = c(789, 395, 198, 99)
  ntree=300
  #ntree=2;nfold=2
  
  # i.total <- 0
  # for (i in 1:300) {
    cat("\nreplicate i = ", i, "\n")  
    a <- Sys.time()
    # 0.4 hold-out sample
    set.seed(i)
    # samp1 = train set index
    samp1 <- sample(1:aval.complete.n, samp.size[1]) %>% sort
    # getting smaller samples as well.
    samp1.list <- lapply(samp.size, function(s) sample(samp1, s) %>% sort)
    min.set <- aval.complete[samp1.list[[4]],]
    lack <- any(! act %in% unique(min.set$GroupActivity)) || any(! cnt %in% unique(min.set$Country))
    if (lack) {
      print("Not all levels are not present in train set. Go to the next i.")
      saveRDS(aval.eval.array.i, paste0(out_path, "/avalanche_eval_", i,".rds"))
      stop("done.")
      # next
    } # else {
      # i.total = i.total + 1
      # if (i.total > 100) break
    #}
    
    # samp2 = test set index
    samp2 <- which(!1:aval.complete.n %in% samp1)
    samp2.n <- samp2 %>% length # 339
    
    # train and test sets
    aval.train.list <- lapply(samp1.list, function(s) aval.complete[s, ])
    aval.test  <- aval.complete[samp2, ]
    
    
    for (j in 1:4) {
      cat("\ni = ", i, ", samp.size[", j, "] = ", samp.size[j], "\n")
      aval.train <- aval.train.list[[j]]
      # 1.0 models
      
      ## 1.1 ICRF-WRS
      cat("ICRF - Quasi honest\n")
      set.seed(i + 100 * j)
      aval.icrf.H <- 
        icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                     L = log(aval.train$L + 1), R = log(aval.train$R + 1),
                     t0 = 0, tau = log(tau + 1), proximity = F, importance = TRUE, nPerm = 10,
                     nfold = nfold, ntree = ntree, nodesize = 6, mtry = 2,         # tree structure
                     replace = F, sampsize = samp.size[j] * 0.95,    # resampling 
                     method = "Wilcoxon", ERT = TRUE, uniformERT = TRUE,      # splitting rules
                     quasihonesty = TRUE, timeSmooth = aval.Grid)
      
      cat("ICRF - Exploitative\n")
      set.seed(i + 100 * j)
      aval.icrf.E <- 
        icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                     L = log(aval.train$L + 1), R = log(aval.train$R + 1),
                     t0 = 0, tau = log(tau + 1), proximity = F, importance = TRUE, nPerm = 10,
                     nfold = nfold, ntree = ntree, nodesize = 6, mtry = 2,         # tree structure
                     replace = F, sampsize = samp.size[j] * 0.95,    # resampling 
                     method = "Wilcoxon", ERT = TRUE, uniformERT = TRUE,      # splitting rules
                     quasihonesty = FALSE, timeSmooth = aval.Grid)
      
      
      bestForest <- data.frame(honest = NA, iter = NA, imse = NA)
      choice <- rbind(aval.icrf.H$bestFold, aval.icrf.E$bestFold)
      bestForest$honest <- choice$imse.best[1] < choice$imse.best[2]
      bestForest$iter <- choice[2 - bestForest$honest, "bestFold"]
      bestForest$imse <- choice[2 - bestForest$honest, "imse.best"]
      
      ## 1.2.1 Fu's Tree 1
      cat("FuTR1\n")
      set.seed(i + 100 * j)
      aval.FuTR1 <- 
        icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                     L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                     t0 = 0, tau = log(tau + 1), proximity = F, importance = TRUE, nPerm = 10,
                     nfold = 1, ntree = 1, nodesize = 20, mtry = 3,          # tree structure
                     replace = F, sampsize = samp.size[j],          # resampling 
                     method = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,# splitting rules
                     quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = 0)
      
      ## 1.2.2 Fu's Tree 2
      cat("FuTR2\n")
      set.seed(i + 100 * j)
      aval.FuTR2 <- 
        icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                     L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                     t0 = 0, tau = log(tau + 1), proximity = F, importance = TRUE, nPerm = 10,
                     nfold = 1, ntree = 1, nodesize = 20, mtry = 3,          # tree structure
                     replace = F, sampsize = samp.size[j],          # resampling 
                     method = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,# splitting rules
                     quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = NULL)
      
      ## 1.3.1 Fu's RF 1
      cat("FuRF1\n")
      set.seed(i + 100 * j)
      aval.FuRF1 <- 
        icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                     L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                     t0 = 0, tau = log(tau + 1), proximity = F, importance = TRUE, nPerm = 10,
                     nfold = 1, ntree = ntree, nodesize = 6, mtry = 2,               # tree structure
                     replace = T, sampsize = ceiling(.632 * samp.size[j]),# resampling 
                     method = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,      # splitting rules
                     quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = 0)
      
      ## 1.3.2 Fu's RF 2
      cat("FuRF2\n")
      set.seed(i + 100 * j)
      aval.FuRF2 <- 
        icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                     L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                     t0 = 0, tau = log(tau + 1), proximity = F, importance = TRUE, nPerm = 10,
                     nfold = 1, ntree = ntree, nodesize = 6, mtry = 2,               # tree structure
                     replace = T, sampsize = ceiling(.632 * samp.size[j]),# resampling 
                     method = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,      # splitting rules
                     quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = NULL)
      
      
      ## 1.4 Cox models
      cat("Cox models\n")
      ## 1.4.1 Cox
      aval.cox <- icenReg::ic_sp(Surv(log(L + 1), log(R + 1), type = 'interval2') ~ Country + GroupActivity + BurialDepth, 
                                 model = "ph", data = aval.train)
      # aval.cox.hat <-
      #   lapply(1:samp.size[j], function(i){
      #     sp_curves <- icenReg::getSCurves(aval.cox, aval.train[i, ])
      #     sp_int    <- t(sp_curves$Tbull_ints[-1, ])
      #     sp_curve  <- sp_curves$S_curves[[1]]
      #     nn        <- length(sp_curve)
      #     sp_curve[is.na(sp_curve)] <- 0  # force NaN to be zero (S(near end) = 0)
      #     npmle <- list(
      #       intmap = sp_int,                           # First row is redundant.
      #       pf = sp_curve[-nn] - sp_curve[-1]   # Likewise, first row is redundant
      #     )
      #     
      #     s.hat1 <- 1 - isdSm(LR = matrix(c(L = log(aval.train$L + 1), R = log(aval.train$R + 1)), ncol=2), 
      #                         grid.smooth = aval.Grid, btt = c(0, 0, tau),
      #                         npmle = npmle)
      #     s.hat2 <- 1 - isdSm(LR = matrix(c(L = log(aval.train$L + 1), R = log(aval.train$R + 1)), ncol=2), 
      #                         grid.smooth = aval.Grid, btt = c(NA, 0, tau),
      #                         npmle = npmle)
      #     matrix(c(s.hat1, s.hat2), ncol = 2)
      #   })
      # 
      # aval.cox.hat.nonsmooth <- sapply(1:samp.size[j], function(s) aval.cox.hat[[s]][,1]) %>% t
      # aval.cox.hat.smooth <- sapply(1:samp.size[j], function(s) aval.cox.hat[[s]][,2]) %>% t
      # 
      # aval.cox.imse.nonsmooth <- 
      #   measure(aval.cox.hat.nonsmooth, aval.Grid, t0 = 0, tau = log(1 + tau), 
      #           method = "imse", L = log(aval.train$L + 1), R = log(aval.train$R + 1))
      # aval.cox.imse.smooth <- 
      #   measure(aval.cox.hat.smooth, aval.Grid, t0 = 0, tau = log(1 + tau), 
      #           method = "imse", L = log(aval.train$L + 1), R = log(aval.train$R + 1))
      # aval.cox.list <- 
      #   list(cox = aval.cox,
      #        predictedNO = aval.cox.hat.nonsmooth,
      #        predictedNO.Sm = aval.cox.hat.smooth,
      #        imse.oob = matrix(c(imse.type1 = NaN, imse.type2 = NaN), ncol = 2),
      #        imse.NO = matrix(aval.cox.imse.nonsmooth, ncol = 2),
      #        imse.NO.Sm = matrix(aval.cox.imse.smooth, ncol = 2),
      #        test = NULL)
      
      if (i == 1 && j == 1) {
        saveRDS(aval.icrf.H, paste0(out_path, "/avalanche_ICRF.H.rds"))
        saveRDS(aval.icrf.E, paste0(out_path, "/avalanche_ICRF.E.rds"))
        saveRDS(aval.FuTR1, paste0(out_path, "/avalanche_FuTR1.rds"))
        saveRDS(aval.FuTR2, paste0(out_path, "/avalanche_FuTR2.rds"))
        saveRDS(aval.FuRF1, paste0(out_path, "/avalanche_FuRF1.rds"))
        saveRDS(aval.FuRF2, paste0(out_path, "/avalanche_FuRF2.rds"))
        #saveRDS(aval.cox.list, paste0(out_path, "/avalanche_cox.rds"))
        saveRDS(aval.cox, paste0(out_path, "/avalanche_cox.rds"))
      }
      
      ## 2. Test set evaluation - Trees
      # prediction for trees
      aval.pred <- lapply(list(aval.icrf.H, aval.icrf.E, aval.FuTR1, aval.FuTR2, aval.FuRF1, aval.FuRF2),
                          function(s) predict(s, newdata = aval.test[, c("Country", "GroupActivity", "BurialDepth")]))
      names(aval.pred) <- c("ICRF.H", "ICRF.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2")
      # test set evaluation for trees
      aval.eval <- sapply(aval.pred, function(s) measure(surv.hat = s, t0 = 0, tau = log(1 + tau), 
                                                         timepoints = aval.icrf.H$time.points.smooth, method = "imse", 
                                                         L = log(aval.test$L + 1), R = log(aval.test$R + 1))) %>% t
      
      
      ## 3. Cox model test prediction
      aval.cox.hat.pred <-
        lapply(1:samp2.n, function(i){
          sp_curves <- icenReg::getSCurves(aval.cox, aval.test[i, ])
          sp_int    <- t(sp_curves$Tbull_ints[-1, ])
          sp_curve  <- sp_curves$S_curves[[1]]
          nn        <- length(sp_curve)
          sp_curve[is.na(sp_curve)] <- 0  # force NaN to be zero (S(near end) = 0)
          npmle <- list(
            intmap = sp_int,                           # First row is redundant.
            pf = sp_curve[-nn] - sp_curve[-1]   # Likewise, first row is redundant
          )

          s.hat1 <- 1 - isdSm(LR = matrix(c(L = log(aval.test$L + 1), R = log(aval.test$R + 1)), ncol=2),
                              grid.smooth = aval.Grid, btt = c(0, 0, tau),
                              npmle = npmle)
          s.hat2 <- 1 - isdSm(LR = matrix(c(L = log(aval.test$L + 1), R = log(aval.test$R + 1)), ncol=2),
                              grid.smooth = aval.Grid, btt = c(NA, 0, tau),
                              npmle = npmle)
          matrix(c(s.hat1, s.hat2), ncol = 2)
        })
      aval.cox.hat.pred.nonsmooth <- sapply(1:samp2.n, function(s) aval.cox.hat.pred[[s]][,1]) %>% t
      aval.cox.hat.pred.smooth <- sapply(1:samp2.n, function(s) aval.cox.hat.pred[[s]][,2]) %>% t

      # evaluation
      aval.cox.imse.pred.nonsmooth <- 
        measure(aval.cox.hat.pred.nonsmooth, aval.Grid, t0 = 0, tau = log(1 + tau), 
                method = "imse", L = log(aval.test$L + 1), R = log(aval.test$R + 1))
      aval.cox.imse.pred.smooth <- 
        measure(aval.cox.hat.pred.smooth, aval.Grid, t0 = 0, tau = log(1 + tau), 
                method = "imse", L = log(aval.test$L + 1), R = log(aval.test$R + 1))
      # putting cox all in one
      aval.cox.pred.list <- 
        list(cox = aval.cox,
             predictedNO = aval.cox.hat.pred.nonsmooth,
             predictedNO.Sm = aval.cox.hat.pred.smooth,
             imse.oob = matrix(c(imse.type1 = NaN, imse.type2 = NaN), ncol = 2),
             imse.NO = matrix(aval.cox.imse.pred.nonsmooth, ncol = 2),
             imse.NO.Sm = matrix(aval.cox.imse.pred.smooth, ncol = 2),
             test = NULL)
      
      aval.eval <- 
        rbind(ICRF.best = aval.eval[2 - bestForest$honest, ],
              aval.eval,
              Cox1 = aval.cox.imse.pred.nonsmooth,
              Cox2 = aval.cox.imse.pred.smooth)
      #aval.eval.array[ , , j, i] <- aval.eval
      aval.eval.array.i[ , , j] <- aval.eval
      print(aval.eval)
    }
    print(Sys.time() - a)
    saveRDS(aval.eval.array.i, paste0(out_path, "/avalanche_eval_", i,".rds"))
    
    # print(apply(aval.eval.array, 1:3, mean, na.rm = TRUE))
    # print(apply(aval.eval.array, 1:3, sd, na.rm = TRUE))
    print(aval.eval.array.i)
  #}
  
  if (i == 1) {
    #saveRDS(aval.eval.array, paste0(out_path, "/eval_array_final.rds"))
    #print(aval.eval.mean <- apply(aval.eval.array, 1:3, mean, na.rm = TRUE))
    #print(aval.eval.sd <- apply(aval.eval.array, 1:3, sd, na.rm = TRUE))
    # xtable(aval.eval.mean, digits = 3)
    
    
    
    # reading back the first replicate
    aval.icrf.H <- readRDS(paste0(out_path, "/avalanche_ICRF.H.rds"))
    aval.icrf.E <- readRDS(paste0(out_path, "/avalanche_ICRF.E.rds"))
    aval.FuTR1 <- readRDS(paste0(out_path, "/avalanche_FuTR1.rds"))
    aval.FuTR2 <- readRDS(paste0(out_path, "/avalanche_FuTR2.rds"))
    aval.FuRF1 <- readRDS(paste0(out_path, "/avalanche_FuRF1.rds"))
    aval.FuRF2 <- readRDS(paste0(out_path, "/avalanche_FuRF2.rds"))
    aval.cox <- readRDS(paste0(out_path, "/avalanche_cox.rds"))
    
    
    ## 4. variable importance
    print(aval.icrf.H$importance)
    print(aval.icrf.E$importance)
    c(aval.icrf.H$importance[,10,"%IncIMSE1"], aval.icrf.E$importance[,10,"%IncIMSE1"]) %>% {./max(.)}
    c(aval.icrf.H$importance[,10,"%IncIMSE2"], aval.icrf.E$importance[,10,"%IncIMSE2"]) %>% {./max(.)}
    
    ## 5. survival prediction
    data.grid <-
      expand.grid(Country = aval.complete$Country %>% levels, 
                  GroupActivity = aval.complete$GroupActivity %>% levels, 
                  BurialDepth = seq(1, 700, length.out = 100))
    data.grid.pred.icrf.H <- 
      predict(aval.icrf.H, data.grid, smooth = TRUE)
    
    data.grid.pred.icrf.E <- 
      predict(aval.icrf.E, data.grid, smooth = TRUE)
    
    data.grid.pred.FuRF <- 
      predict(aval.FuRF2, data.grid, smooth = TRUE)
    
    data.grid.pred.FuTR <- 
      predict(aval.FuTR2, data.grid, smooth = TRUE)
    
    all(data.grid.pred.icrf.H[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    all(data.grid.pred.icrf.E[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    all(data.grid.pred.FuTR[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    all(data.grid.pred.FuRF[, 201] < 0.01) # survival probability after tau is all cloase to zero.
    
    data.grid$ICRF.H <-
      (as.matrix(data.grid.pred.icrf.H)[, -201] - as.matrix(data.grid.pred.icrf.H)[, -1]) %*% 
      aval.icrf.H$time.points.smooth[-201] + 
      as.matrix(data.grid.pred.icrf.H)[, 201] * aval.icrf.H$time.points.smooth[200] # virtually zero.
    
    data.grid$ICRF.E <-
      (as.matrix(data.grid.pred.icrf.E)[, -201] - as.matrix(data.grid.pred.icrf.E)[, -1]) %*% 
      aval.icrf.E$time.points.smooth[-201] + 
      as.matrix(data.grid.pred.icrf.E)[, 201] * aval.icrf.H$time.points.smooth[200] # virtually zero.
    
    data.grid$STIC <-
      (as.matrix(data.grid.pred.FuTR)[, -201] - as.matrix(data.grid.pred.FuTR)[, -1]) %*% 
      aval.icrf.H$time.points.smooth[-201] + 
      as.matrix(data.grid.pred.FuTR)[, 201] * aval.icrf.H$time.points.smooth[200] # virtually zero.
    
    
    data.grid$SFIC <-
      (as.matrix(data.grid.pred.FuRF)[, -201] - as.matrix(data.grid.pred.FuRF)[, -1]) %*% 
      aval.icrf.H$time.points.smooth[-201] + 
      as.matrix(data.grid.pred.FuRF)[, 201] * aval.icrf.H$time.points.smooth[200] # virtually zero.
    
    
    data.grid$Cox <-
      predict(aval.cox, data.grid, type = "response")  
    
    data.grid.long <-
      data.grid %>% tidyr::gather(key = "method", value = "ET", ICRF.H, ICRF.E, STIC, SFIC, Cox) %>% 
      mutate(method = factor(method, levels = c("ICRF.H", "ICRF.E", "STIC", "SFIC", "Cox"),
                             labels = c("ICRF (quasi honest)", "ICRF (exploitative)", "STIC (smooth)", "SFIC (Smooth)", "Cox (Smooth)")))
    
    lab.country <- c("Switzerland", "Canada")
    names(lab.country) <- c("CH", "CND")
    
    lab.activity <- c("Backcountry Skiing", "Out-of-Bounds Skiing", "Mountaineering/Ice Climbing",
                      "Mechanized Skiing", "Snowmobiling", "Other Recreational", "Non-Recreational")
    
    data.grid.long %>% 
      mutate(GroupActivity = factor(GroupActivity, levels = lab.activity)) %>% 
      ggplot(aes(BurialDepth, y = ET, col = GroupActivity)) +
      geom_line() +
      facet_grid(method ~ Country, labeller = labeller(Country = lab.country)) +
      geom_count(data = avalanche, mapping = aes(x = BurialDepth, y = 1, col = GroupActivity), 
                 alpha = 0.5, position = ggstance::position_dodgev(height=0.3)) +
      xlab("burial depth (in cm)") +
      ylab("expected log survival time (in minutes)") +
      theme_bw() +
      theme(legend.position = "bottom", legend.box = "vertical") +
      guides(size = guide_legend("number of data points"),
             color = guide_legend("group activity"))
    ggsave(paste0(fig_path, "/fig_avalanche_ET.png"), width = 20, height = 30, units = "cm")
    
    data.grid.long %>% 
      mutate(GroupActivity = factor(GroupActivity, levels = lab.activity)) %>% 
      ggplot(aes(BurialDepth, y = ET, col = GroupActivity)) +
      geom_line() +
      facet_grid(method ~ Country, labeller = labeller(Country = lab.country)) +
      geom_count(data = avalanche, mapping = aes(x = BurialDepth, y = 1, col = GroupActivity), 
                 alpha = 0.5, position = ggstance::position_dodgev(height=0.3)) +
      xlab("burial depth (in cm)") +
      ylab("expected log survival time (in minutes)") +
      theme_bw() +
      theme(legend.box = "vertical") +
      guides(size = guide_legend("number of data points"),
             color = guide_legend("group activity"))
    ggsave(paste0(fig_path, "/fig_avalanche_ET_.png"), width = 30, height = 19, units = "cm")
    
    # ## 5.2 IMSE
    #   aval.test$imse1.icrf <- 
    #     sapply(1:samp2.n, function(i)
    #       measure(aval.pred$ICRF[i, , drop = FALSE], timepoints = aval.icrf.H$time.points.smooth, 
    #               tau = log(1 + tau), method = "imse", 
    #               L = log(1 + aval.test$L[i]), R = log(1 + aval.test$R[i]))["imse.type1"])
    #   aval.test$imse1.cox <- 
    #     sapply(1:samp2.n, function(i)
    #       measure(aval.cox.hat.pred.smooth[i, , drop = FALSE], timepoints = aval.icrf.H$time.points.smooth, 
    #               tau = log(1 + tau), method = "imse", 
    #               L = log(1 + aval.test$L[i]), R = log(1 + aval.test$R[i]))["imse.type1"])
    #   aval.test.long <-
    #     aval.test %>% tidyr::gather(key = "method", value = "IMSE1", imse1.icrf, imse1.cox) %>% 
    #     mutate(method = factor(method, levels = c("imse1.icrf", "imse1.cox"), labels = c("ICRF", "Cox")))
    #   
    #   aval.test.long %>% 
    #     ggplot(aes(BurialDepth, y = IMSE1, col = GroupActivity)) +
    #     facet_grid(method ~ Country, labeller = labeller(Country = lab.country)) +
    #     geom_smooth(method = "loess", se = FALSE, fullrange = TRUE) +
    #     geom_point(aes(x = BurialDepth, y = IMSE1, col = GroupActivity), alpha = 0.5) +
    #     xlab("burial depth") +
    #     ylab("expected survival time")
    #    
  }
  
  
  if (FALSE) {
    aval.eval.array <-
      array(NA, dim = c(8, 2, 4, 300),
            dimnames = list(method  = c("ICRF.H", "icrf.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "Cox1", "Cox2"),
                            measure = c("imse.type1", "imse.type2"),
                            size    = samp.size, # j. = c(789, 395, 198, 99)
                            replicate = 1:300))  # i.
    i.total = 0
    for (i in 1:120) {
      tmp <- readRDS(paste0(out_path, "/avalanche_eval_", i, ".rds"))
      if (all(is.na(tmp))) next
      cat(i, " ")
      i.total = i.total + 1
      aval.eval.array[ , , , i] <- readRDS(paste0(out_path, "/avalanche_eval_", i, ".rds"))
      if (i.total >= 100) break
    }
    # mean
    aval.eval.m <- apply(aval.eval.array, 1:3, mean, na.rm = TRUE)
    # sd
    aval.eval.sd <- apply(aval.eval.array, 1:3, sd, na.rm = TRUE)
    
    ## WRS312 plot
    lvs1 <- c("Cox1", "Cox2", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "ICRF.H", "icrf.E")
    lbs1 <- c("Cox", "Cox (smooth)", "STIC", "STIC (smooth)", "SFIC", "SFIC (smooth)", "ICRF (quasi-honest)", "ICRF (exploitative)")
    lvs3 <- c("imse.type1", "imse.type2")
    lbs3 <- c("IMSE1", "IMSE2")
    lvs5 <- lbs5 <- samp.size[4:1]
    
    aval.eval.summary <-   
      data.frame(expand.grid(c(dimnames(aval.eval.m)))) %>% 
      mutate(mean = as.vector(aval.eval.m),
             sd = as.vector(aval.eval.sd)) %>% 
      mutate(method = factor(method, levels = lvs1, labels = lbs1),
             measure = factor(measure, levels = lvs3, labels = lbs3),
             size = factor(size, levels = lvs5, labels = lbs5))
    pd <- position_dodge(0.3)
    ggplot(aval.eval.summary, aes(size, mean, col = method, shape = method, group = method)) +
      geom_line(position = pd) +
      geom_point(position = pd) +
      geom_errorbar(aes(ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd), width = 0.1, 
                    position = pd, alpha = 0.4) + 
      facet_grid(. ~ measure) +
      ylab ("Mean error with 95% confindence intervals") +
      xlab ("training sample size")
    ggsave(paste0(fig_path, "/fig_avalanche_size.png"), width = 20, height = 15, units = "cm")
  }
    
    