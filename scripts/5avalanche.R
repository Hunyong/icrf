# Avalanche data analysis
#install.packages("/Users/hunyongcho/Documents/1Research/Rpackages/icrf_1.0.2.tar.gz", repos = NULL, INSTALL_opts = c('--no-lock'))
library(readxl)
library(dplyr)
library(magrittr)
library(icrf)
require(icenReg)
library(xtable)
library(ggplot2)
library(ggstance)

out_path = "output/avalanche"
if (!dir.exists(out_path)) dir.create(out_path)
fig_path = "figure/avalanche"
if (!dir.exists(fig_path)) dir.create(fig_path)

i = as.numeric(commandArgs(trailingOnly=TRUE))[1]  # passed from script
nfold = as.numeric(commandArgs(trailingOnly=TRUE))[2]
if (is.na(nfold)) nfold = 10
if (file.exists(paste0(out_path, "/avalanche_eval_", i,".rds"))) stop("Already done.")

if (is.na(i)) {
  warning("i (the cross-validation index) was not provided. A full loop will be implemented.")
  i = 1; rng = 1:300; parallel = FALSE
} else {
  rng = i:i; parallel = TRUE
  cat("i = ", i, ", nfold = ", nfold)
}

# 0.0 reading
  avalanche = read_xlsx("Canada_Swiss_Avalanche Data_092112.xlsx", sheet = 2)

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
  
  
  samp.size <- ceiling(aval.complete.n * 0.7)
  ntree=300
  
  for (i in rng) {
    cat("\nreplicate i = ", i, "\n")  
    a <- Sys.time()
    # 0.4 hold-out sample
    set.seed(i)
    # samp1 = train set index
    samp1 <- sample(1:aval.complete.n, samp.size[1]) %>% sort
    lack <- 
      any(! act %in% unique(aval.complete$GroupActivity[samp1])) || 
      any(! cnt %in% unique(aval.complete$Country[samp1]))
    if (lack) {
      print("Not all levels are not present in train set. Go to the next i.")
      # saveRDS(NA, paste0(out_path, "/avalanche_eval_", i,"_NULL.rds"))
      if (parallel) stop("done.") else next
    }
    
    # samp2 = test set index
    samp2 <- which(!1:aval.complete.n %in% samp1)
    samp2.n <- samp2 %>% length # 339
    
    # train and test sets
    aval.train <- aval.complete[samp1, ]
    aval.test  <- aval.complete[samp2, ]
    
    
    # 1.0 models
    
    ## 1.1 ICRF-WRS
    cat("ICRF - Quasi honest\n")
    set.seed(i)
    aval.icrf.H <- 
      icrf:::icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                          L = log(aval.train$L + 1), R = log(aval.train$R + 1),
                          tau = log(tau + 1), proximity = F, importance = (i==1), nPerm = 10,
                          nfold = nfold, ntree = ntree, nodesize = 6, mtry = 2,         # tree structure
                          replace = F, sampsize = samp.size * 0.95,    # resampling 
                          split.rule = "Wilcoxon", ERT = TRUE, uniformERT = TRUE,      # splitting rules
                          quasihonesty = TRUE, timeSmooth = aval.Grid)
    
    cat("ICRF - Exploitative\n")
    set.seed(i)
    aval.icrf.E <- 
      icrf:::icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                          L = log(aval.train$L + 1), R = log(aval.train$R + 1),
                          tau = log(tau + 1), proximity = F, importance = (i==1), nPerm = 10,
                          nfold = nfold, ntree = ntree, nodesize = 6, mtry = 2,         # tree structure
                          replace = F, sampsize = samp.size * 0.95,    # resampling 
                          split.rule = "Wilcoxon", ERT = TRUE, uniformERT = TRUE,      # splitting rules
                          quasihonesty = FALSE, timeSmooth = aval.Grid)
    
    
    bestForest <- data.frame(honest = NA, iter = NA, imse = NA)
    choice <- rbind(aval.icrf.H$bestFold, aval.icrf.E$bestFold)
    bestForest$honest <- choice$imse.best[1] < choice$imse.best[2]
    bestForest$iter <- choice[2 - bestForest$honest, "bestFold"]
    bestForest$imse <- choice[2 - bestForest$honest, "imse.best"]
    
    ## 1.2.1 Fu's Tree 1
    cat("FuTR1\n")
    set.seed(i)
    aval.FuTR1 <- 
      icrf:::icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                          L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                          tau = log(tau + 1), proximity = F, importance = (i==1), nPerm = 10,
                          nfold = 1, ntree = 1, nodesize = 20, mtry = 3,          # tree structure
                          replace = F, sampsize = samp.size,          # resampling 
                          split.rule = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,# splitting rules
                          quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = 0)
    
    ## 1.2.2 Fu's Tree 2
    cat("FuTR2\n")
    set.seed(i)
    aval.FuTR2 <- 
      icrf:::icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                          L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                          tau = log(tau + 1), proximity = F, importance = (i==1), nPerm = 10,
                          nfold = 1, ntree = 1, nodesize = 20, mtry = 3,          # tree structure
                          replace = F, sampsize = samp.size,          # resampling 
                          split.rule = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,# splitting rules
                          quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = NULL)
    
    ## 1.3.1 Fu's RF 1
    cat("FuRF1\n")
    set.seed(i)
    aval.FuRF1 <- 
      icrf:::icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                          L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                          tau = log(tau + 1), proximity = F, importance = (i==1), nPerm = 10,
                          nfold = 1, ntree = ntree, nodesize = 6, mtry = 2,               # tree structure
                          replace = T, sampsize = ceiling(.632 * samp.size),# resampling 
                          split.rule = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,      # splitting rules
                          quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = 0)
    
    ## 1.3.2 Fu's RF 2
    cat("FuRF2\n")
    set.seed(i)
    aval.FuRF2 <- 
      icrf:::icrf.default(x = aval.train[, c("Country", "GroupActivity", "BurialDepth")], 
                          L = log(aval.train$L + 1), R = log(aval.train$R + 1), timeSmooth = aval.Grid,
                          tau = log(tau + 1), proximity = F, importance = (i==1), nPerm = 10,
                          nfold = 1, ntree = ntree, nodesize = 6, mtry = 2,               # tree structure
                          replace = T, sampsize = ceiling(.632 * samp.size),# resampling 
                          split.rule = "PetoLogrank", ERT = FALSE, uniformERT = FALSE,      # splitting rules
                          quasihonesty = TRUE, initialSmoothing = FALSE, bandwidth = NULL)
    
    
    ## 1.4 Cox models
    cat("Cox models\n")
    ## 1.4.1 Cox
    aval.cox <- icenReg::ic_sp(Surv(log(L + 1), log(R + 1), type = 'interval2') ~ Country + GroupActivity + BurialDepth, 
                               model = "ph", data = aval.train)
    
    if (i %in% 1:10) {
      saveRDS(aval.icrf.H, sprintf("%s/%savalanche_ICRF.H.%s.rds", out_path, i))
      saveRDS(aval.icrf.E, sprintf("%s/%savalanche_ICRF.E.%s.rds", out_path, i))
      saveRDS(aval.FuTR1, sprintf("%s/%savalanche_FuTR1.%s.rds", out_path, i))
      saveRDS(aval.FuTR2, sprintf("%s/%savalanche_FuTR2.%s.rds", out_path, i))
      saveRDS(aval.FuRF1, sprintf("%s/%savalanche_FuRF1.%s.rds", out_path, i))
      saveRDS(aval.FuRF2, sprintf("%s/%savalanche_FuTR2.%s.rds", out_path, i))
      #saveRDS(aval.cox.list, paste0(out_path, "/avalanche_cox.rds"))
      saveRDS(aval.cox, sprintf("%s/%savalanche_Cox.%s.rds", out_path, i))
    }
    
    ## 2. Test set evaluation - Trees
    # prediction for trees
    aval.pred <- lapply(list(aval.icrf.H, aval.icrf.E, aval.FuTR1, aval.FuTR2, aval.FuRF1, aval.FuRF2),
                        function(s) predict(s, newdata = aval.test[, c("Country", "GroupActivity", "BurialDepth")]))
    names(aval.pred) <- c("ICRF.H", "ICRF.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2")
    # test set evaluation for trees
    aval.eval <- sapply(aval.pred, function(s) measure(surv.hat = s, tau = log(1 + tau), 
                                                       timepoints = aval.icrf.H$time.points.smooth, method = "imse", 
                                                       L = log(aval.test$L + 1), R = log(aval.test$R + 1))) %>% t
    runtime = c(aval.icrf.H$runtime$elapsed, aval.icrf.E$runtime$elapsed, 
                aval.FuTR1$runtime$elapsed, aval.FuTR2$runtime$elapsed, 
                aval.FuRF1$runtime$elapsed, aval.FuRF2$runtime$elapsed)
    aval.eval = cbind(aval.eval, runtime = runtime)
    
    
    ## 3. Cox model test prediction
    aval.cox.hat.pred <-
      lapply(1:samp2.n, function(i){
        sp_curves <- icenReg::getSCurves(aval.cox, aval.test[i, ])
        sp_int    <- t(sp_curves$Tbull_ints[-1, ])
        # When the final interval should be unbounded as long as there is at lest one unbounded interval in data.
        if (is.infinite(max(aval.train$R))) sp_int[2, dim(sp_int)[2]] = Inf 
        sp_curve  <- sp_curves$S_curves[[1]]
        nn        <- length(sp_curve)
        sp_curve[is.na(sp_curve)] <- 0  # force NaN to be zero (S(near end) = 0)
        npmle <- list(
          intmap = sp_int,                           # First row is redundant.
          pf = sp_curve[-nn] - sp_curve[-1]   # Likewise, first row is redundant
        )

        s.hat1 <- 1 - icrf:::isdSm(LR = matrix(c(L = log(aval.test$L + 1), R = log(aval.test$R + 1)), ncol=2),
                                   grid.smooth = aval.Grid, btt = c(0, 0, tau),
                                   npmle = npmle)
        s.hat2 <- 1 - icrf:::isdSm(LR = matrix(c(L = log(aval.test$L + 1), R = log(aval.test$R + 1)), ncol=2),
                                   grid.smooth = aval.Grid, btt = c(NA, 0, tau),
                                   npmle = npmle)
        matrix(c(s.hat1, s.hat2), ncol = 2)
      })
    aval.cox.hat.pred.nonsmooth <- sapply(1:samp2.n, function(s) aval.cox.hat.pred[[s]][,1]) %>% t
    aval.cox.hat.pred.smooth <- sapply(1:samp2.n, function(s) aval.cox.hat.pred[[s]][,2]) %>% t

    # evaluation
    aval.cox.imse.pred.nonsmooth <- 
      measure(aval.cox.hat.pred.nonsmooth, aval.Grid, tau = log(1 + tau), 
              method = "imse", L = log(aval.test$L + 1), R = log(aval.test$R + 1))
    aval.cox.imse.pred.smooth <- 
      measure(aval.cox.hat.pred.smooth, aval.Grid, tau = log(1 + tau), 
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
    print(aval.eval)
    
    
    # evaluation
    aval.cox.imse.pred.nonsmooth <- 
      measure(aval.cox.hat.pred.nonsmooth, aval.Grid, tau = log(1 + tau), 
              method = "imse", L = log(aval.test$L + 1), R = log(aval.test$R + 1))
    aval.cox.imse.pred.smooth <- 
      measure(aval.cox.hat.pred.smooth, aval.Grid, tau = log(1 + tau), 
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
            Cox1 = c(aval.cox.imse.pred.nonsmooth, runtime = NA),
            Cox2 = c(aval.cox.imse.pred.smooth, runtime = NA))
    print(aval.eval)
  
    print(Sys.time() - a)
    saveRDS(aval.eval, paste0(out_path, "/avalanche_eval_", i,".rds"))
    
    # print(apply(aval.eval.array, 1:3, mean, na.rm = TRUE))
    # print(apply(aval.eval.array, 1:3, sd, na.rm = TRUE))
    print(aval.eval)
  }
  
  if (i == 1) {
    #saveRDS(aval.eval.array, paste0(out_path, "/eval_array_final.rds"))
    #print(aval.eval.mean <- apply(aval.eval.array, 1:3, mean, na.rm = TRUE))
    #print(aval.eval.sd <- apply(aval.eval.array, 1:3, sd, na.rm = TRUE))
    # xtable(aval.eval.mean, digits = 3)
    
    
    
    # reading back the first replicate
    aval.icrf.H <- readRDS(paste0(out_path, "/avalanche_ICRF.H.1.rds"))
    aval.icrf.E <- readRDS(paste0(out_path, "/avalanche_ICRF.E.1.rds"))
    aval.FuTR1 <- readRDS(paste0(out_path, "/avalanche_FuTR1.1.rds"))
    aval.FuTR2 <- readRDS(paste0(out_path, "/avalanche_FuTR2.1.rds"))
    aval.FuRF1 <- readRDS(paste0(out_path, "/avalanche_FuRF1.1.rds"))
    aval.FuRF2 <- readRDS(paste0(out_path, "/avalanche_FuRF2.1.rds"))
    aval.cox <- readRDS(paste0(out_path, "/avalanche_cox.1.rds"))
    
    
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
                             labels = c("ICRF (quasi honest)", "ICRF (exploitative)", "Fu (*)", "Yao (*)", "Cox (*)")))
    
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
      # theme(legend.position = "bottom", legend.box = "vertical") +
      theme(legend.position = "right") +
      guides(size = guide_legend("number of data points"),
             color = guide_legend("group activity"))
    ggsave(paste0(fig_path, "/fig_avalanche_ET.png"), width = 32, height = 18, units = "cm")
    
    # data.grid.long %>% 
    #   mutate(GroupActivity = factor(GroupActivity, levels = lab.activity)) %>% 
    #   ggplot(aes(BurialDepth, y = ET, col = GroupActivity)) +
    #   geom_line() +
    #   facet_grid(method ~ Country, labeller = labeller(Country = lab.country)) +
    #   geom_count(data = avalanche, mapping = aes(x = BurialDepth, y = 1, col = GroupActivity), 
    #              alpha = 0.5, position = ggstance::position_dodgev(height=0.3)) +
    #   xlab("burial depth (in cm)") +
    #   ylab("expected log survival time (in minutes)") +
    #   theme_bw() +
    #   theme(legend.box = "vertical") +
    #   guides(size = guide_legend("number of data points"),
    #          color = guide_legend("group activity"))
    # ggsave(paste0(fig_path, "/fig_avalanche_ET_.png"), width = 30, height = 19, units = "cm")
    
  }
  
  
  if (FALSE) {
    method.nm = c("ICRF.H", "ICRF.E", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "Cox1", "Cox2")
    aval.eval.array <-
      array(NA, dim = c(8, 3, 300),
            dimnames = list(method  = method.nm,
                            measure = c("imse.type1", "imse.type2", "runtime"),
                            replicate = 1:300))  # i.
    
    for (i in 1:300) {
      tmp.nm = paste0(out_path, "/avalanche_eval_", i, ".rds")
      if (file.exists(tmp.nm)) {
        if (file.info(tmp.nm)$size < 1) {
          cat(i, " zero file size.\n ")
          next
        } 
        tmp <- readRDS(tmp.nm)
        if (all(is.na(tmp))) next
        aval.eval.array[ , , i] <- tmp[method.nm, ]
      } else {
        cat(i, " not available.\n ")
        next
      }
      
    }
    # mean
    aval.eval.m <- apply(aval.eval.array, 1:2, mean, na.rm = TRUE)
    # sd
    aval.eval.sd <- apply(aval.eval.array, 1:2, sd, na.rm = TRUE)
    
    ## WRS312 plot
    lvs1 <- c("Cox1", "Cox2", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "ICRF.H", "icrf.E")
    lbs1 <- c("Cox", "Cox (*)", "Fu", "Fu (*)", "Yao", "Yao (*)", "ICRF (quasi-honest)", "ICRF (exploitative)")
    lvs3 <- c("imse.type1", "imse.type2")
    lbs3 <- c("IMSE1", "IMSE2")
    # lvs5 <- lbs5 <- samp.size[4:1]
    
    aval.eval.summary <-   
      data.frame(expand.grid(c(dimnames(aval.eval.m)))) %>% 
      mutate(mean = as.vector(aval.eval.m),
             sd = as.vector(aval.eval.sd)) %>% 
      dplyr::filter(measure != "runtime") %>% 
      mutate(method = factor(method, levels = lvs1, labels = lbs1),
             measure = factor(measure, levels = lvs3, labels = lbs3))
    pd <- position_dodge(0.3)
    ggplot(aval.eval.summary, aes(method, mean, col = method, shape = method, group = method)) +
      geom_line(position = pd) +
      geom_point(position = pd) +
      geom_errorbar(aes(ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd), width = 0.1, 
                    position = pd, alpha = 0.4) + 
      facet_grid(. ~ measure) +
      ylab ("Mean error with 95% confindence intervals") +
      xlab ("training sample size")
    ggsave(paste0(fig_path, "/fig_avalanche_size.png"), width = 20, height = 15, units = "cm")
  }
    
    