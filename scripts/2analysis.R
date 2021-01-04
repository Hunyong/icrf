date = "2019-09-13"
ntree = 300; n.sim=100

library(ggplot2); library(dplyr); library(purrr)
library(icenReg); library(MASS)
source("0functions.R")
source("2analysis-functions.R")
path_output   <- paste0("../output/", date, "/")
path_figure    <- paste0("../figure/", date, "/")
if (!dir.exists(path_output)) {
  dir.create(path_output)
  warning("No such folder. Created one.")
}
if (!dir.exists(path_figure)) {
  dir.create(path_figure)
  warning("No such folder. Created one.")
}

fn_eval_tmp   <- paste0(path_output, "eval_scenario_", "scn", "-n.m_", "nm",
                        "-nT_", ntree, "-rep_", "sm", ".rds")
fn_pred_tmp   <- paste0(path_output, "sim_scenario_", "scn", "-n.m_", "nm",
                        "-nT_", ntree, "-rep_", "sm", ".rds")

  tmp <- fnfun(1, 1, 1, fn = fn_eval_tmp)
  tmp[,,]<- NA
  dm <- dim(tmp)
  # tmp <- fnfun(1000,12,12)  #NULL example
  dmn <- dimnames(tmp)
  #dmn[[3]] <- c(dmn[[3]], "w135")
  # Column vectors
  grpVec <- do.call(expand.grid, dimnames(tmp))
  names(grpVec) <- c("fold", "measure", "method")
  # grpVec2 <- do.call(expand.grid, dmn)
  # names(grpVec2) <- c("fold", "measure", "method")
  lvs2 <- c(paste0(c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2"), "-1"), 
            paste0("ICRF", "-", c(1:10, "A", "B", "C")))
  lbs2 <- c("Cox", "Cox (smooth)", "STIC", "STIC (smooth)", "SFIC", "SFIC (smooth)", 
            paste0("ICRF-", c(1:10, "A", "B", "C")))
  pd <- position_dodge(0.2)
  
for (n.monitor in c(1, 3)) {
  fn_fig = paste0(path_figure, "fig_split_n_m_", n.monitor,".png")
  for (scenario in 1:6) {
print(scenario)  
    #lapply(1:n.sim, function(s) {
    result <- 
      lapply(1:n.sim, function(s) {
cat(s, " ")
        tmp1 <- fnfun(scenario, n.monitor, sim = s, fn = fn_eval_tmp)
        if (is.null(tmp1)) stop("null")
        # tmp2 <- fnfun(scenario, n.monitor, sim = s, fn = gsub("eval", "NonHonestyEval", fn_eval_tmp))
        # if (is.null(tmp1)| is.null(tmp2)) stop("null")
        if (!identical(dim(tmp1), dm)) {
          tmp1 <- tmp1[,,dimnames(tmp)[[3]]]     # in case we only ran some of the methods.
          if (!identical(dim(tmp1), dm)) stop("wrong dimension")
        }
        tmp1 <- as.vector(tmp1)
        #tmp1 <- c(as.vector(tmp1), as.vector(tmp2))
        tmp1
      }) %>% do.call(c, .)
    
    # result.best <- # best among quasihonest and exploitative 10 folds
    #   data.frame(value = result, grpVec, sim = rep(1:n.sim, each = dim(grpVec)[1])) %>%
    #   tidyr::spread(key = measure, value = value) %>%
    #   mutate(method2 = gsub("(132|135)", "", method)) %>% # in order to group by the same methods
    #   group_by(sim, method2) %>%
    #   dplyr::filter(`imse.type2 (oob)`== min(`imse.type2 (oob)`)) %>% # find the best imse1.
    #   mutate(method = gsub("(132|135)", 130, method), fold = "A", value = int.error, measure = "int.error") %>%
    #   ungroup %>%
    #   dplyr::select(value, fold, measure, method, sim)
    
    result.A <- # best for each of quasihonest and exploitative
      data.frame(value = result, grpVec, sim = rep(1:n.sim, each = dim(grpVec)[1])) %>%
      tidyr::spread(key = measure, value = value) %>%
      dplyr::filter(!is.na(`ibs.type1 (oob)`)) %>%
      group_by(sim, method) %>%
      dplyr::filter(`ibs.type1 (oob)`== min(`ibs.type1 (oob)`)) %>% # find the best imse1.
      mutate(fold = "A", value = int.error, measure = "int.error") %>%
      ungroup %>%
      dplyr::select(value, fold, measure, method, sim)
    
    result <- 
      data.frame(value = result, grpVec, sim = rep(1:n.sim, each = dim(grpVec)[1])) %>% 
    # result <- data.frame(value = result, grpVec2, sim = rep(1:n.sim, each = dim(grpVec2)[1])) %>% 
      na.omit %>% 
      dplyr::filter(grepl("int\\.error", measure)) %>% 
      rbind(result.A) %>% 
      # rbind(result.best) %>% 
      mutate(methods = paste0(method, "-", fold))
    result.mean <- aggregate(value ~ methods + measure + method, data = result, mean, na.rm = T)
    result.min <- aggregate(value ~ measure, data = result.mean, min, na.rm = TRUE)
    names(result.min)[2] = "min"
    result.mean <- left_join(result.mean, result.min, by = "measure")
    
    result$scenario = scenario
    result.mean$scenario = scenario
    result.min$scenario = scenario
    
    if (scenario == 1) {
      grandResult <- result
      grandResult.mean <- result.mean
      grandResult.min <- result.min
    } else {
      grandResult <- rbind(grandResult, result)
      grandResult.mean <- rbind(grandResult.mean, result.mean)
      grandResult.min <- rbind(grandResult.min, result.min)
    }
  }
  grandResult <-
    grandResult %>% 
    mutate(scenario = factor(scenario, levels = lvs4, labels = lbs4))
  grandResult.mean <-
    grandResult.mean %>% 
    mutate(scenario = factor(scenario, levels = lvs4, labels = lbs4))
  grandResult.min <-
    grandResult.min %>% 
    mutate(scenario = factor(scenario, levels = lvs4, labels = lbs4))
  
  
  grandResult.gg <-
    grandResult %>% 
    #dplyr::filter(n.monitor == 1) %>% 
    dplyr::filter(fold %in% c(1, 2, 3, 5, 10, "A")) %>%
    dplyr::filter(grepl("(p|w|l)", method)) %>%
    dplyr::filter(grepl("int\\.error", measure)) %>%
    mutate(honesty = ifelse(grepl("132", method), "quasi-honest", 
                            ifelse(grepl("135", method), "exploitative",
                                   "best IMSE1"))) %>% 
    mutate(split = ifelse(grepl("^w", method), "GWRS", 
                          ifelse(grepl("^l", method), "GLR",
                                 ifelse(grepl("^pw", method), "SWRS", "SLR")))) %>% 
    mutate(methods = paste0("ICRF-", fold)) %>% 
    mutate(methods = factor(methods, levels = lvs2, labels = lbs2 )) %>% 
    group_by(fold, measure, split, methods, honesty, scenario) %>% 
    summarize(mean  = mean(value, na.rm = TRUE), 
              Q1 = quantile(value, 0.25, na.rm = TRUE), 
              Q3 = quantile(value, 0.75, na.rm = TRUE)) %>% 
    ungroup 
    
    # dplyr::filter(grepl("(cox|w)", method)) %>% 
    # dplyr::filter(grepl("int", measure)) %>% 
    #dplyr::filter(grepl("(p|w)", method) | grepl("Fu", method)| grepl("cox", method)) %>% 
    # dplyr::filter(!grepl("oob", measure)) %>% 
  grandResult.gg %>%   
    ggplot(aes(x = methods, y = mean, col = honesty, group = honesty)) +
    #geom_line(alpha = 0.3) +
    geom_point(aes(shape=honesty)) +
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.1, alpha = 0.3, position = pd) + 
    # stat_summary(fun.y=mean, geom="point", shape = 15, size = .5) +
    # geom_boxplot(outlier.shape = NA) + 
    # geom_jitter(size = 0.1, width = 0.3, height = 0, size= 0.1, alpha = 0.3) +
    # stat_summary(fun.y=meanPlus1, geom="point", shape = 15, size = .5) +
    #position=position_dodge(width=0.75)) +
    # stat_summary(fun.data = meanFunction, geom="text", size = 4, vjust=1.3) +
    geom_hline (data = grandResult.mean, aes(yintercept = min), color = "black") +
    facet_grid(scenario ~ split, scale = "free") + guides(label = FALSE) +
    geom_line(data = grandResult.gg %>% filter(fold != "A"), alpha = 0.3) +
    # geom_text(data = result.mean, size = 2, angle = 90, mapping = aes(label = round(value, 3), x = methods, y = 1)) +
    theme_bw() +  xlab("iterations") +
    ylab("integrated error (1Q, mean, 3Q)") + 
    theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") # +
    # ggtitle(paste0("scenario ", scenario, ", with n.monitor =", n.monitor))
  ggsave(fn_fig,  width = 6, height = 5)
  gc()
  
}
  
  
  # result.icrf.conv <-
  #   result %>%
  #   dplyr::filter(method == "w135", measure == "imse.type1 (oob)") %>%
  #   #dplyr::filter(method == "w132", measure == "imse.type1 (oob)") %>%
  #   group_by(method, sim) %>%
  #   summarize(convA = convMonitor(value, type = "global"),
  #             convB = convMonitor(value, type = "local"),
  #             convC = convMonitor(value, type = "first")) %>% 
  #   left_join(result %>% dplyr::filter(method == "w135"), 
  #             #left_join(result %>% dplyr::filter(method == "w132"), 
  #             ., 
  #             by = c("method", "sim"))
  # result.A <- 
  #   result.icrf.conv %>% 
  #   dplyr::filter(fold == convA) %>% 
  #   mutate(fold = "A", methods = paste0(method, "-A")) %>% 
  #   dplyr::select(-convA, - convB, - convC)
  # 
  # result.B <- 
  #   result.icrf.conv %>% 
  #   dplyr::filter(fold == convB) %>% 
  #   mutate(fold = "B", methods = paste0(method, "-B")) %>% 
  #   dplyr::select(-convA, - convB, - convC)
  # 
  # result.C <- 
  #   result.icrf.conv %>% 
  #   dplyr::filter(fold == convC) %>% 
  #   mutate(fold = "C", methods = paste0(method, "-C")) %>% 
  #   dplyr::select(-convA, - convB, - convC)
  # 
  # result <-
  #   rbind(result, result.A, result.B, result.C)
  # 
  
## WRS312 plot
type = "w132"
wcoxFu = paste0("(^", type, "|cox|Fu)")
lvs1 <- c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2", type)
lbs1 <- c("Cox", "Cox (smooth)", "STIC", "STIC (smooth)", "SFIC", "SFIC (smooth)", "ICRF")
lvs2 <- c(paste0(c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2"), "-1"), paste0(type, "-", c(1:10, "A", "B", "C")))
lbs2 <- c("Cox", "Cox (smooth)", "STIC", "STIC (smooth)", "SFIC", "SFIC (smooth)", paste0("ICRF-", c(1:10, "A", "B", "C")))
lvs3 <- c("int.error", "sup.error", "ibs.type1 (oob)", "ibs.type2 (oob)", "imse.type1", "imse.type2")
lbs3 <- c("integrated error", "supremum error", "IMSE1 (out of bag)", "IMSE2 (out of bag)", "IMSE1", "IMSE2")
lvs4 <- 1:6
lbs4 <- paste0("scenario ", lvs4)
p <- list()
for (n.monitor in c(1, 3)) {
  fn_fig3 = paste0(path_figure, "figWRS_total_", type, "_n_m_", n.monitor,".png")
  for (scenario in 1:6) {
    cat("n.monitor: ", n.monitor, ", scenario = ", scenario, "\n")
    #lapply(1:n.sim, function(s) {
    result <- 
      lapply(1:n.sim, function(s) {
        tmp1 <- fnfun(scenario, n.monitor, sim = s, fn = fn_eval_tmp)
        if (!is.null(tmp1) && !identical(dim(tmp1), dm)) {
          tmp1 <- tmp1[,,dimnames(tmp)[[3]]]     # in case we only ran some of the methods.
          if (!identical(dim(tmp1), dm)) stop("wrong dimension")
        } else if (is.null(tmp1)) {
          warning(paste0(s, " is not available."))
          tmp1 <- tmp
        }
        tmp1 <- as.vector(tmp1)
      }) %>% do.call(c, .)
    
    result <- data.frame(value = result, grpVec, sim = rep(1:n.sim, each = dim(grpVec)[1])) %>% 
      na.omit %>% 
      dplyr::filter(grepl("(type1 \\(oob\\)|error)", measure)) %>% 
      dplyr::filter(grepl(wcoxFu, method)) %>%
      mutate(methods = paste0(method, "-", fold))
    
    
    result.icrf.conv <-
      result %>%
      dplyr::filter(method == type, measure == "ibs.type1 (oob)") %>%
      group_by(method, sim) %>%
      summarize(convA = convMonitor(value, type = "global"),
                convB = convMonitor(value, type = "local"),
                convC = convMonitor(value, type = "first")) %>% 
      left_join(result %>% dplyr::filter(method == type), 
                ., 
                by = c("method", "sim"))
    result.A <- 
      result.icrf.conv %>% 
      dplyr::filter(fold == convA) %>% 
      mutate(fold = "A", methods = paste0(method, "-A")) %>% 
      dplyr::select(-convA, - convB, -convC)
    
    result.B <- 
      result.icrf.conv %>% 
      dplyr::filter(fold == convB) %>% 
      mutate(fold = "B", methods = paste0(method, "-B")) %>% 
      dplyr::select(-convA, - convB, -convC)
    
    result.C <- 
      result.icrf.conv %>% 
      dplyr::filter(fold == convC) %>% 
      mutate(fold = "C", methods = paste0(method, "-C")) %>% 
      dplyr::select(-convA, - convB, -convC)
    
    result <-
      rbind(result, result.A, result.B, result.C)
    
    result$method = factor(result$method, levels = lvs1, labels = lbs1)
    result$methods = factor(result$methods, levels = lvs2, labels = lbs2)
    result$measure = factor(result$measure, levels = lvs3, labels = lbs3)
    result.mean <- aggregate(value ~ methods + measure + method, data = result, mean, na.rm = T)
    result.min <- aggregate(value ~ measure, data = result.mean, min, na.rm = TRUE)
    names(result.min)[2] = "min"
    result.mean <- left_join(result.mean, result.min, by = "measure")
    result$scenario = scenario
    result.mean$scenario = scenario
    result.min$scenario = scenario
    
    if (scenario == 1) {
      grandResult <- result
      grandResult.mean <- result.mean
      grandResult.min <- result.min
    } else {
      grandResult <- rbind(grandResult, result)
      grandResult.mean <- rbind(grandResult.mean, result.mean)
      grandResult.min <- rbind(grandResult.min, result.min)
    }
  }
  grandResult$scenario = factor(grandResult$scenario, levels = lvs4, labels = lbs4)
  grandResult.mean$scenario = factor(grandResult.mean$scenario, levels = lvs4, labels = lbs4)
  grandResult.min$scenario = factor(grandResult.min$scenario, levels = lvs4, labels = lbs4)
  imse.multiplier = 4
  grandResult %>% 
    #dplyr::filter(fold %in% c(1:3, 5, 10, "A", "B", "C")) %>% 
    #dplyr::filter(fold %in% c(1:3, 5, 10)) %>% 
    dplyr::filter(fold %in% c(1:3, 5, 10, "A")) %>% 
    mutate(value = ifelse(grepl("IMS", measure), value * imse.multiplier, value)) %>% 
    ggplot(aes(x = methods, y = value, col = method)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(size = 0.1, width = 0.3, height = 0, size= 0.1, alpha = 0.3) +
    stat_summary(fun.y=mean, geom="point", shape = 22, size = 0.5, col = 'deeppink4', fill = "white") +
    #position=position_dodge(width=0.75)) +
    # stat_summary(fun.data = meanFunction, geom="text", size = 4, vjust=1.3) +
    geom_hline (data = grandResult.mean %>% 
                  mutate(value = ifelse(grepl("IMSE", measure), value * imse.multiplier, value)) %>% 
                  mutate(min = ifelse(grepl("IMSE", measure), min * imse.multiplier, min)), 
                aes(yintercept = min), color = "black") +
    scale_y_continuous("integrated and supremum error", 
                       sec.axis = sec_axis(~ . / imse.multiplier, name = "IMSE")) +
    facet_grid(scenario ~ measure, scales = "free") + guides(label = FALSE, col = FALSE) +
    # geom_text(data = result.mean, size = 2, angle = 90, mapping = aes(label = round(value, 3), x = methods, y = 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) -> p.tmp
    # ggtitle(paste0(# monitoring times =", n.monitor))
  ggsave(fn_fig3, p.tmp, width = 6, height = 9)
  p[[n.monitor]] <- p.tmp
  gc()
}
p <- gridExtra::marrangeGrob(p[c(1, 3)], nrow=1, ncol=2, top = NULL)
ggsave(gsub("n\\_m.*png", "com.png", fn_fig3), p, width = 12, height = 7)


# ### prediction
# setting(scenario = 1, sim = 1, n.monitor, ntree, 0)
# pred <- fnfun_pred(scenario, n.monitor, sim = 1, fn = fn_output)
# i=4
# plotRF(pred$w132, i = i, pred = "test", truth = s.test[i, ], truth.time = Grid, tau = tau)
# # i=1
# # plotRF(pred$w132, i = i, pred = "train", truth = s.train[i, ], truth.time = Grid, tau = tau)
# plotRF(pred$FuRF1, i = i, pred = "test", truth = s.test[i, ], truth.time = Grid, tau = tau)
# plotRF(pred$FuRF2, i = i, pred = "test", truth = s.test[i, ], truth.time = Grid, tau = tau)
# plotRF(pred$cox, i = i, pred = "test", truth = s.test[i, ], truth.time = Grid, tau = tau)
# 
# method.labs <- list(levs = c("cox", "fu", "fuRF", "fuRF.sm", 
#                              "icrf_local_best", "icrf_global_best", 
#                              "icrf_local_best_imse2", "icrf_global_best_imse2",
#                              paste0("icrf", 1:10)),
#                     labs = c("Cox", "Fu", "FuRF", "FuRF.smooth", 
#                              "icrf_l1", "icrf_g1", "icrf_l2", "icrf_g2", paste0("icrf", 1:10)))
# 
# for (scenario in 1:6) {
#   for (n.monitor in c(1,3)) {
#     perf.list <- list()
#     png.name <- paste0(path_figure, test.name, "_sim_scenario_", scenario, "-n_m-", n.monitor, "-eval.png")
#     cat(png.name, "\n")
#     for (sim in 1:100) {
#       cat(sim, " ")
#       tmp.file <- paste0(output.path, "sim_scenario-", scenario,"-n_m-", n.monitor, "-", test.name, "-eval-", sim,".rds")
#       if (file.exists(tmp.file)) tmp2 <- readRDS(tmp.file) else tmp2 <- tmp
#       if (is.null(tmp2)) tmp2 <- tmp
#       perf.list[[sim]] <- tmp2
#     }
#     # perf.list
#     
#     # into a long form
#     perf.list <- 
#       lapply(1:length(perf.list), function(x) {
#         perf.list[[x]] %>% 
#           add_rownames("method") %>% 
#           mutate(sim = x) %>% 
#           tidyr::gather(key = "metric", value = "value", - sim, -method )
#       })
#     
#     # rbind
#     perf.gathered <- do.call(rbind, perf.list)
#     perf.gathered$methodClass <- NA
#     perf.gathered$methodClass[perf.gathered$method == "cox"] <- "Cox"
#     perf.gathered$methodClass[perf.gathered$method == "fu"] <- "Fu-tree"
#     perf.gathered$methodClass[perf.gathered$method == "fuRF"] <- "Fu-RF"
#     perf.gathered$methodClass[perf.gathered$method == "fuRF.sm"] <- "Fu-RF-sm"
#     perf.gathered$methodClass[grep("icrf[0-9]", perf.gathered$method)] <- "ICRF"
#     perf.gathered$methodClass[grep("icrf\\_", perf.gathered$method)] <- "ICRF-best"
#     
#     perf.gathered %>%
#       # dplyr::filter(! method %in% c("(reserved1)", "(reserved2)")) %>% 
#       # dplyr::filter(metric != "int.sq.err") %>% 
#       dplyr::filter(metric != "best.fold") %>% 
#       dplyr::mutate(method = factor(method, levels = method.labs$levs, labels = method.labs$labs)) %>% 
#       ggplot(aes(method, value, col = methodClass)) +
#       geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.1, width = 0.3, height = 0, size= 0.1, alpha = 0.3) + 
#       facet_grid(~metric) +
#       theme(axis.text.x = element_text(angle = 90))
#     ggsave(png.name, width = 20, height = 15, units = "cm")  
#     gc()
#   }
# }
# 
# # readRDS("../output_tmp/sim_scenario-1-eval-4.rds")
# # readRDS("../output/sim_scenario-1-rep-4.rds")[[1]]
# 
# ################################ prediction ################################
# ################################ prediction ################################
# 
# ### example
# scenario = 1
# sim = 1
# icrf <- readRDS(paste0(output.path, "sim_scenario-", scenario, "-n_m-", n.monitor,  "-icrf-", test.name, "-rep-", sim, ".rds"))
# cox <- readRDS(paste0(output.path, "sim_scenario-", scenario, "-n_m-", n.monitor,  "-cox-rep-", sim, ".rds"))
# ltrc <- readRDS(paste0(output.path, "sim_scenario-", scenario, "-n_m-", n.monitor,  "-ltrc-rep-", sim, ".rds"))
# ltrcForest <- readRDS(paste0(output.path, "sim_scenario-", scenario, "-n_m-", n.monitor,  "-ltrcForest-rep-", sim, ".rds"))
# 
# test; s.test  # get s.test from somewhere in C04-sim-run.R
# icrf.hat <- readRDS(paste0(output.path, "sim_scenario-", scenario, "-n_m-", n.monitor,  "-icrf-", test.name, "-pred-rep-", sim, ".rds"))
# cox.hat <- predict.icensReg(cox, newdata = test, grid = grid)
# ltrc.hat <- predict.ltrc(ltrc, test, grid = grid)
# ltrcForest.hat <- predict.ltrc(ltrcForest, test, grid = grid, smooth = 1, bandwidth = tau * n^{-1/5}/2)
# 
# i = 2 # subject id in test set.
# icrf:::predict.icrf(icrf, Test[i,], time.points = seq(0, tau, by = 0.01)) %>%
#   # lapply(1:10, function(x) data.frame(time = icrf.hat[[x]]$time_interest, fold = x, S = icrf.hat[[x]]$Surv_predict[i, ])) %>% 
#   #   {do.call(rbind, .)} %>% 
#   rbind(data.frame(time = grid, fold = "truth", S = s.test[i, ])) %>% 
#   rbind(data.frame(time = grid, fold = "cox", S = cox.hat[i, ])) %>% 
#   rbind(data.frame(time = grid, fold = "Fu", S = ltrc.hat[i, ])) %>% 
#   rbind(data.frame(time = grid, fold = "FuRF", S = ltrcForest.hat$S[i, ])) %>% 
#   rbind(data.frame(time = grid, fold = "FuRF-Sm", S = ltrcForest.hat$S.smooth[i, ])) %>% 
#   mutate(class = ifelse(!is.na(as.numeric(fold)), "ICRF", fold),
#          fold = factor(fold, levels = c("truth", "cox", "Fu", "FuRF", "FuRF-Sm", 1:10),
#                        labels = c("truth", "cox", "Fu", "FuRF", "FuRF-Sm", paste0("icrf", 1:10)))) -> pred.tmp 
# 
# output.path = "../output_lr2"; scenario = 3; n.monitor = 3; sim = 1; test.name = "logRank2"; test.2sample = 3
# icrf.hat <- readRDS(paste0(output.path, "sim_scenario-", scenario, "-n_m-", n.monitor,  "-icrf-", test.name, "-pred-rep-", sim, ".rds"))
# lapply(1:10, function(x) data.frame(time = icrf.hat[[x]]$time_interest, fold = paste0("icrf", x), S = icrf.hat[[x]]$Surv_predict[i, ])) %>% 
# {do.call(rbind, .)} %>% 
#   rbind(data.frame(time = grid, fold = "truth", S = s.test[i, ])) -> pred.tmp
# 
# pred.tmp %>%
#   # dplyr::filter(!fold %in% c(2:3, 5:6, 8:9)) %>%
#   dplyr::filter(!fold %in% paste0("icrf", c(2:4, 6:9))) %>%
#   ggplot(aes(time, S, color = fold)) + geom_line() + ylab("survival") +
#   # scale_linetype_manual(values=c(truth = 1, cox = 1,  Fu = 1, FuRF = 1, "FuRF-Sm" = 1, 
#   #                                "icrf1" = 1, "icrf2" = 2, "icrf3" = 2, "icrf4" = 2, "icrf5" = 3, 
#   #                                "icrf6" = 3, "icrf7" = 3, "icrf8" = 4, "icrf9" = 4, "icrf10" =4)) +
#   scale_size_manual(values=c(truth = 1, cox = 1,  Fu = 1, FuRF = 1, "FuRF-Sm" = 1, 
#                              "icrf1" = 1, "icrf2" = .95, "icrf3" = .9, "icrf4" = .85, "icrf5" = .8, 
#                              "icrf6" = .75, "icrf7" = .7, "icrf8" = .65, "icrf9" = .6, "icrf10" = .55)) +
#   scale_color_manual(name = "method (fold)",
#                      values=c(truth = "maroon2", cox = "red3",  Fu = "khaki4", FuRF = "springgreen4", "FuRF-Sm" = "springgreen1",  
#                               "icrf1" = "grey70", "icrf2" = "grey60", "icrf3" = "grey55", "icrf4" = "grey45", 
#                               "icrf5" = "grey40", "icrf6" = "grey30", "icrf7" = "grey25", "icrf8" = "grey15", 
#                               "icrf9" = "grey10", "icrf10" = "grey0")) +
#   guides(linetype = FALSE, size = FALSE) +
#   ggtitle(paste0(titles, ", train ", sim, ", ", "predicted values for New subject ", i)) + 
#   geom_segment(aes(x = 0, xend = test[i, "L"], y = 1, yend = 1), col = "gray", linetype = "dotted") +
#   geom_segment(aes(xend = max(tau, test[i, "R"]), x = test[i, "R"], y = 0, yend = 0), col = "gray", linetype = "dotted") +
#   geom_point(aes(x = test[i, "L"], y = 1), col = "gray", fill = "white", shape = 21) +
#   geom_point(aes(x = test[i, "R"], y = 0), col = "gray", fill = "white", shape = 21)
# 
# ggsave(paste0(path_figure, "sim_scenario-", scenario,"-rep-", sim, "-pred-", i, ".png"))  
# gc()
