date = "2021-05-10"
ntree = 100; n.sim=300; scenario = 1; n.monitor = 1

library(ggplot2); library(dplyr); library(purrr)
library(icenReg); library(MASS)
source("scripts/0functions.R")
source("scripts/2analysis-functions.R")
path_output   <- paste0("output/", date, "/")
path_figure    <- paste0("figure/", date, "/")
if (!dir.exists(path_output)) {
  dir.create(path_output)
  warning("No such folder. Created one.")
}
if (!dir.exists(path_figure)) {
  dir.create(path_figure)
  warning("No such folder. Created one.")
}

fn_eval_tmp   <- paste0(path_output, "sizeEval_size_", "sz", "-scenario_", "scn", "-n.m_", "nm",
                        "-nT_", ntree, "-rep_", "sm", ".rds")
fn_pred_tmp   <- paste0(path_output, "sizeSim_size_", "sz", "-scenario_", "scn", "-n.m_", "nm",
                        "-nT_", ntree, "-rep_", "sm", ".rds")

tmp <- fnfun(1, 1, 1, size = 100, fn = fn_eval_tmp)
tmp[,,]<- NA
dm <- dim(tmp)
# tmp <- fnfun(1000,12,12)  #NULL example

# Column vectors
grpVec <- do.call(expand.grid, dimnames(tmp))
names(grpVec) <- c("fold", "measure", "method")
fn_fig3 = paste0(path_figure, "figSize_total-scenario_", scenario, "-n_m_", n.monitor,".png")

## GWRS-Quasi_honesty plot
lvs1 <- c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "w132" , "w135")
lbs1 <- c("Cox", "Cox (*)", "Fu", "Fu (*)", "Yao", "Yao (*)", "ICRF", "ICRF_E")
lvs2 <- c(paste0(c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2"), "-1"), 
          paste0("w132-", c(1:10, "A", "B", "C")), paste0("w135-", c(1:10, "A", "B", "C")))
lbs2 <- c("Cox", "Cox (*)", "Fu", "Fu (*)", "Yao", "Yao (*)", 
          paste0("ICRF-", c(1:10, "A", "B", "C")), paste0("ICRF_E-", c(1:10, "A", "B", "C")))
shp_method = c(0, 15, 1, 16, 2, 17, 5, 18)
names(shp_method) = c(lbs2[1:6], "ICRF-1", "ICRF-10")
lty_method = rep(c("dotted", "solid"), 4)
names(lty_method) = c(lbs2[1:6], "ICRF-1", "ICRF-10")

lvs3 <- c("int.error", "sup.error", "imse.type1 (oob)", "imse.type2 (oob)", "imse.type1", "imse.type2")
lbs3 <- c("integrated error", "supremum error", "IMSE1 (out of bag)", "IMSE2 (out of bag)", "IMSE1", "IMSE2")
lvs4 <- 1:6
lbs4 <- paste0("scenario ", lvs4)
lvs5 <- lbs5 <- c(100, 200, 400, 800, 1600)

grandResult <- grandResult.mean <- grandResult.min <- data.frame()
# for (n.monitor in c(1, 3)) {
  for (size in lvs5) {
    cat("n.monitor: ", n.monitor, ", n = ", size,", scenario = ", scenario, "\n")
    result <- 
      lapply(1:n.sim, function(s) {
        tmp1 <- fnfun(scenario, n.monitor, sim = s, fn = fn_eval_tmp, size = size)
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
      dplyr::filter(grepl("(^w|cox|Fu)", method)) %>%
      mutate(methods = paste0(method, "-", fold))
    
    if (dim(result)[1] == 0) next
    
    result$method = factor(result$method, levels = lvs1, labels = lbs1)
    result$methods = factor(result$methods, levels = lvs2, labels = lbs2)
    result$measure = factor(result$measure, levels = lvs3, labels = lbs3)
    result.mean <- aggregate(value ~ methods + measure + method + fold, data = result, mean, na.rm = T)
    result.1Q <- aggregate(value ~ methods + measure + method + fold, data = result, quantile, 0.25, na.rm = T)
    result.3Q <- aggregate(value ~ methods + measure + method + fold, data = result, quantile, 0.75, na.rm = T)
    result.min <- aggregate(value ~ measure, data = result.mean, min, na.rm = TRUE)
    names(result.mean)[names(result.mean) == "value"] = "mean"
    names(result.1Q)[names(result.1Q) == "value"] = "Q1"
    names(result.3Q)[names(result.3Q) == "value"] = "Q3"
    names(result.min)[names(result.min) == "value"] = "min"
    result.mean <- 
      left_join(result.mean, result.1Q, by = c("measure", "methods", "method", "fold")) %>% 
      left_join(result.3Q, by = c("measure", "methods", "method", "fold")) %>% 
      left_join(result.min, by = "measure")
      
    result$size = size
    result.mean$size = size
    result.min$size = size
    
    grandResult <- rbind(grandResult, result)
    grandResult.mean <- rbind(grandResult.mean, result.mean)
    grandResult.min <- rbind(grandResult.min, result.min)
  }
#}
grandResult$size = factor(grandResult$size, levels = lvs5, labels = lbs5)
grandResult.mean$size = factor(grandResult.mean$size, levels = lvs5, labels = lbs5)
grandResult.min$size = factor(grandResult.min$size, levels = lvs5, labels = lbs5)
pd <- position_dodge(0.2)
grandResult.mean %>% 
  dplyr::filter(fold %in% c(1, 10), method != "ICRF_E") %>%
  dplyr::filter(grepl("error", measure)) %>% 
  ggplot(aes(x = size, y = mean, col = methods, group = methods)) +
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.1, position = pd) + 
  geom_line(position = pd, alpha = 0.3) +  
  geom_point(aes(shape = methods), size = 2, position = pd) +
  scale_shape_manual(values = shp_method) +
  scale_linetype_manual(values = lty_method) +
  facet_grid(.~ measure , scales = "free") + guides(label = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("sample size") + ylab("prediction errors (1Q, mean, 3Q)")

ggsave(fn_fig3, width = 8, height = 5)
gc()

### 1.2. computation time comparison per sample sizes

  fn_fig_time = paste0(path_figure, "figSize_time_n_m_", n.monitor,".png")
  grandResultTime <- data.frame()
  for (size in lvs5) {
    print(size)  
    resultTime <- 
      sapply(1:n.sim, function(s) {
        cat(s, " ")
        tmp1 <- fnfun(scenario, n.monitor, sim = s, fn = fn_eval_tmp, size = size)
        if (is.null(tmp1)) {
          warning(paste0(s, " is not available."))
          tmp1 <- NULL
        }
        tmp1 <- attr(tmp1, "runtime")["w132"]
        if (is.null(tmp1)) tmp1 <- NA
        tmp1
      }) %>% data.frame(time = ., scenario = scenario, sim = 1:n.sim, size = size)
    grandResultTime <- rbind(grandResultTime, resultTime)
  }
  
  grandResultTime.gg <-
    grandResultTime %>% 
    #mutate(size = size) %>% 
    group_by(scenario, size) %>% 
    summarize(time  = mean(time, na.rm = TRUE)) %>% 
    ungroup 
  
  grandResultTime.gg %>% 
    ggplot(aes(size, time)) + 
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = FALSE) +
    # geom_abline(slope = 7.79/800, intercept = 0, col = "red") +
    # geom_point(size = 3, col = "magenta") +
    geom_bar(stat = "identity", mapping = aes(fill = (size))) +
    guides(fill = FALSE) +
    xlab("sample size") +
    scale_x_continuous(breaks = c(100, 200, 400, 800, 1600)) +
    # scale_y_log10() +
    theme_bw() + ylab("time (in minutes)")
  ggsave(fn_fig_time,  width = 8, height = 6)

