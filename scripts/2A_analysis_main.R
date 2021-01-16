date = "2021-01-13"
ntree = 300; n.sim=300

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

fn_eval_tmp   <- paste0(path_output, "eval_scenario_", "scn", "-n.m_", "nm",
                        "-nT_", ntree, "-rep_", "sm", ".rds")
fn_pred_tmp   <- paste0(path_output, "sim_scenario_", "scn", "-n.m_", "nm",
                        "-nT_", ntree, "-rep_", "sm", ".rds")

  tmp <- fnfun(1, 1, 1, fn = fn_eval_tmp)
  tmp[,,]<- NA
  dm <- dim(tmp)
  dmn <- dimnames(tmp)
  grpVec <- do.call(expand.grid, dimnames(tmp))
  names(grpVec) <- c("fold", "measure", "method")
  lvs2 <- c(paste0(c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2"), "-1"), 
            paste0("ICRF", "-", c(1:10, "A", "B", "C")))
  lbs2 <- c("Cox", "Cox (smooth)", "STIC", "STIC (smooth)", "SFIC", "SFIC (smooth)", 
            paste0("ICRF-", c(1:10, "A", "B", "C")))
  pd <- position_dodge(0.2)
  
for (n.monitor in c(1, 3)) {
  fn_fig = paste0(path_figure, "fig_split_n_m_", n.monitor,".png")
  for (scenario in 1:6) {
    print(scenario)  
    result <- 
      lapply(1:n.sim, function(s) {
        cat(s, " ")
        tmp1 <- fnfun(scenario, n.monitor, sim = s, fn = fn_eval_tmp)
        if (is.null(tmp1)) stop("null")
        if (!identical(dim(tmp1), dm)) {
          tmp1 <- tmp1[,,dimnames(tmp)[[3]]]     # in case we only ran some of the methods.
          if (!identical(dim(tmp1), dm)) stop("wrong dimension")
        }
        tmp1 <- as.vector(tmp1)
        tmp1
      }) %>% do.call(c, .)
    
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
      na.omit %>% 
      dplyr::filter(grepl("int\\.error", measure)) %>% 
      rbind(result.A) %>% 
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
    dplyr::filter(fold %in% c(1, 2, 3, 5, 10, "A")) %>%
    dplyr::filter(grepl("(p|w|l)", method)) %>%
    dplyr::filter(grepl("int\\.error", measure)) %>%
    mutate(honesty = ifelse(grepl("H", method), "quasi-honest", 
                            ifelse(grepl("E", method), "exploitative",
                                   "best IMSE1"))) %>% 
    mutate(split = ifelse(grepl("^w", method), "GWRS", 
                          ifelse(grepl("^l", method), "GLR",
                                 ifelse(grepl("^sw", method), "SWRS", "SLR")))) %>% 
    mutate(methods = paste0("ICRF-", fold)) %>% 
    mutate(methods = factor(methods, levels = lvs2, labels = lbs2 )) %>% 
    group_by(fold, measure, split, methods, honesty, scenario) %>% 
    summarize(mean  = mean(value, na.rm = TRUE), 
              Q1 = quantile(value, 0.25, na.rm = TRUE), 
              Q3 = quantile(value, 0.75, na.rm = TRUE)) %>% 
    ungroup 
    
  grandResult.gg %>%   
    ggplot(aes(x = methods, y = mean, col = honesty, group = honesty)) +
    geom_point(aes(shape=honesty)) +
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.1, alpha = 0.3, position = pd) + 
    geom_hline (data = grandResult.mean, aes(yintercept = min), color = "black") +
    facet_grid(scenario ~ split, scale = "free") + guides(label = FALSE) +
    geom_line(data = grandResult.gg %>% filter(fold != "A"), alpha = 0.3) +
    theme_bw() +  xlab("iterations") +
    ylab("integrated error (1Q, mean, 3Q)") + 
    theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") # +
  ggsave(fn_fig,  width = 6, height = 5)
  gc()
  
}
  
  
## GWRS-Quasi_honesty plot
type = "wH"
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
    dplyr::filter(fold %in% c(1:3, 5, 10, "A")) %>% 
    mutate(value = ifelse(grepl("IMS", measure), value * imse.multiplier, value)) %>% 
    ggplot(aes(x = methods, y = value, col = method)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(size = 0.1, width = 0.3, height = 0, size= 0.1, alpha = 0.3) +
    stat_summary(fun.y=mean, geom="point", shape = 22, size = 0.5, col = 'deeppink4', fill = "white") +
    geom_hline (data = grandResult.mean %>% 
                  mutate(value = ifelse(grepl("IMSE", measure), value * imse.multiplier, value)) %>% 
                  mutate(min = ifelse(grepl("IMSE", measure), min * imse.multiplier, min)), 
                aes(yintercept = min), color = "black") +
    scale_y_continuous("integrated and supremum error", 
                       sec.axis = sec_axis(~ . / imse.multiplier, name = "IMSE")) +
    facet_grid(scenario ~ measure, scales = "free") + guides(label = FALSE, col = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) -> p.tmp
  ggsave(fn_fig3, p.tmp, width = 6, height = 9)
  p[[n.monitor]] <- p.tmp
  gc()
}
p <- gridExtra::marrangeGrob(p[c(1, 3)], nrow=1, ncol=2, top = NULL)
ggsave(gsub("n\\_m.*png", "com.png", fn_fig3), p, width = 12, height = 7)
