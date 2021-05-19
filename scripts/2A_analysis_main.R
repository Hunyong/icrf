date = "2021-05-10"
ntree = 300; n.sim=300

### 0.0 libraries and paths
  library(ggplot2); library(dplyr); library(purrr)
  library(icenReg); library(MASS)
  library(ggpattern); #remotes::install_github("coolbutuseless/ggpattern")
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

### 0.1 file names, labels
  fn_eval_tmp   <- paste0(path_output, "eval_scenario_", "scn", "-n.m_", "nm",
                          "-nT_", ntree, "-rep_", "sm", ".rds")
  fn_pred_tmp   <- paste0(path_output, "sim_scenario_", "scn", "-n.m_", "nm",
                          "-nT_", ntree, "-rep_", "sm", ".rds")
  lvs1 <- c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2", "wH")     # wH is for ICRF with G[W]RS + quasi[H]onesty
  lbs1 <- c("Cox", "Cox (*)", "Fu", "Fu (*)", "Yao", "Yao (*)", "ICRF")
  lvs2 <- c(paste0(c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2"), "-1"), 
            paste0("ICRF", "-", c(1:10, "A", "B", "C")))
  lbs2 <- c("Cox", "Cox (*)", "Fu", "Fu (*)", "Yao", "Yao (*)", 
            paste0("ICRF-", c(1:10, "A", "B", "C")))
  lvs2B <- c(paste0(c("cox", "cox.sm", "FuTR1", "FuTR2", "FuRF1", "FuRF2"), "-1"), 
            paste0("wH", "-", c(1:10, "A", "B", "C")))
  lbs2B <- c("Cox", "Cox (*)", "Fu", "Fu (*)", "Yao", "Yao (*)", 
            paste0("ICRF-", c(1:10, "A", "B", "C")))
  lvs3 <- c("int.error", "sup.error", "imse.type1 (oob)", "imse.type2 (oob)", "imse.type1", "imse.type2")
  lbs3 <- c("integrated error", "supremum error", "IMSE1 (out of bag)", "IMSE2 (out of bag)", "IMSE1", "IMSE2")
  lvs4 <- 1:6
  lbs4 <- paste0("scenario ", lvs4)
  
  # computation
  lvs5 = c("wH", "lH", "wE", "lE")
  lvs5 = c(lvs5, paste0("s", lvs5))
  lbs5 = c("GWRS", "GLR", "GWRS", "GLR", "SWRS", "SLR", "SWRS", "SLR")
  lbs5b = rep(c(rep("quasi-honest", 2), rep("exploitative", 2)), 2)


### 0.2 skeleton and constants
  tmp <- fnfun(1, 1, 1, fn = fn_eval_tmp)
  tmp[,,]<- NA
  dm <- dim(tmp)
  dmn <- dimnames(tmp)
  grpVec <- do.call(expand.grid, dimnames(tmp))
  names(grpVec) <- c("fold", "measure", "method")
  pd <- position_dodge(0.45)
  
  
  
  
### 1 Results
### 1.1. Splitting rules comparison
  for (n.monitor in c(1, 3)) {
    fn_fig = paste0(path_figure, "fig_split_n_m_", n.monitor,".png")
    grandResult <- grandResult.mean <- grandResult.min  <- data.frame()
    for (scenario in 1:6) {
      print(scenario)  
      result <- 
        lapply(1:n.sim, function(s) {
          cat(s, " ")
          tmp1 <- fnfun(scenario, n.monitor, sim = s, fn = fn_eval_tmp)
          if (is.null(tmp1)) {
            warning(paste0(s, " is not available."))
            tmp1 <- tmp
            # stop("null")
          }
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
        dplyr::filter(!is.na(`imse.type1 (oob)`)) %>%
        group_by(sim, method) %>%
        dplyr::filter(`imse.type1 (oob)`== min(`imse.type1 (oob)`)) %>% # find the best imse1.
        mutate(fold = "A", value = int.error, measure = "int.error") %>%
        ungroup %>%
        dplyr::select(value, fold, measure, method, sim)
      
      
      result <- 
        data.frame(value = result, grpVec, sim = rep(1:n.sim, each = dim(grpVec)[1])) %>% 
        na.omit %>% 
        dplyr::filter(grepl("int\\.error", measure)) %>% 
        rbind(result.A) %>% 
        mutate(methods = paste0(method, "-", fold))
      
      if (dim(result)[1] == 0) next
      result.mean <- aggregate(value ~ methods + measure + method, data = result, mean, na.rm = T)
      result.min <- aggregate(value ~ measure, data = result.mean, min, na.rm = TRUE)
      names(result.min)[2] = "min"
      result.mean <- left_join(result.mean, result.min, by = "measure")
      
      result$scenario = scenario
      result.mean$scenario = scenario
      result.min$scenario = scenario
      
      grandResult <- rbind(grandResult, result)
      grandResult.mean <- rbind(grandResult.mean, result.mean)
      grandResult.min <- rbind(grandResult.min, result.min)
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
      mutate(split = factor(split, levels = c("GWRS", "GLR", "SWRS", "SLR"))) %>% 
      mutate(methods = paste0("ICRF-", fold)) %>% 
      mutate(methods = factor(methods, levels = lvs2, labels = lbs2 )) %>% 
      group_by(fold, measure, split, methods, honesty, scenario) %>% 
      summarize(mean  = mean(value, na.rm = TRUE), 
                Q1 = quantile(value, 0.25, na.rm = TRUE), 
                Q3 = quantile(value, 0.75, na.rm = TRUE)) %>% 
      ungroup
      
    # grandResult.gg %>% 
    #   # left_join(grandResultTime.gg) %>% 
    #   ggplot(aes(x = methods, y = mean, col = honesty, group = honesty)) +
    #   geom_point(aes(shape=honesty)) +
    #   geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.1, alpha = 0.3, position = pd) + 
    #   geom_hline (data = grandResult.mean, aes(yintercept = min), color = "black") +
    #   facet_grid(scenario ~ split, scale = "free") + guides(label = FALSE) +
    #   geom_line(data = grandResult.gg %>% filter(fold != "A"), alpha = 0.3) +
    #   theme_bw() +  xlab("iterations") +
    #   ylab("integrated error (1Q, mean, 3Q)") + 
    #   # geom_text(aes(label = ifelse(is.na(time), NA, paste0("time = ", round(time + 0.15,1))))) +
    #   theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    
    grandResult.gg %>% 
      ggplot(aes(x = methods, y = mean, col = split, group = split, linetype = split)) +
      geom_point(aes(shape=split), position = pd, size = 2) +
      geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.1, alpha = 0.8, position = pd) + 
      geom_hline (data = grandResult.mean, aes(yintercept = min), color = "black") +
      geom_line(data = grandResult.gg %>% filter(fold != "A")) +
      facet_grid(scenario ~ honesty, scale = "free") + guides(label = FALSE, linetype = "legend") +
      theme_bw() +  xlab("iterations") +
      ylab("integrated error (1Q, mean, 3Q)") + 
      # geom_text(aes(label = ifelse(is.na(time), NA, paste0("time = ", round(time + 0.15,1))))) +
      theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
    ggsave(fn_fig,  width = 8, height = 6)
    gc()
  }

  
  
### 1.2. Splitting rules computation time comparison
  for (n.monitor in c(1, 3)) {
    fn_fig_time = paste0(path_figure, "fig_split_time_n_m_", n.monitor,".png")
    grandResultTime <- data.frame()
    for (scenario in 1:6) {
      print(scenario)  
      resultTime <- 
        sapply(1:n.sim, function(s) {
          cat(s, " ")
          tmp1 <- fnfun(scenario, n.monitor, sim = s, fn = fn_eval_tmp)
          if (is.null(tmp1)) stop("null")
          tmp1 <- attr(tmp1, "runtime")[lvs5]
          tmp1
        }) %>% t %>% data.frame(scenario = scenario, sim = 1:n.sim)
      grandResultTime <- rbind(grandResultTime, resultTime)
    }
    
    grandResultTime <-
      grandResultTime %>% 
      mutate(scenario = factor(scenario, levels = lvs4, labels = lbs4)) %>% 
      tidyr::gather(key = "method", value = "time", -scenario, -sim) %>% 
      mutate(honesty = factor(method, levels = lvs5, lbs5b),
             split = factor(method, levels = lvs5, lbs5))
    
    grandResultTime.gg <-
      grandResultTime %>% 
      group_by(split, honesty, scenario) %>% 
      summarize(methods = factor("ICRF-10", levels = lvs2, labels = lbs2),
                time  = mean(time, na.rm = TRUE)) %>% 
      ungroup %>% 
      arrange(scenario)
    
    grandResultTime.gg %>% 
      ggplot(aes(honesty, time, fill = split)) +
      # geom_bar(position = "dodge", stat='identity') +
      geom_bar_pattern(position = "dodge", stat='identity',
                       pattern = c("none", "stripe", "circle", "crosshatch")[
                         as.numeric(grandResultTime.gg$split)],
                       pattern_density = 0.1,
                       pattern_spacing = 0.05,
                       pattern_angle = 45,
                       pattern_colour = c("none", "white", "black", "white")[
                         as.numeric(grandResultTime.gg$split)]) +
      facet_grid(scenario ~ .) + 
      theme_bw() + ylab("time (in minutes)") +
      guides(fill = guide_legend(override.aes = 
                                   list(pattern = c("none", "stripe", "circle", "crosshatch"),
                                        pattern_colour = c("none", "white", "black", "white"),
                                        pattern_spacing = .02,
                                        pattern_angle = 45)
      )) +
      theme(axis.title.x = element_blank())
    
    ggsave(fn_fig_time,  width = 8, height = 6)
  }


    
### 1.3 the main plot (based on GWRS-Quasi_honesty)
  p <- list()
  for (n.monitor in c(1, 3)) {
    fn_fig3 = paste0(path_figure, "figWRS_total_wH_n_m_", n.monitor,".png")
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
        dplyr::filter(grepl("(^wH|cox|Fu)", method)) %>%
        mutate(methods = paste0(method, "-", fold))
      
      
      result.icrf.conv <-
        result %>%
        dplyr::filter(method == "wH", measure == "imse.type1 (oob)") %>%
        group_by(method, sim) %>%
        summarize(convA = convMonitor(value, type = "global"),
                  convB = convMonitor(value, type = "local"),
                  convC = convMonitor(value, type = "first")) %>% 
        left_join(result %>% dplyr::filter(method == "wH"), 
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
      result$methods = factor(result$methods, levels = lvs2B, labels = lbs2B)
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
      geom_jitter(size = 0.1, width = 0.3, height = 0, size= 0.1, alpha = 0.2) +
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
    # ggsave(fn_fig3, p.tmp, width = 6, height = 6)
    p[[n.monitor]] <- p.tmp
    gc()
  }
  p <- gridExtra::marrangeGrob(p[c(1, 3)], nrow=1, ncol=2, top = NULL)
  ggsave(gsub("n\\_m.*png", "com.png", fn_fig3), p, width = 12, height = 5.5)
