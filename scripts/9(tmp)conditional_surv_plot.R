install.packages("icrf")
library(icrf)
?icrf
data(rat2)
data(rat2)

set.seed(1)
rats.icrf <-
  icrf(~ dose.lvl + weight + male + cage.no, data = rat2,
       data.type = "currentstatus", currentstatus.label = c("survtime", "tumor"),
       returnBest = TRUE, ntree=10, nfold=3)
survplot(rats.icrf, c(1,3,5))


tmp.surv <- as.numeric(rats.icrf$predicted[5,])
tmp.surv <- pmin(1, pmax(0, tmp.surv))
tmp.time <- rats.icrf$time.points
tmp.surv <- tmp.surv[is.finite(tmp.time)]
tmp.surv <- (tmp.surv - 0.5)*2
tmp.time <- tmp.time[is.finite(tmp.time)]
tmp.cond <- tmp.surv[tmp.time %in% c(80, 100)]
tmp.surv.cond <- pmin(1, pmax(0, (tmp.surv - tmp.cond[2])/(tmp.cond[1] - tmp.cond[2])))
ggplot2::ggplot(data.frame(time = rep(tmp.time, 2), 
                           s = c(tmp.surv, tmp.surv.cond),
                           type = factor(rep(c("full", "conditional"), each = length(tmp.time)), 
                                         levels = c("full", "conditional")))) + 
  ggplot2::geom_line(ggplot2::aes(time, s, col = type, linetype = type)) + 
  ggplot2::ylab("survival probability") + 
  ggplot2::theme_bw() + ggplot2::ylim(0:1) + ggplot2::ggtitle("Full and interval-conditional survival probabilities")


