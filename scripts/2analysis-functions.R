fnfun <- function(scenario, n.monitor, sim, fn, size = 200) {
  fn <- gsub("scn", scenario, fn)
  fn <- gsub("nm", n.monitor, fn)
  fn <- gsub("sm", sim, fn)
  fn <- gsub("sz", size, fn)
  if (!file.exists(fn)) {
    warning(paste0("No such file: ", fn))
    return(NULL)
  }
  readRDS(fn)
}
fnfun_pred <- function(scenario, n.monitor, sim, fn, size = 200) {
  fn <- gsub("scn", scenario, fn)
  fn <- gsub("nm", n.monitor, fn)
  fn <- gsub("sm", sim, fn)
  fn <- gsub("sz", size, fn)
  if (!file.exists(fn)) {
    warning(paste0("No such file: ", fn))
    return(NULL)
  }
  readRDS(fn)
}
meanFunction <- function(x)
  data.frame(y=round(mean(x), 3), label=round(mean(x, na.rm=T), 3))
meanPlus1 <- function(x) mean(x, na.rm = TRUE) + 1

convMonitor <- function(vec, type = c("global", "local", "first")) {
  if (length(vec) < 2) return(length(vec))
  type = match.arg(type)
  if (type == "global") {
    return(which.min(vec))
  } else if (type == "local") { # argmin of the first decreasing sequence.
    inc <- c(vec[-1], Inf) - vec
    inc <- ifelse(inc <= 0, 1, 0)
    return(max(1, sum(cumprod(inc))))
  } else { #first: if there is hike at second iteration, stop at first. O/W go through all the iteratons.
    return(ifelse(vec[2] >= vec[1], 1, length(vec)))
  }
}
