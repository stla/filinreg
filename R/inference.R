inference <- function(fidsamples, param, alpha = 0.05){
  fipred <- inherits(fidsamples, "filinreg.pred")
  out <- numeric(4L)
  names(out) <- c("mean", "median", "lwr", "upr")
  sample <-
    if(fipred) fidsamples[["FPD"]][[param]] else fidsamples[["Theta"]][[param]]
  weights <- fidsamples[["weight"]]
  out[1L] <- sum(sample * weights) # mean
  h <- cbind(sample, weights)
  hsort <- h[order(h[,1L]), ]
  hsum <- cumsum(hsort[, 2L])
  ci_u <- min(which(hsum >= 1-alpha/2))
  ci_l <- min(which(hsum >= alpha/2))
  ci_m <- min(which(hsum >= 0.5))
  out[3L] <- hsort[ci_l, 1L] # lower bound
  out[4L] <- hsort[ci_u, 1L] # upper bound
  out[2L] <- hsort[ci_m, 1L] # estimate (median)
  out
}
