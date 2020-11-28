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

#' Title
#'
#' @param fidsamples xx
#' @param conf confidence level
#'
#' @return xx
#' @export
#'
#' @examples xx
fiSummary <- function(fidsamples, conf = 0.95){
  sims <- if(inherits(gfi, "filinreg.pred")){
    fidsamples[["FPD"]]
  }else fidsamples[["Theta"]]
  seq_ <- 1L:ncol(sims)
  names(seq_) <- names(sims)
  out <-
    t(vapply(seq_, function(x) inference(fidsamples, x, 1-conf), numeric(4L)))
  attr(out, "confidence level") <- conf
  out
}

#' Title
#'
#' @param parameter xx
#' @param fidsamples xx
#' @param conf xx
#'
#' @return xx
#'
#' @importFrom lazyeval f_eval_rhs
#' @importFrom spatstat ewcdf quantile.ewcdf
#' @export
#'
#' @examples xx
fiConfInt <- function(parameter, fidsamples, conf = 0.95){
  dataName <- ifelse(inherits(fidsamples, "filinreg.pred"), "FPD", "Theta")
  data <- fidsamples[[dataName]]
  fsims <- f_eval_rhs(parameter, data = data)
  fcdf <- ewcdf(fsims, weights = fidsamples[["weight"]])
  alpha <- 1 - conf
  quantile.ewcdf(fcdf, c(alpha/2, 1-alpha/2))
}
