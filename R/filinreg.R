#' Fiducial sampler for linear regression model
#' @description ddd
#'
#' @param y xx
#' @param X xx
#' @param L xx
#' @param lucky xx
#'
#' @examples xxx
#'
#' @importFrom arrangements icombinations
#' @importFrom EigenR Eigen_rank Eigen_inverse
#' @importFrom utils head
#' @export
filinreg <- function(
  y, X = as.matrix(rep(1,length(y))), df = Inf, L = 10L, lucky = FALSE
){
  qdistr <- function(x, ...) qt(x, df=df, ...)
  ddistr <- function(x, ...) dt(x, df=df, ...)
  X <- unname(X)
  n <- nrow(X)
  p <- ncol(X)
  q <- p + 1L
  # centers of hypercubes (volume 1/L^p)
  centers <- as.matrix(
    do.call(
      expand.grid, rep(list(seq(0, 1, length.out = L+1L)[-1L] - 1/(2*L)), q)
    )
  )
  # remove centers having equal coordinates (H'H is not invertible)
  centers <-
    centers[apply(centers, 1L, function(row) length(unique(row)) > 1L),]
  # outputs
  M <- (L^q - L) / 2L # number of centers yielding sigma>0
  J <-  rep(NA_real_, M)
  Theta <- matrix(NA_real_, nrow = M, ncol = q)
  # algorithm
  Iiterator <- icombinations(n, q)
  I <- Iiterator$getnext()
  XI <- X[I, ]
  while(Eigen_rank(XI) < p){ #TODO check X full rank at the beginning
    I <- Iiterator$getnext()
    XI <- X[I, ]
  }
  XmI <- X[-I,]
  yI <- y[I]
  ymI <- y[-I]
  counter <- 0L
  if(lucky){
    for(m in 1L:nrow(centers)){
      H <- cbind(XI, qdistr(centers[m, ]))
      theta <- Eigen_inverse(crossprod(H)) %*% t(H) %*% yI
      if(theta[q] > 0){ # sigma>0
        counter <- counter + 1L
        J[counter] <-
          sum(ddistr((ymI - XmI %*% head(theta, -1L))/theta[q], log = TRUE)) -
          (n-q) * log(theta[q])
        Theta[counter,] <- theta
      }
    }
  }else{
    for(m in 1L:nrow(centers)){
      H <- cbind(XI, qdistr(centers[m, ]))
      if(Eigen_rank(H) < q){
        Theta <- head(Theta, -1L)
        J <- head(J, -1L)
        next
      }
      theta <- Eigen_inverse(crossprod(H)) %*% t(H) %*% yI
      if(theta[q] > 0){ # sigma>0
        counter <- counter + 1L
        J[counter] <-
          sum(ddistr((ymI - XmI %*% head(theta, -1L))/theta[q], log = TRUE)) -
          (n-q) * log(theta[q])
        Theta[counter,] <- theta
      }
    }
  }
  J <- exp(J)
  list(Beta = Theta[, -q], sigma = Theta[, q], W = J/sum(J))
}
