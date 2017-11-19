#' Single Variable Rejection Sampler
#'
#' This function implements single variable rejection sampling for rvs with bounded support and which have bounded pdf.
#'
#'
#' @param fun The pdf that we are sampling from, input as a string. The variable sampled from must be input using the variable "x".
#'
#' @param a The lower bound of the pdf.
#'
#' @param b The upper bound of the pdf.
#'
#' @param N The number of samples output by the sampling function.
#'
#' @return A vector containing the samples from the pdf.
#' @export
#'
#' @examples
#'
#' hist(samp1D("2*x", 0, 1, 10000))
#'

samp1D <- function(fun,a,b,N) {
  g <- parse(text = fun)
  f <- function(x) {
    ifelse(a <= x & x <= b, eval(g[[1]]), 0)
  }
  maxf <- optimize(f,c(a,b), tol = 0.0001, maximum = TRUE)
  samples <- rep(0,N)
  i <- 0
  while ( i < N) {
    potsamp <- runif(1,a,b)
    testsamp <- runif(1, 0, maxf$objective + .1)
    if ( testsamp < f(potsamp) ) {
      samples[i+1] = potsamp
      i = i + 1
    }
  }
  samples
}
