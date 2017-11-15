#' Single Variable Rejection Sampler
#'
#' This function immplents single variable rejection sampling for rvs with bounded support and which have bounded pdf.
#'
#' This paragrpah will show up somewhere else for additional information.
#'
#' @param fun the pdf that we are sampling from, input as a string
#'
#' @param a the lower bound on the function f
#'
#' @param b the upper bound on the function f
#'
#' @param N the number of attempted samples.
#'
#' @return A vector containg the samples form pdf
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
  ones <- runif(N, a, b)
  unis <- runif(N, 0, maxf$objective + .1)
  ones[unis < f(ones)]
}
