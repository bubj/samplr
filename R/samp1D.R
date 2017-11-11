#' Single Variable Rejection Sampler
#'
#' This function immplents single variable rejection sampling for rvs with bounded support and which have bounded pdf.
#'
#' This paragrpah will show up somewhere else for additional information.
#'
#' @param fun the pdf that we are sampling from, input as an expression
#'
#' @param a the lower bound on the function f
#'
#' @param b the upper bound on the function f
#'
#' @param lb the lower bound of support f
#'
#' @param ub upper bound of support of f
#'
#' @param N the number of attempted samples.
#'
#' @return A vector containg the samples form pdf
#' @export
#'
#' @examples
#'
#' hist(samp1D(expression(2*x), 0, 1, 0, 1, 10000))
#'

samp1D <- function(fun,a,b,lb,ub,N) {
  f <- function(x) {
    ifelse(a < x & x < b, eval(fun[[1]]), 0)
  }
  maxf <- max(f(seq(a,b,.01))+1)
  ones <- runif(N, lb, ub)
  unis <- runif(N, 0, maxf)
  ones[unis < f(ones)]
}
