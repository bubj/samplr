#' Single Variable Rejection Sampler
#'
#' This function immplents single variable rejection sampling for rvs with bounded support and which have bounded pdf.
#'
#' This paragrpah will show up somewhere else for additional information.
#'
#' @param f the pdf that we are sampling from
#'
#' @param N the number of attempted samples.
#'
#' @param lb the lower bound of support f
#'
#' @param ub upper bound of support of f
#'
#' @param maxf bound on f
#'
#' @return A vector containg the samples form pdf
#' @export
#'
#' @examples
#' betaPDF <- function(x) {
#'     ifelse(0 < x & x < 1, 2*x, 0)
#' }
#' hist(samp1D(betaPDF, 10000, 0, 1, 2))
#'

samp1D <- function(f,N,lb,ub,maxf) {
  ones <- runif(N, lb, ub)
  unis <- runif(N, 0, maxf)
  ones[unis < f(ones)]
}
