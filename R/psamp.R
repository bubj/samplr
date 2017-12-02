#' Single and Two Variable Probability Function
#'
#' This function implements one and two variable rejection sampling to find probabilities.
#'
#' This function is in an implimentation of a distribution function for single and two variable pdfs using the \code{samp1D} and \code{samp2D} rejection samplers.
#'
#' For one dimensional variables, this function uses \code{samp1D} to find \code{P(X < q)}.
#'
#' For two dimensional random variables, this function uses \code{samp2D} to find \code{P(g(X,Y) < q)}, where \code{X} and \code{Y} are found using \code{samp2D(f,10000)}.
#'
#' @param f The pdf you wish to use to find the probability of. For 2D probability density functions, the argument must be a vector of the two parameters of the pdf.
#'
#' @param q The value you want the samples to be less than. For 2D, the value you want your function g to be less than.
#'
#' @param g The function of random variables x and y for 2D that you wish to find the probability associated with.
#'
#' @return A number between 0 and 1 associated with the distribution function.
#' @export
#'
#' @examples
#'
#' One dimensional
#'
#' f <- function(x) {
#'   ifelse(0 < x & x < 1, 2*x, 0)
#' }
#' psamp(f,.5)
#'
#' f <- function(x) {
#'   ifelse(0 < x & x < 2, 1/2*x, 0)
#' }
#' psamp(f,1)
#'
#' f <- function(x) {
#'   ifelse(0 < x & x < 6.2832, 1/2/pi*(sin(x) + 1), 0)
#' }
#' psamp(f,pi)
#'
#' Two dimensional
#'
#' f <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   ifelse(0 <= x & x <= 1 & 0 <= y & y <= 1, x + y, 0)
#' }
#' g <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   2*x + y
#' }
#' psamp(f,1,f)
#' psamp(f,2,g)
#'

psamp <- function(f, q, g = NULL) {
  if (is.null(g)) {
    samples <- samp1D(f,10000)
    mean(samples < q)
  } else {
    samples <- samp2D(f,10000)
    values <- rep(0,10000)
    for ( i in 1:10000) {
      values[i] <- g(c(samples[i,1],samples[i,2]))
    }
    mean(values < q) # change this
  }
}
