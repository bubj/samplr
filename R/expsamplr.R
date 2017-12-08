#' Single and Two Variable Expectation
#'
#' This function implements one and two variable rejection sampling to find expected values for one and two dimensional pdfs.
#'
#' For one dimensional variables, this function uses \code{samplr} to find \code{E[X]}
#'
#' For two dimensional random variables, this function uses \code{samplr} to find \code{E[g(X,Y)]}, where \code{X} and \code{Y} are found using \code{samplr(f, 10000, twod = TRUE)}.
#'
#' @param f The pdf you wish to use to find the expectation of. For 2D probability density functions, the argument must be a vector of the two parameters of the pdf.
#'
#' @param g The function of random variables x and y for 2D that you wish to find the expected value of.
#'
#' @return For one dimensional pdfs the ouput is \code{E[X]} and for two dimensional pdfs the ouput is \code{E[g(X,Y)]}.
#' @export
#'
#' @examples
#'
#' One dimensional
#'
#' f <- function(x) {
#'   ifelse(0 < x & x < 1, 2*x, 0)
#' }
#' expsamplr(f)
#'
#' f <- function(x) {
#'   ifelse(0 < x & x < 2, 1/2*x, 0)
#' }
#' expsamplr(f)
#'
#' f <- function(x) {
#'   ifelse(0 < x & x < 6.2832, 1/2/pi*(sin(x) + 1), 0)
#' }
#' expsamplr(f)
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
#'   x
#' }
#' expsamplr(f,f)
#' expsamplr(f,g)
#'


expsamplr <- function(f, g = NULL) {
  if (is.null(g)) {
    samples <- samplr(10000,f)
    mean(samples)
  }
  else {
    samples <- samplr(10000,f, twod = TRUE)
    values <- rep(0,10000)
    for (i in 1:10000) {
      values[i] <- g(c(samples[i,1],samples[i,2]))
    }
    mean(values)
  }
}
