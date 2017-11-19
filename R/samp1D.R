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
  g <- parse(text = fun) # changes the input fun string into an expression
  f <- function(x) {
    ifelse(a <= x & x <= b, eval(g[[1]]), 0) # defines the pdf as a function by evaluating the expression g
  }
  maxf <- optimize(f,c(a,b), tol = 0.0001, maximum = TRUE) # uses the optimize function to find the maximum of the pdf
  samples <- rep(0,N) # creates a vector for storing the sample values
  i <- 0
  while ( i < N) { # while loop that runs until N samples are selected
    potsamp <- runif(1,a,b) # creating a potential sample
    testsamp <- runif(1, 0, maxf$objective + .1) # creating a test sample. Uses maxf$objective to get the max of the function, and adds .1 to account for errors in tolerance
    if ( testsamp < f(potsamp) ) { # tests if the potential sample should be rejected. If it is kept, it is kept in the samples vector
      samples[i+1] = potsamp
      i = i + 1
    }
  }
  samples # output the samples
}
