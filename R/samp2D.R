#' Two Variable Rejection Sampler
#'
#' This function implements two variable rejection sampling for rvs with bounded support and which have bounded pdf.
#'
#' @param fun The pdf that we are sampling from, input as a string. The two variables sampled from must be input using the variables "x" and "y".
#'
#' @param a The lower bound of the variable "x" in the pdf.
#'
#' @param b The upper bound of the variable "x" in the pdf.
#'
#' @param c The lower bound of the variable "y" in the pdf.
#'
#' @param d The upper bound of the variable "y" in the pdf.
#'
#' @param N the number of samples output by the sampling function
#'
#' @return A matrix with N rows and 2 columns containing the samples from the pdf.
#' @export
#'
#' @examples
#'
#' samps <- samp2D("x + y", 0, 1, 0, 1, 10000)
#' samps <- data.frame(samps)
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'     geom_density_2d()
#'


samp2D <- function(fun,a,b,c,d,N) {
  g <- parse(text = fun)
  f <- function(z) {
    x <- z[1]
    y <- z[2]
    ifelse(a <= x & x <= b & c <= y & y <= d, eval(g[[1]]), 0)
  }
  maxf <- optim(c(b/2,d/2),f,control = list(fnscale = -1))
  samples <- matrix(rep(0,2*N), nrow = N, ncol = 2)
  i <- 0
  while( i < N) {
    psx <- runif(1,a,b)
    psy <- runif(1,c,d)
    testsamp <- runif(1,0, maxf$value + .1)
    if (testsamp < f(c(psx,psy))) {
      samples[i+1,1] <- psx
      samples[i+1,2] <- psy
      i <- i +1
    }
  }
  samples
}
