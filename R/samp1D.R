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

samp1D <- function(f,N) {
  text <- toString(body(f))
  leftboundpos <- rep(0,10)
  rightboundpos <- rep(0,10)

  # find the left bound of the function
  i <- 1
  pleftpos <- gregexpr('< x',text)
  for (j in 1:length(pleftpos[[1]])) {
    leftboundpos[i] <- pleftpos[[1]][j]
    i <- i + 1
  }
  pleftpos <- gregexpr('<= x',text)
  for (j in 1:length(pleftpos[[1]])) {
    leftboundpos[i] <- pleftpos[[1]][j]
    i <- i + 1
  }
  pleftpos <- gregexpr('x > ',text)
  for (j in 1:length(pleftpos[[1]])) {
    leftboundpos[i] <- pleftpos[[1]][j]
    i <- i + 1
  }
  pleftpos <- gregexpr('x >= ',text)
  for (j in 1:length(pleftpos[[1]])) {
    leftboundpos[i] <- pleftpos[[1]][j]
    i <- i + 1
  }
  leftboundpos <- leftboundpos[which(leftboundpos > 0)]
  leftbounds <- rep(0,length(leftboundpos))
  for (j in 1:length(leftboundpos)) {
    leftbounds[j] <- as.numeric(gsub("[^0-9\\-]","",(substr(text,leftboundpos[j]-5,leftboundpos[j]+5))))
  }
  leftbound <- min(leftbounds)
  rm(leftboundpos,leftbounds,pleftpos)

  #find the right bound of the function
  i <- 1
  prightpos <- gregexpr('> x',text)
  for (j in 1:length(prightpos[[1]])) {
    rightboundpos[i] <- prightpos[[1]][j]
    i <- i + 1
  }
  prightpos <- gregexpr('>= x',text)
  for (j in 1:length(prightpos[[1]])) {
    rightboundpos[i] <- prightpos[[1]][j]
    i <- i + 1
  }
  prightpos <- gregexpr('x < ',text)
  for (j in 1:length(prightpos[[1]])) {
    rightboundpos[i] <- prightpos[[1]][j]
    i <- i + 1
  }
  prightpos <- gregexpr('x <= ',text)
  for (j in 1:length(prightpos[[1]])) {
    rightboundpos[i] <- prightpos[[1]][j]
    i <- i + 1
  }
  rightboundpos <- rightboundpos[which(rightboundpos > 0)]
  rightbounds <- rep(0,length(rightboundpos))
  for (j in 1:length(rightboundpos)) {
    rightbounds[j] <- as.numeric(gsub("[^0-9\\-]","",(substr(text,rightboundpos[j]-5,rightboundpos[j]+5))))
  }
  rightbound <- max(rightbounds)
  rm(rightboundpos,rightbounds,prightpos)

  maxf <- optimize(f,c(leftbound,rightbound),maximum = TRUE)
  maxf <- maxf$objective + .1
  samples <- rep(0,N) # creates a vector for storing the sample values
  i <- 0
  while ( i < N) { # while loop that runs until N samples are selected
    potsamp <- runif(1,leftbound,rightbound) # creating a potential sample
    testsamp <- runif(1, 0, maxf) # creating a test sample. Uses maxf$objective to get the max of the function, and adds .1 to account for errors in tolerance
    if ( testsamp < f(potsamp) ) { # tests if the potential sample should be rejected. If it is kept, it is kept in the samples vector
      samples[i+1] = potsamp
      i = i + 1
    }
  }
  samples # output the samples
}
