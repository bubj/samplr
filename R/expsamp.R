expsamp <- function(f, g = NULL) {
  if (is.null(g)) {
    samples <- samp1D(f,10000)
    mean(samples)
  }
  else {
    samples <- samp1D(f,10000)
    values <- rep(0,10000)
    for (i in 1:10000) {
      values[i] <- g(c(samples[i,1],samples[i,2]))
    }
    mean(values)
  }
}
