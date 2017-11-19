psamp <- function(fun, a, b, q, c = NULL, d = NULL) {
  if (is.null(c) & is.null(d)) {
    samples <- samp1D(fun, a, b, 10000)
    mean(samples < q)
  } else {
    samples <- samp2D(fun, a, b, c, d, 10000)
    mean(samples[,1] + samples[,2] < q)
  }
}
