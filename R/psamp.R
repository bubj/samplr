psamp <- function(fun, a, b, q, c = NULL, d = NULL) {
  if (is.null(c) & is.null(d)) {
    samples <- samplr(fun,a,b,10000)
    mean(samples < q)
  } else {
    samples <- samplr(fun, a = a, b = b, c = c, d = d, N = 10000)
    mean(samples[,1] + samples[,2] < q)
  }
}
