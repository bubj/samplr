samp2D <- function(fun,a,b,c,d,N) {
  g <- parse(text = fun)
  f <- function(z) {
    x <- z[1]
    y <- z[2]
    ifelse(a <= x & x <= b & c <= y & y <= d, eval(g[[1]]), 0)
  }
  maxf <- optim(c(b/2,d/2),f,control = list(fnscale = -1))
  maxf$value + .1
}
